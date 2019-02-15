from arithmetic import *
from case import Case, Const
import ipopt
import numpy as np
from scipy.sparse import *
from scipy.optimize import *
from pdb import *

class opf_mdl(object):
    
    def __init__(self, c):
        self.case = c

    def objective(self, x):
        return costfcn(x, self.case)
        
    def gradient(self, x):
        return costfcn_jac(x, self.case)

    def constraints(self, x):       
        const = Const() 
        ii = get_var_idx(self.case)
        tload = sum(self.case.bus[:,const.PD]) / self.case.mva_base
        return sum(x[ii['i1']['pg']:ii['iN']['pg']]) - tload

    def jacobian(self, x):        
        const = Const() 
        ii = get_var_idx(self.case)
        simple_powerbalance = np.zeros_like(x)
        simple_powerbalance[ii['i1']['pg']:ii['iN']['pg']] = 1
        return simple_powerbalance

    def intermediate(
        self,
        alg_mod,
        iter_count,
        obj_value,
        inf_pr,
        inf_du,
        mu,
        d_norm,
        regularization_size,
        alpha_du,
        alpha_pr,
        ls_trials
        ):
        
        #
        # Example for the use of the intermediate callback.
        #
        print "Objective value at iteration #%d is - %g" % (iter_count, obj_value)


def runcopf(c, flat_start):
    const = Const()
    
    nb     = c.bus.shape[0]
    ng     = c.gen.shape[0]
    nbr    = c.branch.shape[0]
    neq    = 2 * nb
    niq    = 2 * ng + nb + nbr
    neqnln = 2 * nb
    niqnln = nbr

    ii = get_var_idx(c)

    if flat_start:
        x0 = np.concatenate((deg2rad(c.bus.take(const.VA, axis=1)), \
            c.bus.take([const.VMAX, const.VMIN], axis=1).mean(axis=1), \
            c.gen.take([const.PMAX, const.PMIN], axis=1).mean(axis=1) / c.mva_base, \
            c.gen.take([const.QMAX, const.QMIN], axis=1).mean(axis=1) / c.mva_base), axis=0)
    else:
        x0   = np.genfromtxt(c.path+"x0.csv", delimiter=',')
        
    xmin = np.concatenate((-np.inf * np.ones(nb), \
                           c.bus.take(const.VMIN, axis=1), \
                           c.gen.take(const.PMIN, axis=1) / c.mva_base, \
                           c.gen.take(const.QMIN, axis=1) / c.mva_base), axis=0)
    xmax = np.concatenate((np.inf * np.ones(nb), \
                           c.bus.take(const.VMAX, axis=1), \
                           c.gen.take(const.PMAX, axis=1) / c.mva_base, \
                           c.gen.take(const.QMAX, axis=1) / c.mva_base), axis=0)

    xmin[(c.bus.take(const.BUS_TYPE, axis=1) == 3).nonzero()] = 0
    xmax[(c.bus.take(const.BUS_TYPE, axis=1) == 3).nonzero()] = 0

    ####################################################################
    # Test Environment
    #################################################################### 
    cl = [0.]
    cu = [0.]

    nlp = ipopt.problem(n=len(x0), m=len(cl), lb=xmin, ub=xmax, cl=cl, cu=cu, \
        problem_obj=opf_mdl(c))

    res = nlp.solve(x0)

    # ii = get_var_idx(c)
    # res_va = rad2deg(res.x[ii['i1']['va']:ii['iN']['va']])
    # res_vm = res.x[ii['i1']['vm']:ii['iN']['vm']]
    # res_pg = res.x[ii['i1']['pg']:ii['iN']['pg']] * c.mva_base
    # res_qg = res.x[ii['i1']['qg']:ii['iN']['qg']] * c.mva_base

    # float_fmtr = {'float_kind': lambda x: "%7.3f" % x}

    # print('___________')  
    # # print('     Statue | Exit mode %d' % res.status)
    # # print('    Message | %s' % res.message)
    # # print('       Iter | %d' % res.nit)
    # print('  Objective | %10.3f $/hr' % res.fun)
    # print('  VA (deg)  | %s' % np.array2string(res_va[0:7], formatter=float_fmtr))
    # print('  VM (pu)   | %s' % np.array2string(res_vm[0:7], formatter=float_fmtr))
    # print('  PG (MW)   | %s' % np.array2string(res_pg, formatter=float_fmtr))
    # print('  QG (MVAR) | %s' % np.array2string(res_qg, formatter=float_fmtr))
    # print('___________ | ')  
    
    set_trace()
    

# region [ Cost-Related Functions ]

def costfcn(x, c):
    ng = c.gen.shape[0]
    ii = get_var_idx(c)

    pg = c.mva_base * x[ii['i1']['pg']:ii['iN']['pg']]
    qg = c.mva_base * x[ii['i1']['qg']:ii['iN']['qg']]

    gencost = np.zeros(ng)
    for gi in range(ng):
        gencost[gi] = polycost(c.gencost[gi], pg[gi])

    return gencost.sum()

def costfcn_jac(x, c):
    ng = c.gen.shape[0]
    ii = get_var_idx(c)

    pg = c.mva_base * x[ii['i1']['pg']:ii['iN']['pg']]
    qg = c.mva_base * x[ii['i1']['qg']:ii['iN']['qg']]

    gencost_jac = np.zeros(ng)
    for gi in range(ng):
        gencost_jac[gi] = c.mva_base * polycost_jac(c.gencost[gi], pg[gi])

    df = np.zeros(x.shape[0])
    df[ii['i1']['pg']:ii['iN']['pg']] = gencost_jac

    return df

def costfcn_hess(x, c):
    ng = c.gen.shape[0]
    ii = get_var_idx(c)

    pg = c.mva_base * x[ii['i1']['pg']:ii['iN']['pg']]
    qg = c.mva_base * x[ii['i1']['qg']:ii['iN']['qg']]

    gencost_hess = np.zeros(ng)
    for gi in range(ng):
        gencost_hess[gi] = c.mva_base ** 2 * polycost_hess(c.gencost[gi], pg[gi])

    d2f = np.zeros((x.shape[0], x.shape[0]))
    d2f[ii['i1']['pg']:ii['iN']['pg'], ii['i1']['pg']:ii['iN']['pg']] = np.diag(gencost_hess)

    return d2f

def polycost(cost_metrics, pg):
    const = Const()
    cost = 0.
    pn = int(cost_metrics[const.NCOST])
    for pi in range(pn):
        cost += cost_metrics[-(1+pi)] * pg ** pi
    
    return cost

def polycost_jac(cost_metrics, pg):
    const = Const()
    cost = 0.
    pn = int(cost_metrics[const.NCOST])
    for pi in range(1, pn):
        cost += pi * cost_metrics[-(1+pi)] * pg ** (pi - 1)

    return cost

def polycost_hess(cost_metrics, pg):
    const = Const()
    cost = 0.
    pn = int(cost_metrics[const.NCOST])
    for pi in range(2, pn):
        cost += pi * cost_metrics[-(1+pi)] * pg ** (pi - 2)

    return cost

# endregion


# region [ Constraint Functions ]

def build_bound_cons(xmin, xmax):
    b = ()
    for vi in range(len(xmin)):
        b += ((xmin[vi], xmax[vi]),)
    return b
    
def acpf_consfcn(x, c):
    const = Const()

    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    nbr = c.branch.shape[0]

    ii = get_var_idx(c)
    va = x[ii['i1']['va']:ii['iN']['va']]
    vm = x[ii['i1']['vm']:ii['iN']['vm']]
    pg = x[ii['i1']['pg']:ii['iN']['pg']]
    qg = x[ii['i1']['qg']:ii['iN']['qg']]

    vcplx = vm * np.exp(1j * va)
    c.gen[:, const.PG] = c.mva_base * pg
    c.gen[:, const.QG] = c.mva_base * qg

    Ybus, _, _ = makeYbus(c)
    Sbus = makeSbus(c.mva_base, c.bus, c.gen)
    mis = - Sbus + \
          vcplx * np.asarray(np.conj(Ybus * np.matrix(vcplx).T)).flatten() 

    return np.concatenate((np.real(mis), np.imag(mis)))

# def acpf_consfcn_jac(x, c):
#     const = Const()

#     nb = c.bus.shape[0]
#     ng = c.gen.shape[0]
#     nbr = c.branch.shape[0]
#     nx = 2 * (nb + ng)

#     ii = get_var_idx(c)
#     va = x[ii['i1']['va']:ii['iN']['va']]
#     vm = x[ii['i1']['vm']:ii['iN']['vm']]
#     pg = x[ii['i1']['pg']:ii['iN']['pg']]
#     qg = x[ii['i1']['qg']:ii['iN']['qg']]

#     vcplx = vm * np.exp(1j * va)
#     c.gen[:, const.PG] = c.mva_base * pg
#     c.gen[:, const.QG] = c.mva_base * qg

#     Ybus = makeYbus(c)
#     Sbus = makeSbus(c.mva_base, c.bus, c.gen)

#     dSdVa = dSbus_dVa(Ybus, vcplx)
#     dSdVm = dSbus_dVm(Ybus, vcplx)

#     # dg = csr_matrix(, shape=(2*nb, nx))

#     return 0

# def dSbus_dVa(Ybus, V):

#     return 0

# def dSbus_dVm(Ybus, V):

#     return 0


def linerating_consfcn(x, c):
    const = Const()

    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    nbr = c.branch.shape[0]

    ii = get_var_idx(c)
    va = x[ii['i1']['va']:ii['iN']['va']]
    vm = x[ii['i1']['vm']:ii['iN']['vm']]

    vcplx = vm * np.exp(1j * va)

    _, Yf, Yt = makeYbus(c)

    fbus_idx = np.array(c.branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(c.branch[:, const.T_BUS] - 1, dtype=int)

    flow_max = (c.branchrate / c.mva_base ) ** 2
    Sf = vcplx[fbus_idx] * np.conj(Yf * vcplx)
    St = vcplx[tbus_idx] * np.conj(Yf * vcplx)
    
    Sfreal_sq = np.real(Sf) ** 2
    Sfimag_sq = np.imag(Sf) ** 2
    Streal_sq = np.real(St) ** 2
    Stimag_sq = np.imag(St) ** 2

    return np.concatenate((flow_max - Sfreal_sq - Sfimag_sq, \
                           flow_max - Streal_sq - Stimag_sq))

# endregion


# region [ Powerflow-related Functions ]

def makeSbus(mva_base, bus, gen):
    const = Const()

    nb = bus.shape[0]
    ng = gen.shape[0]

    g_idx = np.array(range(ng), dtype=int)
    g_busnum = np.array(gen[:, const.GEN_BUS] - 1, dtype=int)

    Sbusg = np.ones(nb) * (0 + 0j)
    Sbusg[g_busnum] = (gen[:, const.PG] + 1j * gen[:, const.QG]) / mva_base
    Sbusd = (bus[:, const.PD] + 1j * bus[:, const.QD]) / mva_base
    
    return Sbusg - Sbusd

def makeYbus(c):
    const = Const()

    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    nbr = c.branch.shape[0]

    Ys = 1 / (c.branch[:, const.BR_R] + 1j * c.branch[:, const.BR_X])
    Bc = c.branch[:, const.BR_B]
    tap = np.ones(nbr)
    tap_idx = (c.branch.take(const.TAP, axis=1) != 0).nonzero()
    tap[tap_idx] = c.branch[tap_idx, const.TAP]

    Ytt = Ys + 1j * Bc/2
    Yff = Ytt / (tap * np.conj(tap))
    Yft = - Ys / np.conj(tap)
    Ytf = - Ys / tap

    ysh = (c.bus[:, const.GS] + 1j * c.bus[:, const.BS]) / c.mva_base

    b_idx = np.array(range(nb), dtype=int)
    br_idx = np.array(range(nbr), dtype=int)
    fbus_idx = np.array(c.branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(c.branch[:, const.T_BUS] - 1, dtype=int)

    Cf = csr_matrix((np.ones(nbr), (br_idx, fbus_idx)), shape=(nbr,nb))
    Ct = csr_matrix((np.ones(nbr), (br_idx, tbus_idx)), shape=(nbr,nb))


    Yf = csr_matrix((np.concatenate((Yff, Yft)), \
                     (np.concatenate((br_idx, br_idx)), np.concatenate((fbus_idx, tbus_idx)))), \
                    shape=(nbr,nb))
    Yt = csr_matrix((np.concatenate((Ytf, Ytt)), \
                     (np.concatenate((br_idx, br_idx)), np.concatenate((fbus_idx, tbus_idx)))), \
                    shape=(nbr,nb))
    Ysh = csr_matrix((ysh, (b_idx, b_idx)), shape=(nb,nb))

    Ybus = Cf.T * Yf + Ct.T * Yt + Ysh
    
    return Ybus, Yf, Yt

# endregion

def get_var_idx(c):
    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    ii = {'N': {'va': nb, 'vm': nb, 'pg': ng, 'qg': ng}, \
          'i1': {'va': 0, 'vm': nb, 'pg': 2*nb, 'qg': 2*nb+ng}, \
          'iN': {'va': nb, 'vm': 2*nb, 'pg': 2*nb+ng, 'qg': 2*(nb+ng)}}
    return ii