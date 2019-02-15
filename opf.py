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
        return np.concatenate((acpf_consfcn(x, self.case), linerating_consfcn(x, self.case)))

    def jacobian(self, x):
        return np.concatenate((acpf_consfcn_jac(x, self.case), linerating_consfcn_jac(x, self.case)))

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
        
        print('Objective value at iteration #%d is - %g' % (iter_count, obj_value))


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
        # x0 = np.concatenate((deg2rad(c.bus.take(const.VA, axis=1)), \
        #     c.bus.take([const.VMAX, const.VMIN], axis=1).mean(axis=1), \
        #     c.gen.take([const.PMAX, const.PMIN], axis=1).mean(axis=1) / c.mva_base, \
        #     c.gen.take([const.QMAX, const.QMIN], axis=1).mean(axis=1) / c.mva_base), axis=0)
        x0 = np.concatenate((np.zeros(nb), \
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
    cl = np.concatenate((np.zeros(2 * nb), -np.inf * np.ones(2 * nbr)))
    cu = np.concatenate((np.zeros(2 * nb), np.zeros(2 * nbr)))

    nlp = ipopt.problem(n=len(x0), m=len(cl), lb=xmin, ub=xmax, cl=cl, cu=cu, \
        problem_obj=opf_mdl(c))

    res_x, res_info = nlp.solve(x0)

    ii = get_var_idx(c)
    res_va = rad2deg(res_x[ii['i1']['va']:ii['iN']['va']])
    res_vm = res_x[ii['i1']['vm']:ii['iN']['vm']]
    res_pg = res_x[ii['i1']['pg']:ii['iN']['pg']] * c.mva_base
    res_qg = res_x[ii['i1']['qg']:ii['iN']['qg']] * c.mva_base

    float_fmtr = {'float_kind': lambda x: "%7.3f" % x}

    print('___________')  
    print('     Status | Exit mode %d' % res_info['status'])
    print('    Message | %s' % res_info['status_msg'])
    # print('       Iter | %d' % res_info['iter_count'])
    print('  Objective | %.3f $/hr' % res_info['obj_val'])
    print('  VA (deg)  | %s' % np.array2string(res_va[0:7], formatter=float_fmtr))
    print('  VM (pu)   | %s' % np.array2string(res_vm[0:7], formatter=float_fmtr))
    print('  PG (MW)   | %s' % np.array2string(res_pg, formatter=float_fmtr))
    print('  QG (MVAR) | %s' % np.array2string(res_qg, formatter=float_fmtr))
    print('___________ | ')  
    
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

def acpf_consfcn_jac(x, c):
    const = Const()

    nb  = c.bus.shape[0]
    ng  = c.gen.shape[0]
    nbr = c.branch.shape[0]
    nx  = 2 * (nb + ng)

    ii       = get_var_idx(c)
    g_idx    = np.array(range(0, ng), dtype=int)
    gbus_idx = np.array(c.gen[:, const.GEN_BUS] - 1, dtype=int)
    cons_idx = np.array(range(2 * nb), dtype=int)
    va_idx   = np.array(range(ii['i1']['va'], ii['iN']['va']), dtype=int)
    vm_idx   = np.array(range(ii['i1']['vm'], ii['iN']['vm']), dtype=int)
    pg_idx   = np.array(range(ii['i1']['pg'], ii['iN']['pg']), dtype=int)
    qg_idx   = np.array(range(ii['i1']['qg'], ii['iN']['qg']), dtype=int)
    x_idx    = np.concatenate((va_idx, vm_idx, pg_idx, qg_idx))

    va = x[va_idx]
    vm = x[vm_idx]
    pg = x[pg_idx]
    qg = x[qg_idx]

    vcplx = vm * np.exp(1j * va)
    c.gen[:, const.PG] = c.mva_base * pg
    c.gen[:, const.QG] = c.mva_base * qg

    Ybus, _, _ = makeYbus(c)
    Sbus = makeSbus(c.mva_base, c.bus, c.gen)

    dSdVa, dSdVm = dSbus_dV(Ybus, vcplx)
    dSdV = hstack((dSdVa, dSdVm))
    neg_Cg = csr_matrix((-np.ones(ng), (gbus_idx, g_idx)), shape=(nb, ng))
    zeros_b_g = csr_matrix(([],([],[])), shape=(nb, ng))

    dpinj = hstack((csr_matrix(np.real(dSdV.toarray())), neg_Cg, zeros_b_g))
    dqinj = hstack((csr_matrix(np.imag(dSdV.toarray())), zeros_b_g, neg_Cg))

    dg = lil_matrix((2*nb, nx), dtype=float)
    dg[:, x_idx] = vstack((dpinj, dqinj))

    return dg.toarray().reshape(-1)

def dSbus_dV(Ybus, V):
    nb = V.shape[0]
    I = Ybus.dot(V)

    b_idx = np.array(range(nb), dtype=int)

    diagV = csr_matrix((V, (b_idx, b_idx)), shape=(nb, nb))
    diagI = csr_matrix((I, (b_idx, b_idx)), shape=(nb, nb))
    diagVnorm = csr_matrix((V/np.abs(V), (b_idx, b_idx)), shape=(nb, nb))

    dSbus_dVa = ((1j) * diagV).dot(np.conj(diagI - Ybus.dot(diagV)))
    dSbus_dVm = diagV.dot(np.conj(Ybus.dot(diagVnorm))) + np.conj(diagI).dot(diagVnorm)

    return dSbus_dVa, dSbus_dVm 

def linerating_consfcn(x, c):
    const = Const()

    nb  = c.bus.shape[0]
    ng  = c.gen.shape[0]
    nbr = c.branch.shape[0]

    ii       = get_var_idx(c)
    va_idx   = np.array(range(ii['i1']['va'], ii['iN']['va']), dtype=int)
    vm_idx   = np.array(range(ii['i1']['vm'], ii['iN']['vm']), dtype=int)
    fbus_idx = np.array(c.branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(c.branch[:, const.T_BUS] - 1, dtype=int)
    x_idx    = np.concatenate((va_idx, vm_idx))

    va    = x[va_idx]
    vm    = x[vm_idx]
    vcplx = vm * np.exp(1j * va)

    _, Yf, Yt = makeYbus(c)

    flow_max = (c.branchrate / c.mva_base ) ** 2
    Sf = vcplx[fbus_idx] * np.conj(Yf * vcplx)
    St = vcplx[tbus_idx] * np.conj(Yf * vcplx)
    
    Sfreal_sq = np.real(Sf) ** 2
    Sfimag_sq = np.imag(Sf) ** 2
    Streal_sq = np.real(St) ** 2
    Stimag_sq = np.imag(St) ** 2

    return np.concatenate((Sfreal_sq + Sfimag_sq - flow_max, \
                           Streal_sq + Stimag_sq - flow_max))

def linerating_consfcn_jac(x, c):
    const = Const()

    nb  = c.bus.shape[0]
    ng  = c.gen.shape[0]
    nbr = c.branch.shape[0]
    nx  = 2 * (nb + ng)

    ii       = get_var_idx(c)
    va_idx   = np.array(range(ii['i1']['va'], ii['iN']['va']), dtype=int)
    vm_idx   = np.array(range(ii['i1']['vm'], ii['iN']['vm']), dtype=int)
    fbus_idx = np.array(c.branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(c.branch[:, const.T_BUS] - 1, dtype=int)
    x_idx    = np.concatenate((va_idx, vm_idx))

    va    = x[va_idx]
    vm    = x[vm_idx]
    vcplx = vm * np.exp(1j * va)

    _, Yf, Yt = makeYbus(c)

    br_idx   = np.array(range(0, nbr), dtype=int)
    fbus_idx = np.array(c.branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(c.branch[:, const.T_BUS] - 1, dtype=int)

    dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft = dSbr_dV(c.branch, Yf, Yt, vcplx)
    df_dVa, df_dVm, dt_dVa, dt_dVm = dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft)

    dh = lil_matrix((2*nbr, nx), dtype=float)
    dh[:, x_idx] = vstack((hstack((df_dVa, df_dVm)),hstack((dt_dVa, dt_dVm))))

    return dh.toarray().reshape(-1)

def dSbr_dV(branch, Yf, Yt, V):
    const = Const()

    nb  = V.shape[0]
    nbr = branch.shape[0]

    b_idx    = np.array(range(nb), dtype=int)
    br_idx   = np.array(range(nbr), dtype=int)
    fbus_idx = np.array(branch[:, const.F_BUS] - 1, dtype=int)
    tbus_idx = np.array(branch[:, const.T_BUS] - 1, dtype=int)

    Vnorm = V / np.abs(V)
    If    = Yf.dot(V)
    It    = Yt.dot(V)
    Sf    = V[fbus_idx] * np.conj(If)
    St    = V[tbus_idx] * np.conj(It)

    diagVf    = csr_matrix((V[fbus_idx], (br_idx, br_idx)), shape=(nbr, nbr))
    diagIf    = csr_matrix((If, (br_idx, br_idx)), shape=(nbr, nbr))
    diagVt    = csr_matrix((V[tbus_idx], (br_idx, br_idx)), shape=(nbr, nbr))
    diagIt    = csr_matrix((It, (br_idx, br_idx)), shape=(nbr, nbr))
    diagV     = csr_matrix((V, (b_idx, b_idx)), shape=(nb, nb))
    diagVnorm = csr_matrix((Vnorm, (b_idx, b_idx)), shape=(nb, nb))

    Cvf     = csr_matrix((V[fbus_idx], (br_idx, fbus_idx)), shape=(nbr, nb))
    Cvt     = csr_matrix((V[tbus_idx], (br_idx, tbus_idx)), shape=(nbr, nb))
    Cvnormf = csr_matrix((Vnorm[fbus_idx], (br_idx, fbus_idx)), shape=(nbr, nb))
    Cvnormt = csr_matrix((Vnorm[tbus_idx], (br_idx, tbus_idx)), shape=(nbr, nb))

    dSf_dVa = (1j) * ( np.conj(diagIf) * Cvf - diagVf * np.conj(Yf * diagV) )
    dSt_dVa = (1j) * ( np.conj(diagIt) * Cvt - diagVt * np.conj(Yt * diagV) ) 
    dSf_dVm = diagVf * np.conj(Yf * diagVnorm) + np.conj(diagIf) * Cvnormf
    dSt_dVm = diagVt * np.conj(Yt * diagVnorm) + np.conj(diagIt) * Cvnormt

    return dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St

def dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St):
    nbr = Sf.shape[0]
    br_idx = np.array(range(nbr), dtype=int)

    dAf_dPf = csr_matrix((2 * np.real(Sf), (br_idx, br_idx)), shape=(nbr, nbr))
    dAf_dQf = csr_matrix((2 * np.imag(Sf), (br_idx, br_idx)), shape=(nbr, nbr))
    dAt_dPt = csr_matrix((2 * np.real(St), (br_idx, br_idx)), shape=(nbr, nbr))
    dAt_dQt = csr_matrix((2 * np.imag(St), (br_idx, br_idx)), shape=(nbr, nbr))

    dAf_dVm = dAf_dPf * np.real(dSf_dVm) + dAf_dQf * np.imag(dSf_dVm)
    dAf_dVa = dAf_dPf * np.real(dSf_dVa) + dAf_dQf * np.imag(dSf_dVa)
    dAt_dVm = dAt_dPt * np.real(dSt_dVm) + dAt_dQt * np.imag(dSt_dVm)
    dAt_dVa = dAt_dPt * np.real(dSt_dVa) + dAt_dQt * np.imag(dSt_dVa)

    return dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm

def simple_powerbalance_consfcn(x, c):
    const = Const() 
    ii = get_var_idx(c)
    tload = sum(c.bus[:,const.PD]) / c.mva_base
    return sum(x[ii['i1']['pg']:ii['iN']['pg']]) - tload

def simple_powerbalance_consfcn_jac(x, c):
    const = Const() 
    ii = get_var_idx(c)
    jac = np.zeros_like(x)
    jac[ii['i1']['pg']:ii['iN']['pg']] = 1
    return jac

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