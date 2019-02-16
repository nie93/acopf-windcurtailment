from case import Case, Const
import numpy as np
from scipy.sparse import *
from scipy.optimize import *
from pdb import *


def runcopf(c):
    const = Const()
    
    nb     = c.bus.shape[0]
    ng     = c.gen.shape[0]
    nbr    = c.branch.shape[0]
    neq    = 2 * nb
    niq    = 2 * ng + nb + nbr
    neqnln = 2*nb
    niqnln = nbr

    ii = get_var_idx(c)

    # x0 = np.concatenate((deg2rad(c.bus.take(const.VA, axis=1)), \
    #     c.bus[:,[const.VMIN, const.VMAX]].mean(axis=1), \
    #     c.gen[:,[const.PMAX, const.PMIN]].mean(axis=1) / c.mva_base, \
    #     c.gen[:,[const.QMAX, const.QMIN]].mean(axis=1) / c.mva_base), axis=0)
    x0 = np.concatenate((np.zeros(nb), \
        c.bus[:, [const.VMIN, const.VMAX]].mean(axis=1), \
        c.gen[:, [const.PMAX, const.PMIN]].mean(axis=1) / c.mva_base, \
        c.gen[:, [const.QMAX, const.QMIN]].mean(axis=1) / c.mva_base), axis=0)    
    
    xmin = np.concatenate((-np.inf * np.ones(nb), \
                           c.bus[:, const.VMIN], \
                           c.gen[:, const.PMIN] / c.mva_base, \
                           c.gen[:, const.QMIN] / c.mva_base), axis=0)
    xmax = np.concatenate((np.inf * np.ones(nb), \
                           c.bus[:, const.VMAX], \
                           c.gen[:, const.PMAX] / c.mva_base, \
                           c.gen[:, const.QMAX] / c.mva_base), axis=0)

    xmin[(c.bus[:, const.BUS_TYPE] == 3).nonzero()] = 0
    xmax[(c.bus[:, const.BUS_TYPE] == 3).nonzero()] = 0
   
    ####################################################################
    # Polynomial Cost Functions (f)
    #################################################################### 
    f_fcn   = lambda x: costfcn(x, c)
    df_fcn  = lambda x: costfcn_jac(x, c)
    d2f_fcn = lambda x: costfcn_hess(x, c)

    ####################################################################
    # Simple Power Balance Constraint (Linear, lossless)
    ####################################################################   
    simple_powerbalance = np.zeros(x0.shape[0])
    tload = sum(c.bus[:,const.PD]) / c.mva_base
    simple_powerbalance[ii['i1']['pg']:ii['iN']['pg']] = 1
    simple_lincons = LinearConstraint(simple_powerbalance, tload, tload)
    # simple_lincons = {'type': 'eq',
    #                   'fun': lambda x : sum(x[ii['i1']['pg']:ii['iN']['pg']]) - tload}
    
    ####################################################################
    # Nonlinear Power Flow Constraints (g: eqcons, h: ineqcons)
    ####################################################################   
    g_fcn   = lambda x: acpf_consfcn(x, c)
    # dg_fcn  = lambda x: acpf_consfcn_jac(x, c)
    # d2g_fcn = lambda x: acpf_consfcn_hess(x, c)

    h_fcn   = lambda x: linerating_consfcn(x, c)
    # dh_fcn  = lambda x: linerating_consfcn_jac(x, c)
    # d2h_fcn = lambda x: linerating_consfcn_hess(x, c)

    eqcons = NonlinearConstraint(g_fcn, 0, 0)
    ineqcons = NonlinearConstraint(h_fcn, -np.inf, 0)

    ####################################################################
    # Test Environment
    ####################################################################    
    res = minimize(f_fcn, x0, jac=df_fcn, hess=d2f_fcn, method='trust-constr', \
        constraints=(eqcons, ineqcons), bounds=Bounds(xmin, xmax), \
        options={'disp': True})

    ii = get_var_idx(c)
    res_va = rad2deg(res.x[ii['i1']['va']:ii['iN']['va']])
    res_vm = res.x[ii['i1']['vm']:ii['iN']['vm']]
    res_pg = res.x[ii['i1']['pg']:ii['iN']['pg']] * c.mva_base
    res_qg = res.x[ii['i1']['qg']:ii['iN']['qg']] * c.mva_base

    float_fmtr = {'float_kind': lambda x: "%7.3f" % x}

    print('___________')  
    print('     Statue | Exit mode %d' % res.status)
    print('    Message | %s' % res.message)
    print('       Iter | %d' % res.niter)
    print('  Objective | %.3f $/hr' % res.fun)
    print('  VA (deg)  | %s' % np.array2string(res_va[0:7], formatter=float_fmtr))
    print('  VM (pu)   | %s' % np.array2string(res_vm[0:7], formatter=float_fmtr))
    print('  PG (MW)   | %s' % np.array2string(res_pg, formatter=float_fmtr))
    print('  QG (MVAR) | %s' % np.array2string(res_qg, formatter=float_fmtr))
    print('___________ | ') 

    # set_trace()
    

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

#     # dg = bsr_matrix(, shape=(2*nb, nx))

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

    return np.concatenate((Sfreal_sq + Sfimag_sq - flow_max, Streal_sq + Stimag_sq - flow_max))

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

    Cf = bsr_matrix((np.ones(nbr), (br_idx, fbus_idx)), shape=(nbr,nb))
    Ct = bsr_matrix((np.ones(nbr), (br_idx, tbus_idx)), shape=(nbr,nb))


    Yf = bsr_matrix((np.concatenate((Yff, Yft)), \
                     (np.concatenate((br_idx, br_idx)), np.concatenate((fbus_idx, tbus_idx)))), \
                    shape=(nbr,nb))
    Yt = bsr_matrix((np.concatenate((Ytf, Ytt)), \
                     (np.concatenate((br_idx, br_idx)), np.concatenate((fbus_idx, tbus_idx)))), \
                    shape=(nbr,nb))
    Ysh = bsr_matrix((ysh, (b_idx, b_idx)), shape=(nb,nb))

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

def deg2rad(d):
    return d / 180 * np.pi

def rad2deg(r):
    return r / np.pi * 180