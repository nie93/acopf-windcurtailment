from case import Case, Const
import numpy as np
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

    x0 = np.concatenate((np.zeros(nb), \
                         c.bus.take([const.VMAX, const.VMIN], axis=1).mean(axis=1), \
                         c.gen.take([const.PMAX, const.PMIN], axis=1).mean(axis=1) / c.mva_base, \
                         c.gen.take([const.QMAX, const.QMIN], axis=1).mean(axis=1) / c.mva_base), axis=0)

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
    # eqcons = NonlinearConstraint()
    # ineqcons = NonlinearConstraint()
    # eqcons = {'type': 'eq',
    #           'fun' : lambda x: eqconsfcn(x, c),
    #           'jac' : lambda x: eqconsfcn_jac(x,c)}
              
    # ineqcons = {'type': 'eq',
    #             'fun' : lambda x: ineqconsfcn(x, c),
    #             'jac' : lambda x: ineqconsfcn_jac(x,c)}


    res = minimize(f_fcn, x0, jac=df_fcn, hess=d2f_fcn, \
                   constraints=[simple_lincons], bounds=Bounds(xmin, xmax))

    print(res)
    

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

def get_var_idx(c):
    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    ii = {'N': {'va': nb, 'vm': nb, 'pg': ng, 'qg': ng}, \
          'i1': {'va': 0, 'vm': nb, 'pg': 2*nb, 'qg': 2*nb+ng}, \
          'iN': {'va': nb, 'vm': 2*nb, 'pg': 2*nb+ng, 'qg': 2*(nb+ng)}}
    return ii
