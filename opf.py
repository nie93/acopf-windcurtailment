from case import Case, Const
import numpy as np
from scipy.optimize import *
# from scipy.optimize import Bounds, LinearConstraint, NonlinearConstraint
import pdb


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
    
    f = costfcn(x0,c)

    bounds = Bounds(xmin, xmax)

    simple_powerbalance = np.zeros(x0.shape[0])
    tload = sum(c.bus[:,const.PD]) / c.mva_base
    simple_powerbalance[ii['i1']['pg']:ii['iN']['pg']] = 1
    simple_lincons = LinearConstraint(simple_powerbalance, tload, tload)

    pdb.set_trace()

    f_fcn = lambda x: costfcn(x, c)

    res = minimize(f_fcn, x0, constraints=[simple_lincons], bounds=bounds)

    print(tload)
    print(res.x[ii['i1']['pg']:ii['iN']['pg']])
    print(sum(res.x[ii['i1']['pg']:ii['iN']['pg']]))


    

def costfcn(x, c):
    ng = c.gen.shape[0]
    ii = get_var_idx(c)

    pg = c.mva_base * x[ii['i1']['pg']:ii['iN']['pg']]
    qg = c.mva_base * x[ii['i1']['qg']:ii['iN']['qg']]

    gencost = np.zeros(ng)
    for gi in range(ng):
        gencost[gi] = polycost(c.gencost[gi], pg[gi])
    f = gencost.sum()
    return f

def polycost(cost_metrics, pg):
    const = Const()
    cost = 0.
    pn = int(cost_metrics[const.NCOST])
    for pi in range(pn):
        cost += cost_metrics[-(1+pi)] * pg ** pi
    return cost

def get_var_idx(c):
    nb = c.bus.shape[0]
    ng = c.gen.shape[0]
    ii = {'N': {'va': nb, 'vm': nb, 'pg': ng, 'qg': ng}, \
          'i1': {'va': 0, 'vm': nb, 'pg': 2*nb, 'qg': 2*nb+ng}, \
          'iN': {'va': nb, 'vm': 2*nb, 'pg': 2*nb+ng, 'qg': 2*(nb+ng)}}
    return ii
