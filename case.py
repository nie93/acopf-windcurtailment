import csv

class Const(object):
    def __init__(self):
        self.BUS_I       = 1    ## bus number (1 to 29997)
        self.BUS_TYPE    = 2    ## bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
        self.PD          = 3    ## Pd, real power demand (MW)
        self.QD          = 4    ## Qd, reactive power demand (MVAr)
        self.GS          = 5    ## Gs, shunt conductance (MW at V = 1.0 p.u.)
        self.BS          = 6    ## Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
        self.BUS_AREA    = 7    ## area number, 1-100
        self.VM          = 8    ## Vm, voltage magnitude (p.u.)
        self.VA          = 9    ## Va, voltage angle (degrees)
        self.BASE_KV     = 10   ## baseKV, base voltage (kV)
        self.ZONE        = 11   ## zone, loss zone (1-999)
        self.VMAX        = 12   ## maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
        self.VMIN        = 13   ## minVm, minimum voltage magnitude (p.u.)      (not in PTI format)
        self.F_BUS       = 1    ## f, from bus number
        self.T_BUS       = 2    ## t, to bus number
        self.BR_R        = 3    ## r, resistance (p.u.)
        self.BR_X        = 4    ## x, reactance (p.u.)
        self.BR_B        = 5    ## b, total line charging susceptance (p.u.)
        self.RATE_A      = 6    ## rateA, MVA rating A (long term rating)
        self.RATE_B      = 7    ## rateB, MVA rating B (short term rating)
        self.RATE_C      = 8    ## rateC, MVA rating C (emergency rating)
        self.TAP         = 9    ## ratio, transformer off nominal turns ratio
        self.SHIFT       = 10   ## angle, transformer phase shift angle (degrees)
        self.BR_STATUS   = 11   ## initial branch status, 1 - in service, 0 - out of service
        self.ANGMIN      = 12   ## minimum angle difference, angle(Vf) - angle(Vt) (degrees)
        self.ANGMAX      = 13   ## maximum angle difference, angle(Vf) - angle(Vt) (degrees)
        self.PF          = 14   ## real power injected at "from" bus end (MW)       (not in PTI format)
        self.QF          = 15   ## reactive power injected at "from" bus end (MVAr) (not in PTI format)
        self.PT          = 16   ## real power injected at "to" bus end (MW)         (not in PTI format)
        self.QT          = 17   ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)
        self.MODEL       = 1    ## cost model, 1 = piecewise linear, 2 = polynomial 
        self.STARTUP     = 2    ## startup cost in US dollars
        self.SHUTDOWN    = 3    ## shutdown cost in US dollars
        self.NCOST       = 4    ## number breakpoints in piecewise linear cost function,
                                ## or number of coefficients in polynomial cost function
        self.COST        = 5    ## parameters defining total cost function begin in this col
                                ## (MODEL = 1) : p0, f0, p1, f1, ..., pn, fn
                                ##      where p0 < p1 < ... < pn and the cost f(p) is defined
                                ##      by the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of
                                ##      the end/break-points of the piecewise linear cost
                                ## (MODEL = 2) : cn, ..., c1, c0
                                ##      n+1 coefficients of an n-th order polynomial cost fcn,
                                ##      starting with highest order, where cost is
                                ##      f(p) = cn*p^n + ... + c1*p + c0
        

class Case(object):

    def __init__(self):
        self.mva_base = 100

    def import_case(self, path):
        self.path = path
        with open(path+"bus.csv", "r") as f:
            rdr = csv.reader(f)
            self.bus = list(list(r) for r in csv.reader(f,delimiter=','))
        with open(path+"gen.csv", "r") as f:
            rdr = csv.reader(f)
            self.gen = list(list(r) for r in csv.reader(f,delimiter=','))
        with open(path+"branch.csv", "r") as f:
            rdr = csv.reader(f)
            self.branch = list(list(r) for r in csv.reader(f,delimiter=','))
        with open(path+"gencost.csv", "r") as f:
            rdr = csv.reader(f)
            self.gencost = list(list(r) for r in csv.reader(f,delimiter=','))

