from case import Case, Const
import opf
from pdb import *


def main():
    const = Const()
    casepath = './case14mod/'
    c = Case()
    c.import_case(casepath)
    # c.set_gen_prop(const.PMAX, [1,2,3], [30, 80, 80])
    # c.set_branch_prop('RATE', [14], [34.999])
    c.scale_branch_prop([const.BR_R, const.BR_X], multi=1.0)
    
    opf.runcopf(c, flat_start=True)


if __name__ == "__main__":
    main()
