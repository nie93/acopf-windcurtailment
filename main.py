from case import Case, Const
import opf
from pdb import *


def main():
    const = Const()
    casepath = r"""C:\\Users\\niezj\\Documents\\wsu\\Research\\ARPAE\\acopf\\acopf_windcurtail_python\\case14mod\\"""
    c = Case()
    c.import_case(casepath)
    c.set_gen_prop(const.PMAX, [1,2,3], [30, 80, 80])
    c.set_branch_prop('RATE', [14], [34.999])
    c.scale_branch_prop([const.BR_R, const.BR_X], multi=1.0)    
    
    for i in range(0,20):
        opf.runcopf(c, flat_start=False)


if __name__ == "__main__":
    main()