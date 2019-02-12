from case import Case, Const
import opf
from pdb import *


def main():
    const = Const()
    casepath = r"""C:\\Users\\niezj\\Documents\\wsu\\Research\\ARPAE\\acopf\\acopf_windcurtail_python\\case14mod\\"""
    c = Case()
    c.import_case(casepath)
    # c.set_gen_prop('PMAX', [1,2,3], [30, 80, 80])
    c.set_branch_prop('RATE', [14], [34.999])
    
    opf.runcopf(c)


if __name__ == "__main__":
    main()