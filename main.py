import os
import numpy
import scipy
from case import Case, Const
import opf
import pdb

def main():
    const = Const()
    casepath = r"""C:\\Users\\niezj\\Documents\\wsu\\Research\\ARPAE\\acopf\\acopf_windcurtail_python\\case14mod\\"""
    c = Case()
    c.import_case(casepath)
    # print(c.branchrate.shape[0])
    
    opf.runcopf(c)



if __name__ == "__main__":
    main()