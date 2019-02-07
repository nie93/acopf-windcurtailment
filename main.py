import os
import numpy
import scipy
from case import Case, Const
import pdb

def main():
    const = Const()
    casepath = r"""C:\\Users\\niezj\\Documents\\wsu\\Research\\ARPAE\\acopf\\acopf_windcurtail_python\\case14mod\\"""
    c = Case()
    c.import_case(casepath)
    print(c.gencost)



if __name__ == "__main__":
    main()