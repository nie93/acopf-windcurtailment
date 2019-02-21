from adapter import C37118InputDataAdapter
from case import Case, Const
from compute_interface import *
import opf
from time import time
from pdb import set_trace
from pprint import pprint


def main():
    c37118_data_adapter = C37118InputDataAdapter()
    c37118_data_adapter.add_connection({"ip": "192.168.1.119", "port": 4712, "id": 9})
    c37118_data_adapter.add_connection({"ip": "192.168.1.119", "port": 4722, "id": 9})
    c37118_data_adapter.connect_pmu()
    c37118_data_adapter.close()
    frame = c37118_data_adapter.get_pmu_measurements()

    pprint(frame)

    pg = calculatePg(frame)

    print(pg)
    set_trace()


    const = Const()
    casepath = './case14mod/'
    c = Case()
    c.import_case(casepath)
    c.set_gen_prop(const.PMAX, [1,2,3], pg[1:4])
    c.set_branch_prop('RATE', [14], [34.999])
    c.scale_branch_prop([const.BR_R, const.BR_X], multi=1.0)    
    
    for i in range(0,1):
        start_time = time()
        opf.runcopf(c, flat_start=False)
        end_time = time()

        print('Optimization execution time: %.8f' % (end_time - start_time))

if __name__ == "__main__":
    main()
