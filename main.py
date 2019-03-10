from adapter import C37118InputDataAdapter, RtdsCmdAdapter
from case import Case, Const
from compute_interface import *
import opf
import time
from pdb import set_trace
from pprint import pprint


def main():
    rtds_cmd_adapter = RtdsCmdAdapter({"ip": "192.168.1.104", "port": 4575})
    rtds_cmd_adapter.connect_rtds()


    c37118_data_adapter = C37118InputDataAdapter()
    c37118_data_adapter.add_connection({"ip": "192.168.1.111", "port": 4712, "id": 1})
    c37118_data_adapter.add_connection({"ip": "192.168.1.111", "port": 4722, "id": 2})
    c37118_data_adapter.connect_pmu()
    c37118_data_adapter.close()

    frame = c37118_data_adapter.get_pmu_measurements()
    currentPg = calculatePg(frame)
    pprint(frame)
    print(currentPg)

    const = Const()
    casepath = './case14mod/'
    c = Case()
    c.import_case(casepath)
    c.set_gen_prop(const.PMAX, [1,2,3], currentPg[1:4])
#     c.set_gen_prop(const.PMAX, [1,2,3], [40.05814367749094, 70.12544088354686, 78.47923733719254])
#     c.set_branch_prop('RATE', [14], [34.999])
    c.scale_branch_prop([const.BR_R, const.BR_X], multi=1.0)    
    
    for i in range(0, 100):
        start_time = time.time()
        res = opf.runcopf(c, flat_start=False)
        end_time = time.time()
        print("Optimal Outputs: %s" % str(res['PG']))
        print('Optimization execution time: %.8f' % (end_time - start_time))
        
        if res['PG'][2] < 0.99*currentPg[2]:
                rtds_cmd_adapter.send_cmd('SetSlider "SL3" = %.4f;' % (res['PG'][2] / 100))

    rtds_cmd_adapter.close()

if __name__ == "__main__":
    main()
