from pypmu.pdc import Pdc
from pypmu.frame import CommandFrame
# from pypmu.frame import DataFrame
import threading
import struct
from pdb import set_trace
from pprint import *


# class C37118Thread(threading.Thread):    

class C37118InputDataAdapter(object):

    def __init__(self):
        self.s = None
        self.pdc = None
        self.phnmr = dict()
        self.dframes = dict()
        self.channel_names = dict()
        self.connection_dicts = list()
        self.lock = threading.Lock()

    def add_connection(self, conn_dict):
        self.connection_dicts.append(conn_dict)

    def connect_pmu(self):
        self.pdc = list()
        for conn_dict in self.connection_dicts:
            pmu = Pdc(pdc_id=conn_dict["id"], pmu_ip=conn_dict["ip"], pmu_port=conn_dict["port"])
            pmu.logger.setLevel("DEBUG")
            try:
                pmu.run()
                self.pdc.append(pmu)
            except:
                pass
        
        for pmu in self.pdc:
            self.parse_config(pmu)
            pmu.start()

        # dframe_thread = threading.Thread(target=self.get_dframes)
        # dframe_thread.start()
        self.get_dframes()

    def get_dframes(self):
        if self.pdc is None:
            return False

        for i in range(0, 30):
            _dframes = list()
            for pmu in self.pdc:
                _dframes.append(pmu.get())

            self.lock.acquire()
            self.dframes = _dframes[:]
            self.lock.release()


    def get_pmu_measurements(self):
        if not self.dframes:
            return None
        
        _dframes = self.dframes[:]
        phasors = self.parse_dframes(_dframes[0])
        
        return phasors

    def parse_config(self, pmu):
        if not pmu:
            return None

        pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
        self.channel_names[pmukey] = list()
        config = pmu.get_config()
        phnmr = struct.unpack('>h', config[40:42])[0]
        self.phnmr[pmukey] = phnmr

        chnam = ""
        for iph in range(phnmr):
            chnam = ""
            for i in range(46 + 16 * iph, 46 + 16 * (iph+1)):
                chnam += struct.unpack('>c', config[i:i+1])[0].decode("utf-8")
            self.channel_names[pmukey].append(chnam.replace(" ", ""))

        return

    def parse_dframes(self, dframes):
        if not dframes:
            return None

        data = dict()
        for pmu in self.pdc:
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            data[pmukey] = dict()

            for iph in range(self.phnmr[pmukey]):
                magkey = self.channel_names[pmukey][iph] + "_mag"
                angkey = self.channel_names[pmukey][iph] + "_ang"
                data[pmukey][magkey] = struct.unpack('>f', dframes[16 + 8*iph:20 + 8*iph])[0]
                data[pmukey][angkey] = struct.unpack('>f', dframes[20 + 8*iph:24 + 8*iph])[0]
                
        return data

    def close(self):
        for pmu in self.pdc:
            pmu.stop()
            pmu.quit()


        





