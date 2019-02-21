from synchrophasor.pdc import Pdc
from synchrophasor.frame import DataFrame
import threading
from pdb import set_trace

# class C37118Thread(threading.Thread):    

class C37118InputDataAdapter(object):

    def __init__(self):
        self.s = None
        self.pdc = None
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
                pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
                self.channel_names[pmukey] = [chkey.replace(" ", "") for chkey in pmu.get_config().get_channel_names()]                     
                self.pdc.append(pmu)
            except:
                pass
        
        for pmu in self.pdc:
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
                data = pmu.get()
                _dframes.append(data.get_measurements()["measurements"][0]["phasors"])

            self.lock.acquire()
            self.dframes = _dframes[:]
            self.lock.release()


    def get_pmu_measurements(self):
        if not self.dframes:
            return None
        
        _dframes = self.dframes[:]
        
        phasors = self.parse_dframes(_dframes)
        
        return phasors

    def parse_dframes(self, dframes):
        if not dframes:
            return None

        data = dict()
        for i, pmu in enumerate(self.pdc):
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            data[pmukey] = dict()
            for chi, chkey in enumerate(self.channel_names[pmukey]):
                magkey = chkey + "_mag"
                angkey = chkey + "_ang"
                data[pmukey][magkey] = dframes[i][chi][0]
                data[pmukey][angkey] = dframes[i][chi][1]

        return data

    def close(self):
        for pmu in self.pdc:
            pmu.stop()
            pmu.quit()


        





