from pypmu.pdc import Pdc
from pypmu.frame import CommandFrame
# from pypmu.frame import DataFrame
import threading
import struct
import socket
from pdb import set_trace
from pprint import *
import time


class RtdsCmdAdapter(object):

    def __init__(self, connection_dict):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.conn_dict = connection_dict
        
    def connect_rtds(self):
        self.s.connect((self.conn_dict['ip'], self.conn_dict['port'])) # RSCAD socket location
        self.send_cmd('Start;')
        self.send_cmd('ListenOnPortHandshake("temp_string");')

        rmsg = self.s.recv(64).decode("utf-8")
        while ('temp_string' not in rmsg):
            rmsg = self.s.recv(64).decode("utf-8")
        print("*    Successfully connected to RSCAD Script TCP Server.")

    def send_cmd(self, cmdstr):
        print("*    Sending command: %s." % cmdstr)
        self.s.send(cmdstr.encode())

    def close(self):
        self.s.shutdown(socket.SHUT_RDWR)
        self.s.close()
        print("*    Successfully closed RSCAD Script TCP Server connection.")



# class C37118Thread(threading.Thread):    

class C37118InputDataAdapter(object):

    def __init__(self):
        self.s = None
        self.pdc = None
        self.phnmr = dict()
        self.dframes = dict()
        self.chnme = dict()
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

        self.get_dframes()
        # set_trace()

    def get_dframes(self):
        if self.pdc is None:
            return None

        _dframes = dict()
        for pmu in self.pdc:
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            _dframes[pmukey] = pmu.get()

        # self.lock.acquire()
        self.dframes = _dframes
        # self.lock.release()


    def get_pmu_measurements(self):
        if not self.dframes:
            return None
        
        return self.parse_dframes()

    def parse_config(self, pmu):
        if not pmu:
            return None

        pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
        self.chnme[pmukey] = list()
        config = pmu.get_config()
        phnmr = struct.unpack('>h', config[40:42])[0]
        self.phnmr[pmukey] = phnmr

        chnam = ""
        for iph in range(phnmr):
            chnam = ""
            for i in range(46 + 16 * iph, 46 + 16 * (iph+1)):
                chnam += struct.unpack('>c', config[i:i+1])[0].decode("utf-8")
            self.chnme[pmukey].append(chnam.replace(" ", ""))

        return

    def parse_dframes(self):
        if not self.dframes:
            return None

        _dframes = self.dframes
        data = dict()
        time.sleep(1)
        for pmu in self.pdc:
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            data[pmukey] = dict()
            
            for iph in range(self.phnmr[pmukey]):
                magkey = self.chnme[pmukey][iph] + "_mag"
                angkey = self.chnme[pmukey][iph] + "_ang"
                data[pmukey][magkey] = struct.unpack('>f', _dframes[pmukey][16 + 8*iph:20 + 8*iph])[0]
                data[pmukey][angkey] = struct.unpack('>f', _dframes[pmukey][20 + 8*iph:24 + 8*iph])[0]
                      
        return data

    def close(self):
        for pmu in self.pdc:
            pmu.stop()
            pmu.quit()


        





