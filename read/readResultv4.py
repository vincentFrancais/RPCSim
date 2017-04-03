from ctypes import *
#import numpy as np


class resultStruct(Structure):
    _fields_ = [('Dx', c_double),
    			('Dt', c_double),
                ('iNstep', c_int),
                ('thrCrossTimeStep', c_int),
                ('avalStatus', c_int),
                ('computeTime', c_double),
                ('streamer', c_int),
                ('nCluster', c_int),
                #('charges_size', c_uint),
                #('chargesTot_size', c_uint),
                #('signal_size', c_uint),
                ('size', c_uint),
                ('charges', c_double * 2000),
                ('chargesTot', c_double * 2000),
                ('signal', c_double * 2000),
                ('nions', c_double * 2000),
                ('pions', c_double * 2000),
                ('nelec', c_double * 2000),
                ('clPos', c_double * 2000),
                ('clNe', c_double * 2000)]


class RPCSimResult():
	
	def __init__(self, filename):
		self.filename = filename
		self.file = open(self.filename, 'rb')
		self.result = resultStruct()
		
	def nextResult(self):
		x = resultStruct()
		if self.file.readinto(x) == sizeof(x):
    			self.result = x
    			return True
    		else:
    			return False
