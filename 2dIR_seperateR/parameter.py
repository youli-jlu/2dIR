import pandas as pd
import numpy as np

def main():
    fpara=open("parameter.dat",'w')
    
    freqfile="energy-tot"
    dt=0.04
    tgap=10
    omega_min=-200
    omega_max=200
    momega1=2710.662
    momega2=2596.495



class spectra:
    def __init__( self,**kwargs):
        self.responseflag= = kwargs.get("responseflag","2dresponse")
        self.fftflag =  kwargs.get("fftflag","2dfft")
        self.freqfile=  kwargs.get("frequencyfile","freq.dat")
        self.tdfile=  kwargs.get("tdfile","")
        self.dimension=  kwargs.get("dimension","2")
        self.t2=  kwargs.get("waitingtime","")


