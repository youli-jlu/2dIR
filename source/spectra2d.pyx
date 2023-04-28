# distutils: libraries = fftw3 
# distutils: language = c++

from IR2d cimport para
from IR2d cimport spectra2d
from IR2d cimport spectra1d
from libcpp.string cimport string
from libcpp.vector cimport vector


#cdef class pypara:
#    cdef para* cparameter
#
#    def __cinit__(self):
#        self.cparameter=new para()
#    def __dealloc__(self):
#        del self.cparameter
#
class pyspectra:
    def __init__(self,**kwargs):
        self.responseflag=kwargs.get("response",b"condon")
        self.spectraflag =kwargs.get("spectra",b"2d")
        self.freqfile    =kwargs.get("freqfile",b"freq.dat")
        self.spectrafile =kwargs.get("spectrafile",b"2dspectra")
        self.tdfile      =kwargs.get("tdfile",b"td.dat")
        self.dimension   =int(kwargs.get("dimension",2))
        self.t2          =int(kwargs.get("t2",0))
        self.rtgap       =int(kwargs.get("readgap",1))
        self.tgap        =int(kwargs.get("dt_time_average",10))
        self.nzp         =int(kwargs.get("zeropodding",0))
        self.tmax        =int(kwargs.get("response_tmax",0))
        self.dt          =float(kwargs.get("dt",0.01))
        self.omega_max   =float(kwargs.get("maxfreq",0))
        self.omega_min   =float(kwargs.get("minfreq",0))
        self.trans1T1    =float(kwargs.get("trans1relaxation",0))
        self.trans2T1    =float(kwargs.get("trans2relaxation",0))


    def calculation1d(self):
        parameter=new para()
        parameter.responseflag=self.responseflag
        parameter.spectraflag =self.spectraflag
        parameter.freqfile    =self.freqfile
        parameter.spectrafile =self.spectrafile 
        parameter.tdfile      =self.tdfile      
        parameter.dimension   =self.dimension   
        parameter.t2          =self.t2
        parameter.dt          =self.dt          
        parameter.tmax        =self.tmax
        parameter.rtgap       =self.rtgap       
        parameter.tgap        =self.tgap        
        parameter.nzp         =self.nzp         
        parameter.omega_max   =self.omega_max   
        parameter.omega_min   =self.omega_min   
        parameter.trans1T1    =self.trans1T1    
        parameter.trans2T1    =self.trans2T1    

        spectra=new spectra1d(parameter)
        try:
            spectra.calculation()
        finally:
            del spectra
            print("end of spectra calculation\n")
    def calculation2d(self):
        parameter=new para()
        parameter.responseflag=self.responseflag
        parameter.spectraflag =self.spectraflag     
        parameter.freqfile    =self.freqfile    
        parameter.spectrafile =self.spectrafile 
        parameter.tdfile      =self.tdfile      

        parameter.dimension   =self.dimension   
        parameter.t2          =self.t2
        parameter.dt          =self.dt          
        parameter.rtgap       =self.rtgap       
        parameter.tmax        =self.tmax
        parameter.tgap        =self.tgap        
        parameter.nzp         =self.nzp         
        parameter.omega_max   =self.omega_max   
        parameter.omega_min   =self.omega_min   
        parameter.trans1T1    =self.trans1T1    
        parameter.trans2T1    =self.trans2T1    

        spectra=new spectra2d(parameter)
        try:
            spectra.calculation()
        finally:
            del spectra
            print("end of spectra calculation\n")

    def calculation(self):
        if (self.dimension==1):
            print(" doing 1d spectra")
            self.calculation1d();
        elif (self.dimension==2):
            print(" doing 2d spectra")
            self.calculation2d();
        else:
            print("wrong dimension")
            exit()
