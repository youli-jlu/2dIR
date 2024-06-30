from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "2dSFG.cc":
    pass
cdef extern from "1dSFG.cc":
    pass
cdef extern from "response2d.cc":
    pass
cdef extern from "response1d.cc":
    pass
cdef extern from "dft.cc":
    pass
cdef extern from "class.cc":
    pass

cdef extern from "class.h" :
    cdef cppclass para:
        para() except +
        #spectra category and response function
        int dimension 
        string responseflag
        string spectraflag
        string freqfile
        string tdfile1
        string tdfile2
        string tpfile1
        string tpfile2
        string spectrafile1
        string spectrafile2
        string responsefile 

        #basic parameter
        int n    #total time step
        double dt#time interval between adjacent step
        double rtgap #read gap to save time
        double omega_min,omega_max #frequency range in 2d-SFG simulation
        double momega1,momega2 #frequency range in 2d-SFG simulation
        double dw #frequency range in 2d-SFG simulation
        double domega #momega2-momega1, with unit trans
        double t2  #t2 time in 2d-SFG simulation
        double tgap#time interval between adjacent ensemble
        double nzp # number of zero padding
        double trans1T1;
        double trans2T1;
        # ensemble parameter
        int tmax
        int np


cdef extern from "2dSFG.h" :
    cdef cppclass spectra2d:
        spectra2d(para*) except +
        para parameter
        void calculation()

cdef extern from "1dSFG.h" :
    cdef cppclass spectra1d:
        spectra1d(para*) except +
        para parameter
        void calculation()
