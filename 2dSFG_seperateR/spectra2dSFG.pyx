# distutils: libraries = fftw3 
# distutils: language = c++

from SFG2d cimport para
from SFG2d cimport spectra2d
from SFG2d cimport spectra1d
from libcpp.string cimport string
from libcpp.vector cimport vector

import pandas as pd


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
        self.tdfile1      =kwargs.get("tdfile1",b"td1.dat")
        self.tdfile2      =kwargs.get("tdfile2",b"td2.dat")
        self.tpfile1      =kwargs.get("tpfile1",b"tp1.dat")
        self.tpfile2      =kwargs.get("tpfile2",b"tp2.dat")
        self.dimension   =int(kwargs.get("dimension",2))
        self.t2          =int(kwargs.get("t2",0))
        self.rtgap       =int(kwargs.get("readgap",1))
        self.tgap        =int(kwargs.get("dt_time_average",10))
        self.nzp         =int(kwargs.get("zeropadding",0))
        self.tmax        =int(kwargs.get("response_tmax",0))
        self.dt          =float(kwargs.get("dt",0.01))
        self.omega_max   =float(kwargs.get("maxfreq",0))
        self.omega_min   =float(kwargs.get("minfreq",0))
        self.trans1T1    =float(kwargs.get("trans1relaxation",0))
        self.trans2T1    =float(kwargs.get("trans2relaxation",0))
        #output
        self.spectrafile1 =kwargs.get("spectrafile1",b"2dspectra")
        self.spectrafile2 =kwargs.get("spectrafile2",b"2dspectra")
        self.responsefile=kwargs.get("responsefile",b"response")


    def calculation1d(self):
        parameter=new para()
        parameter.responseflag=self.responseflag
        parameter.freqfile    =self.freqfile
        parameter.tdfile1      =self.tdfile1      
        parameter.tdfile2      =self.tdfile2      
        parameter.tpfile1      =self.tpfile1      
        parameter.tpfile2      =self.tpfile2      
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
        parameter.spectraflag =self.spectraflag

        #parameter.spectrafile1 =self.spectrafile1 
        #parameter.spectrafile2 =self.spectrafile2
        parameter.spectrafile1 =self.spectrafile1
        parameter.responsefile=self.responsefile

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
        parameter.tdfile1      =self.tdfile1      
        parameter.tdfile2      =self.tdfile2      
        parameter.tpfile1      =self.tpfile1      
        parameter.tpfile2      =self.tpfile2      

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

        parameter.spectrafile1 =self.spectrafile1
        parameter.spectrafile2 =self.spectrafile2
        parameter.responsefile=self.responsefile

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

    def show_SFG(self,**kwargs):
        contour_line=kwargs.get("nline",10)
        import matplotlib
        matplotlib.use(kwargs.get("backend",'TKagg'))  #if user has different backend
        import matplotlib.pyplot as plt
        if (self.dimension==1):
            file1=kwargs.get("SFGfile1",self.spectrafile1.decode())
            spectra_data=pd.read_csv(file1,header=4,sep=r'\s+')
            plt.xlabel('frequency/$cm^{-1}$')
            plt.ylabel('relative intensity')
            plt.plot(spectra_data.w1,spectra_data.I)
            plt.scatter(spectra_data.w1,spectra_data.I,marker='X')
            plt.show()
        if (self.dimension==2):
            def plot_once(spectra_data):
                #plot coutour figure
                plt.xlabel('pump frequency$\omega_{1}/cm^{-1}$')
                plt.ylabel('probe frequency$\omega_{3}/cm^{-1}$')
                plt.title(f"waitting time t2 = {self.t2}")
                #region of our contour 
                if ( (-spectra_data.I_total.min())>=spectra_data.I_total.max()):
                    vmax=-spectra_data.I_total.min()
                else:
                    vmax=spectra_data.I_total.max()
                C=plt.tricontour(spectra_data.w1,spectra_data.w3, spectra_data.I_total , contour_line, linewidths=0.5, colors='k')
                plt.clabel(C, inline=True, fontsize=12)
                plt.tricontourf(spectra_data.w1,spectra_data.w3, spectra_data.I_total,contour_line,cmap='RdBu_r',
                        norm=matplotlib.colors.Normalize(-vmax,vmax)
                        )
                plt.colorbar()
                plt.show()
            file1=kwargs.get("SFGfile1",self.spectrafile1.decode())
            file2=kwargs.get("SFGfile2",self.spectrafile2.decode())
            spectra_data1=pd.read_csv(file1,header=4,sep=r'\s+')
            spectra_data2=pd.read_csv(file2,header=4,sep=r'\s+')
            plot_once(spectra_data1)
            plot_once(spectra_data2)
            
    def show_response(self,**kwargs):
        filename=kwargs.get("rfile",self.responsefile.decode())
        contour_line=kwargs.get("nline",10)
        import matplotlib
        matplotlib.use(kwargs.get("backend",'TKagg'))  #if user has different backend
        import matplotlib.pyplot as plt
        spectra_data=pd.read_csv(filename,header=0,sep=r'\s+')
        if (self.dimension==1):
            plt.xlabel('t / ps ')
            plt.ylabel('response function')
            plt.plot(spectra_data.t,spectra_data.re,label='re')
            plt.scatter(spectra_data.t,spectra_data.re,marker='X')
            plt.plot(spectra_data.t,spectra_data.im,label='im')
            plt.scatter(spectra_data.t,spectra_data.im,marker='X')
            plt.legend()
            plt.show()
        if (self.dimension==2):
            #plot coutour figure
            def plotc(x,y,z,title):
                plt.xlabel('$t_{1}  / ps$')
                plt.ylabel('$t_{3}  / ps$')
                ##region of our contour 
                #if ( (-spectra_data.I_total.min())>=spectra_data.I_total.max()):
                #    vmax=-spectra_data.I_total.min()
                #else:
                #    vmax=spectra_data.I_total.max()
                plt.title(f"{title} waitting time t2 = {self.t2}")
                C=plt.tricontour(x,y,z, contour_line, linewidths=0.5, colors='k')
                plt.clabel(C, inline=True, fontsize=12)
                plt.tricontourf(x,y,z,contour_line,cmap='RdBu_r',
                        #norm=matplotlib.colors.Normalize(-vmax,vmax)
                        )
                plt.colorbar()
                plt.show()

            plotc(spectra_data.t1,spectra_data.t3, spectra_data.r1_re, "response1 re rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.r1_im, "response1 im rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.nr1_re,"response1 re non-rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.nr1_im,"response1 im non-rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.r2_re, "response2 re rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.r2_im, "response2 im rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.nr2_re,"response2 re non-rephasing")
            plotc(spectra_data.t1,spectra_data.t3, spectra_data.nr2_im,"response2 im non-rephasing")
            
