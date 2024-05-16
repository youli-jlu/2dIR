#!/usr/bin/python3.6
from spectra2d import  pyspectra

def test1d():
    spectra = pyspectra(
            dimension=1,
            response=b"cummulant",
    # energy unit: cm-1
            freqfile=b"freq-all.dat",
            tdfile1=b"td1.dat",
            #tdfile2=b"td2.dat",
    # time unit: fs
            dt=10.0,
    # time unit: ps
            trans1relaxation=0.7,
    # int types, in /(dt) unit 
            t2=0,
            readgap=4,
            dt_time_average=10,
            zeropadding=10,
            response_tmax=200,
            maxfreq=300,
            minfreq=-300,
            #output 
            spectrafile=b"IR-cy.dat",
            responsefile=b"1dresponse"
            )
    
    spectra.calculation()
    spectra.show_IR()
    spectra.show_response()

def test2d():
    spectra = pyspectra(
            dimension=2,
            response=b"non-condon",
    # energy unit: cm-1
            freqfile=b"freq-hod.dat",
            tdfile1=b"td1.dat",
            tdfile2=b"td2.dat",
    # time unit: fs
            dt=10.0,
    # time unit: ps
            trans1relaxation=0.7,
            trans2relaxation=0.7,
    # int types, in /(dt) unit 
            t2=0,
            readgap=4,
            dt_time_average=10,
            zeropadding=10,
            response_tmax=20,
            maxfreq=300,
            minfreq=-300,
            #output file
            spectrafile=b"IR-hod.dat",
            responsefile=b"2dresponse"
            )
    
    spectra.calculation()
    spectra.show_IR()
    spectra.show_response()

test2d()
#test1d()
