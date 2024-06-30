#!/usr/bin/python3.6
from spectra2dSFG import  pyspectra

def test1d():
    spectra = pyspectra(
            dimension=1,
            response=b"cumulant",
    # energy unit: cm-1
            freqfile=b"freq-all.dat",
            tdfile1=b"td1.dat",
            tpfile1=b"tp1.dat",
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
            spectrafile1=b"SFG-cy.dat",
            responsefile=b"1dresponse"
            )

    spectra.calculation()
    spectra.show_SFG()
    spectra.show_response()

def test2d():
    spectra = pyspectra(
            dimension=2,
            response=b"non-condon",
    # energy unit: cm-1
            freqfile=b"freq-hod.dat",
            tdfile1=b"td1.dat",
            tdfile2=b"td2.dat",
            tpfile1=b"tp1.dat",
            tpfile2=b"tp2.dat",
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
            spectrafile1=b"SFG1-hod.dat",
            spectrafile2=b"SFG2-hod.dat",
            responsefile=b"2dresponse"
            )

    spectra.calculation()
    spectra.show_SFG()
    spectra.show_response()

test1d()
test2d()
