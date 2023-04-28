#!/usr/bin/python3.6
from spectra2d import  pyspectra


spectra = pyspectra(
        dimension=1,
        responseflag=b"condon",
# energy unit: cm-1
        freqfile=b"freq-hod.dat",
        spectrafile=b"IR-cy.dat",
        tdfile=b"td.dat",
# time unit: fs
        dt=10.0,
        t2=0,
# int types, in /(dt) unit 
        readgap=4,
        dt_time_average=10,
        zeropadding=10,
        response_tmax=20,
        maxfreq=300,
        minfreq=-300,
        trans1relaxation=1.7,
        trans2relaxation=10
        )

spectra.calculation()
