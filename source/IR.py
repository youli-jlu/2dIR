#!/usr/bin/python3.6
from spectra2d import  pyspectra


#x=freq2-2-400-mono2
spectra = pyspectra(
        dimension=1,
        response=b"cummulant",
# energy unit: cm-1
        freqfile=b"freq1-2-400-mono1",
        spectrafile=b"IR-mono1-x2.dat",
        #freqfile=b"freq1-nma2-1d-r9-xtb2",
        #spectrafile=b"IR-nma2-1d-x2.dat",
        tdfile=b"td1.dat",
# time unit: fs
        dt=10.0,
# int types, in /(dt) unit 
        t2=0,
        readgap=1,
        dt_time_avarege=1,
        zeropodding=1000,
        #zeropadding=1500,
        response_tmax=400,
        maxfreq=150,
        minfreq=-150,
        trans1relaxation=0.85,
        )

spectra.calculation()
spectra.show_IR()
spectra.show_response()
