import pandas as pd
import matplotlib
matplotlib.use('TKagg')  #if user has different backend
import matplotlib.pyplot as plt



contour_line=10
t2=0
spectra_data=pd.read_csv("IR-hod.dat",header=4,sep=r'\s+')
vmin=spectra_data.I_total.min()
vmax=spectra_data.I_total.max()
print(vmax,vmin)
plt.xlabel('pump frequency$\omega_{1}/cm^{-1}$')
plt.ylabel('probe frequency$\omega_{3}/cm^{-1}$')
plt.ylabel('relative intensity')
plt.title(f"waitting time t2 = {t2}")

C=plt.tricontour(spectra_data.w1,spectra_data.w3, spectra_data.I_total , contour_line, linewidths=0.5, colors='k')
plt.clabel(C, inline=True, fontsize=12)
plt.tricontourf(spectra_data.w1,spectra_data.w3, spectra_data.I_total,contour_line,cmap='RdBu_r',
        #norm=matplotlib.colors.Normalize(vmin,-vmin)
        norm=matplotlib.colors.Normalize(-vmax,vmax)
        )
plt.colorbar()
plt.show()
