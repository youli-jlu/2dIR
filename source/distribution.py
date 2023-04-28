import pandas as pd
import numpy as np

rawdata= pd.read_csv("freq-hod.dat",delimiter=r'\s+')

print(rawdata.w1.max())
print(rawdata.w1.min())

dis=open("distribution-hod",'w')
dt=2
for omega in np.arange(rawdata.w1.min(),rawdata.w1.max(),0.5):
    raw_condition=rawdata[ (rawdata['w1']<omega)
            &(rawdata['w1']>omega-dt)]
    n1=raw_condition.w1.count()
    dis.write(f"{omega}  {n1}    \n")
dis.close()


