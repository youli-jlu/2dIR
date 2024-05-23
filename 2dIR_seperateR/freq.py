import pandas as pd

raw = pd.read_csv("freq-all.dat",delimiter=r'\s+')
raw.insert(1,'freq2',raw.freq1)

raw.freq1/=0.188
raw.freq2/=0.188
raw.freq1+=0.0
raw.freq2+=-240.0

raw.to_csv("freq-hod.dat",sep=" ",index=False)
