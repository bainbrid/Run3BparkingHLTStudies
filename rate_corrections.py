import numpy as np
import pandas as pd
import json
from common import path

pd.options.display.float_format = '{:8.2f}'.format
pd.set_option('display.width', 100)

################################################################################
# Sebastian's rates

df = pd.read_csv(path+'rates/L1SkimRevised_compare.csv',index_col='Path')

usecols = [ col for col in df.columns if 'seeds active' in col ]
dropcols = [ col for col in df.columns if col not in usecols ]

shortrows = [ x for x in np.arange(4.,10.,0.5) ]
shortcols = [ x for x in np.arange(5.,11.,0.5) ]

df.drop(columns=dropcols,inplace=True)
df.index = shortrows
df.columns = shortcols
print
print("Columns = L1 pT threshold (GeV)")
print("Rows    = HLT pT threshold (GeV)")
print("Entries = HLT trigger rate (Hz)")

print
print("Sebastian's rates:")
print(df)

linst_seb = 1.385
for col in df.columns: df[col].values[:] *= (linst_seb/2.0) # Scale linearly down from 2.0 to 1.385
print
print("Sebastian's rates @ Linst = {:.3f}E34:".format(linst_seb))
print(df)

dct = df.to_dict() # {L1:{HLT:rate,...},...}
#print(dct)
json_string = df.to_json(path+'rates/rates.json') # {L1:{HLT:rate,...},...}
#print(json_string)

################################################################################
# Yuta's rates

npu = 36 #56
infile = open(path+'rates/rates_{:.0f}.json'.format(npu),'r')
dct1 = json.load(infile)
#print(dct1)

df1 = pd.DataFrame.from_dict(dct1)
df1.columns = df1.columns.astype(float)
df1 = df1.reindex(sorted(df1.columns), axis=1)
df1.index = df1.index.astype(float)
df1 = df1.reindex(sorted(df1.index), axis=0)
#df1.drop(columns=[4.0,4.5],inplace=True)
print
print("Yuta's rates (npu={:.0f}):".format(npu))
print(df1)

# Linst = npu * 0.0357338 - 0.0011904
npu_1385 = (linst_seb + 0.0011904)/0.0357338
for col in df1.columns: df1[col].values[:] *= (npu_1385/npu) # Scale to Linst = 1.385
print
print("Yuta's rates @ npu={:.1f}:".format(npu_1385))
print(df1)

################################################################################
# Ratio of rates

df2 = df/df1
print
print("Ratio of rates (Sebastian/Yuta):")
print(df2)

print
print("Summary statistics per column:")
print(df2.describe())

print
print("Summary statistics per row:")
print(df2.T.describe().T)

print
print("Correction factor to apply to Yuta's rates: {:.2f}".format(df2.stack().mean()))
print("std. dev. in correction factor:             {:.2f}".format(df2.stack().std()))

################################################################################
# Corrections to rates

factor = 2.17

df3 = df2.copy()
for col in df3.columns: df3[col].values[:] = factor
print
print("Ratios to be applied for npu = 56 (Linst = 2E34):")
print(df3)

for npu in [56, 48, 42, 36, 30, 25, 17]:
    corr = factor # (factor-1.) * (npu*1./56.) + 1. # <-- const factor, don't scale linearly with Linst
    for col in df3.columns: df3[col].values[:] = corr # scale linearly according to npu
    filename = path+'rates/corrections_{:.0f}.json'.format(npu)
    json_string1 = df3.to_json(filename)
