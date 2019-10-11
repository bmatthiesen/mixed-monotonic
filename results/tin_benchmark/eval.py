import h5py
import numpy as np
import pandas as pd

# ratespace
# wp 143: TO 7d
# wp 155: TO 7d
# wp 173: TO 7d

cutoff = 6

f = h5py.File('result.h5','r')
gee = ('gee_mmp','gee_dinkelbach')
tin = ('dm','mmp','mmpFifo','pa','ratespace','mapel')

def extract(keys, fnp):
    ds = {'iterations': {}, 'runtime': {}, 'memory': {}}
    for k in keys:
        if k == 'mapel':
            sel = np.sum(np.isnan(f[k]['iterations']),axis=-1) > 13 # I'm not exactly sure why these results are not here so we're just going to ignore them (I'm quite sure that they ran longer than 7 days and where cancelled but I'm unable to confirm that)
        else:
            sel = np.sum(np.isnan(f[k]['iterations']),axis=-1) > cutoff

        for d in ds:
            tmp = f[k][d][...]
            if k == "ratespace" and d == "runtime":
                tmp[2,43] = 7 * 24 * 3600
                tmp[2,55] = 7 * 24 * 3600
                tmp[2,73] = 7 * 24 * 3600

            tmp = np.nanmean(tmp, axis=-1)

            tmp[sel] = np.nan
            ds[d][k] = tmp

    for d in ds:
        l = np.min([k.shape[0] for k in ds[d].values()])
        tmp = pd.DataFrame(ds[d], pd.RangeIndex(1, l+1))
        tmp.to_csv(fnp + d + '.dat')

if __name__=="__main__":
    extract(gee,"gee_")
    extract(tin,"tin_")
