import h5py
import numpy as np

f = h5py.File('result3.h5', 'r')

objval = None

for ds in f:
    assert(~np.any(f[ds]['runtime'] == 0))
    print('{}: {:.3f} / {:.3f}'.format(ds, np.mean(f[ds]['runtime']), np.median(f[ds]['runtime'])))

    if objval is None:
        objval = f[ds]['objval'][:]
    else:
        assert(np.allclose(objval, f[ds]['objval'], atol = 1e-3))


print('\n\nAloha 4')
g = h5py.File('result4.h5','r')

objval = None
for ds in ['MMP','MMPred']:
    assert(~np.any(g[ds]['runtime'] == 0))
    print('{}: {:.3f} / {:.3f}'.format(ds, np.mean(g[ds]['runtime']), np.median(g[ds]['runtime'])))

    if objval is None:
        objval = g[ds]['objval'][:]
    else:
        assert(np.allclose(objval, g[ds]['objval'], atol = 1e-3))

gp = g['GP']['runtime']
gp[gp==0] = 7*24*3600 # TO after 7 days
print('{}: {:.3f} / {:.3f}'.format('GP', np.mean(gp), np.median(gp)))
assert(np.allclose(objval, g[ds]['objval'], atol = 1e-3))
