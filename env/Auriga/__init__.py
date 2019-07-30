#origin from https://bitbucket.org/illustris/illustris_python
#modified by Yingjie
__all__ = ["groupcat", "snapshot", "util", "sublink", "sublink_gal", "lhalotree"]

from . import *

import h5py
def loadheader(basePath,snapNum,name):
    '''
    name:
    'snap' or 'group' equal 'subhalo'
    '''
    if name=='snap':
        fname= snapshot.snapPath(basePath,snapNum)
    elif name in ['suhalo','group']:
        fname= groupcat.gcPath(basePath,snapNum)
    with h5py.File(fname,'r') as f:
        header = dict( f['Header'].attrs.items() )
    return header