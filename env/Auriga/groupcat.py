""" Illustris Simulation: Public Data Release.
groupcat.py: File I/O related to the FoF and Subfind group catalogs. """
from __future__ import print_function

import six
from os.path import isfile
import numpy as np
import h5py
import os


def gcPath(basePath, snapNum, chunkNum=0):
    """ Return absolute path to a group catalog HDF5 file (modify as needed). """
    gcPath = basePath + '/groups_%03d/' % snapNum
    filePath1 = gcPath + 'groups_%03d.%d.hdf5' % (snapNum, chunkNum)
    filePath2 = gcPath + 'fof_subhalo_tab_%03d.%d.hdf5' % (snapNum, chunkNum)

    if isfile(filePath1):
        return filePath1
    return filePath2


def snapfileoff(basePath,snapNum):
    result={}
    def snapPath(basePath,snapNum,chunkNum=0):
        """ Return absolute path to a snapshot HDF5 file (modify as needed). """
        snapPath = basePath + '/snapdir_' + str(snapNum).zfill(3) + '/'
        filePath = snapPath + 'snapshot_' + str(snapNum).zfill(3)
        filePath += '.' + str(chunkNum) + '.hdf5'

        return filePath

    with h5py.File(snapPath(basePath,snapNum),'r') as f:

        header = dict( f['Header'].attrs.items() )
        result['count'] = f['Header'].attrs['NumPart_Total']
    # loop over chunks
    num_perfile= np.zeros((header['NumFilesPerSnapshot'],6),dtype=np.int64)
    for i in range(header['NumFilesPerSnapshot']):
        f = h5py.File(snapPath(basePath,snapNum,i),'r')
        this_num= f['Header'].attrs['NumPart_ThisFile']
        for parttype in range(6):
            try:
                num_act= f['PartType'+str(parttype)]['ParticleIDs'].shape[0]
            except KeyError:
                num_act=0
            #print(num_act)
            if this_num[parttype]!= num_act:
                print('num in Header is wrong',i,parttype,this_num[parttype],num_act)
                return None
            
        num_perfile[i]= this_num
    result['num_perfile']= num_perfile
    return result
offset_basePath=None
def set_offset_basePath(path):
    global offset_basePath
    offset_basePath=path
    
def offsetPath(basePath, snapNum,):
    """ Return absolute path to a separate offset file (modify as needed). """
    #offsetPath = basePath + '../postprocessing/offsets/offsets_%03d.hdf5' % snapNum
    global offset_basePath
    if offset_basePath is None:
        offset_basePath= '/data/dell4/userdir/jyj/Auriga_offset/'
    offsetPath= offset_basePath+'/'.join(basePath.rstrip('/').split('/')[-3:-1])+ \
                               '/postprocessing/offsets/offsets_%03d.hdf5' % snapNum
    if os.path.exists(offsetPath):
        #print(offsetPath)
        return offsetPath
    if not os.path.exists(os.path.dirname(offsetPath)):
        print('generating offset file:',offsetPath)
        os.makedirs(os.path.dirname(offsetPath))
    # read lentype
    offgroup= loadHalos(basePath,snapNum,["GroupLenType","GroupNsubs","GroupFirstSub",],gen_off=True)
    offsub= loadSubhalos(basePath, snapNum,["SubhaloLenType","SubhaloMassType"],gen_off=True)
    #gen offset 
    GroupOffset=np.zeros((offgroup['count'],6),dtype=np.int64)
    SubOffset= np.zeros((offsub['count'],6),dtype=np.int64)
    k=0
    for i in range(0, offgroup['count']):
        if (i>0):
            GroupOffset[i] =  GroupOffset[i-1] + offgroup['GroupLenType'][i-1]
        if offgroup['GroupNsubs'][i]>0:
            SubOffset[k] = GroupOffset[i]
            k += 1
            for j in range(1, offgroup['GroupNsubs'][i]):
                SubOffset[k] =  SubOffset[k-1] + offsub['SubhaloLenType'][k-1]
                k+=1
    if k!=offsub['count']:
        print("READHALO: problem with offset table", k, offsub['count'])

    f=h5py.File(offsetPath,'w')
    f['Group/SnapByType']= GroupOffset
    f['Subhalo/SnapByType']= SubOffset
    f['FileOffsets/Group'],f['FileOffsets/Subhalo'],f['FileOffsets/SnapByType']= map(
                       lambda x:np.append(x[0:1]-x[0:1],np.cumsum(x,axis=0)[:-1],axis=0), [offgroup['num_perfile'],
                                                                     offsub['num_perfile'],
                                                       snapfileoff(basePath,snapNum)['num_perfile']])
    f.close()
    return offsetPath


def loadObjects(basePath, snapNum, gName, nName, fields, gen_off=False):
    """ Load either halo or subhalo information from the group catalog. """
    result = {}

    # make sure fields is not a single element
    if isinstance(fields, six.string_types):
        fields = [fields]

    # load header from first chunk
    with h5py.File(gcPath(basePath, snapNum), 'r') as f:

        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header'].attrs['N' + nName + '_Total']

        if not result['count']:
            print('warning: zero groups, empty return (snap=' + str(snapNum) + ').')
            return result

        # if fields not specified, load everything
        if not fields:
            fields = list(f[gName].keys())

        for field in fields:
            # verify existence
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # replace local length with global
            shape = list(f[gName][field].shape)
            shape[0] = result['count']

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=f[gName][field].dtype)

    # loop over chunks
    wOffset = 0
    
    if gen_off:
        num_perfile= np.zeros(header['NumFiles'],dtype=np.int64)
        
    for i in range(header['NumFiles']):
        f = h5py.File(gcPath(basePath, snapNum, i), 'r')

        if not f['Header'].attrs['N'+nName+'_ThisFile']:
            continue  # empty file chunk
        if gen_off:
            num_perfile[i]= f['Header'].attrs['N'+nName+'_ThisFile']
        # loop over each requested field
        for field in fields:
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # shape and type
            shape = f[gName][field].shape

            # read data local to the current file
            if len(shape) == 1:
                result[field][wOffset:wOffset+shape[0]] = f[gName][field][0:shape[0]]
            else:
                result[field][wOffset:wOffset+shape[0], :] = f[gName][field][0:shape[0], :]

        wOffset += shape[0]
        f.close()
    if gen_off:
        result['num_perfile']= num_perfile
    # only a single field? then return the array instead of a single item dict
    #if len(fields) == 1:
    #    return result[fields[0]]

    return result


def loadSubhalos(basePath, snapNum, fields=None, gen_off=False):
    """ Load all subhalo information from the entire group catalog for one snapshot
       (optionally restrict to a subset given by fields). """

    return loadObjects(basePath, snapNum, "Subhalo", "subgroups", fields, gen_off)


def loadHalos(basePath, snapNum, fields=None, gen_off=False):
    """ Load all halo information from the entire group catalog for one snapshot
       (optionally restrict to a subset given by fields). """

    return loadObjects(basePath, snapNum, "Group", "groups", fields, gen_off)


def loadHeader(basePath, snapNum):
    """ Load the group catalog header. """
    with h5py.File(gcPath(basePath, snapNum), 'r') as f:
        header = dict(f['Header'].attrs.items())

    return header


def load(basePath, snapNum):
    """ Load complete group catalog all at once. """
    r = {}
    r['subhalos'] = loadSubhalos(basePath, snapNum)
    r['halos']    = loadHalos(basePath, snapNum)
    r['header']   = loadHeader(basePath, snapNum)
    return r


def loadSingle(basePath, snapNum, haloID=-1, subhaloID=-1):
    """ Return complete group catalog information for one halo or subhalo. """
    if (haloID < 0 and subhaloID < 0) or (haloID >= 0 and subhaloID >= 0):
        raise Exception("Must specify either haloID or subhaloID (and not both).")

    gName = "Subhalo" if subhaloID >= 0 else "Group"
    searchID = subhaloID if subhaloID >= 0 else haloID

    # old or new format
    if 'fof_subhalo' in gcPath(basePath, snapNum):
        # use separate 'offsets_nnn.hdf5' files
        with h5py.File(offsetPath(basePath, snapNum), 'r') as f:
            offsets = f['FileOffsets/'+gName][()]
    else:
        # use header of group catalog
        with h5py.File(gcPath(basePath, snapNum), 'r') as f:
            offsets = f['Header'].attrs['FileOffsets_'+gName]

    offsets = searchID - offsets
    fileNum = np.max(np.where(offsets >= 0))
    groupOffset = offsets[fileNum]

    # load halo/subhalo fields into a dict
    result = {}

    with h5py.File(gcPath(basePath, snapNum, fileNum), 'r') as f:
        for haloProp in f[gName].keys():
            result[haloProp] = f[gName][haloProp][groupOffset]

    return result
