from antares import * 
import h5py
import sys
import os.path

if len(sys.argv)==2:
    filename = sys.argv[1]
else: 
    print('Usage: convert_to_antares.py <params_pf_casename.hdf5>')
    raise RuntimeError('Missing filename to convert')

casename = os.path.basename(filename).replace('params_pf','').replace('.hdf5','')

fparams = h5py.File(filename,'r')
if not fparams.get('surface_name',getclass=True) is None:
    if fparams.get('surface_name',getclass=True) == h5py.Group:
        mesh_group = fparams['surface_name']
        zones = list(mesh_group.keys())
    else:
        raise RuntimeError('Wrong file format')
else:
    raise RuntimeError('Wrong file format, is it a surface file?')

print('Surface meshes to read:',zones)
read_data = False
if not fparams.get('data',getclass=True) is None:
    if fparams.get('data',getclass=True) == h5py.Group:
        data_group = fparams['data']
        instants = list(data_group[zones[0]].keys())
        read_data = True
    else:
        print('No instant to read')
else:
    print('No instant to read')

if read_data:
    print('Frames to read:',instants)
else:
    instants = ['0000',]

b = Base()
b.init(zones=zones,instants=instants)
for zn in zones:
    gsurf = mesh_group[zn]
    coords = gsurf['coordinates']
    for icoor,coor in enumerate(['x','y','z']):
        b[zn].shared[coor] = coords[:,icoor]
    b[zn].shared.connectivity['tri'] = gsurf['tri_conn'][()]
    b[zn].shared.connectivity['qua'] = gsurf['qua_conn'][()]
    for inst in instants:
        b[zn][inst]['pressure'] = data_group[f'{zn}/{inst}/static_pressure'][()]

w = Writer('hdf_antares')
w['base'] = b
w['filename'] = f'Antares_{casename}'
w.dump()