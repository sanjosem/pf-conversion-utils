#!/usr/bin/env python
from importlib import reload  # If on Python 3

import yaml
import numpy as np
import sys
import os.path
import pandas as pd
from os import makedirs

import pftools.snc_conv as pfsnct
reload(pfsnct)

if len(sys.argv)>1:
    input_file = sys.argv[1]
else:
    input_file = "../input_post_surface.yaml"

if not os.path.exists(input_file):
    print('Usage: script_clean_surface.py [input_post_surface.yaml]')
    sys.exit('Input file \'{0:s}\' cannot be found'.format(input_file))

fio = open(input_file,'r')
uinfo = yaml.load(fio,Loader=yaml.Loader)
fio.close()

for required_field in ['directory','snc_filename']:
    if required_field not in uinfo.keys():
        sys.exit('Missing required parameter {0:s} in \'{1:s}\''.format(required_field,input_file))
if 'output_directory' not in uinfo.keys():
    uinfo['output_directory'] = uinfo['directory']
if 'face_names' not in uinfo.keys():
    uinfo['face_names'] = None

if not os.path.isdir(uinfo['directory']):
    sys.exit('Directory {0:s} cannot be found'.format(uinfo['directory']))
if not os.path.isdir(uinfo['output_directory']):
    makedirs(uinfo['output_directory'],exist_ok=True)

sncfile = os.path.join(uinfo['directory'],uinfo['snc_filename'])

if not os.path.exists(sncfile):
    sys.exit('File {0:s} cannot be found'.format(sncfile))
if 'casename' not in uinfo.keys():
    uinfo['casename'] = uinfo['snc_filename'].replace('.snc','')

if 'verbose' in uinfo.keys():
    sncConv = pfsnct.sncConversion(sncfile,uinfo['verbose'])
else:
    sncConv = pfsnct.sncConversion(sncfile)

if 'save_format' in uinfo.keys():
    save_format = uinfo['save_format']
else:
    save_format = 'vtk'

sncConv.read_conversion_parameters()

sncConv.read_surface_names()
list_av_face_names = sncConv.surface_list['names']

if uinfo['face_names'] is not None:
    for name in uinfo['face_names'].keys():
        for ptn in uinfo['face_names'][name]:
            if not ptn in list_av_face_names:
                print('Face {0:s} is not available in {1:s}'.format(ptn,uinfo['snc_filename']))
else:
    uinfo['face_names'] = dict()
    for name in list_av_face_names:
        uinfo['face_names'][name] = [name,]

if 'frame' in uinfo.keys():
    frame_list = uinfo['frame']
elif 'framelist' in uinfo.keys():
    fl = uinfo['framelist']
    frame_list = np.arange(fl['start'],fl['end'],fl['step'],dtype='i').tolist()
else:
    frame_list = [0,]
print(frame_list)

sncConv.read_coordinates()
sncConv.read_connectivity()

# Probes
if 'probes' in uinfo.keys():
    for probe_name in uinfo['probes'].keys():
        sncConv.extract_probe(probe_name,uinfo['probes'][probe_name])
    sncConv.export_temporal_data(uinfo['casename'],uinfo['output_directory'],extension='txt')
elif 'probes_file'  in uinfo.keys():
    ldf = pd.read_excel(uinfo['probes_file'],sheet_name=None)
    if 'probes_file_sheet' in uinfo.keys():
        work_list = (uinfo['probes_file_sheet'],)
    else:
        work_list = ldf.keys()
    for pbcol in work_list:
        for ip in range(ldf[pbcol].shape[0]):
            probe_name = '{0:s}-{1:s}'.format(pbcol,ldf[pbcol].loc[ip,'probe_name'])
            print('extracting probe {0:s}'.format(probe_name))
            sncConv.extract_probe(probe_name,(ldf[pbcol].loc[ip,'x'],ldf[pbcol].loc[ip,'y'],ldf[pbcol].loc[ip,'z']))
    # sncConv.export_temporal_data(uinfo['casename'],uinfo['output_directory'],delimiter=',')
    sncConv.export_temporal_data(uinfo['casename'],uinfo['output_directory'],extension='pkl')

else:

    if 'surface_min' in uinfo.keys():
        for surface_name in uinfo['face_names'].keys():
            sncConv.triangulate_surface(surface_name,uinfo['face_names'][surface_name],
                uinfo['surface_min'])
    else:
        for surface_name in uinfo['face_names'].keys():
            sncConv.triangulate_surface(surface_name,uinfo['face_names'][surface_name])

    if save_format == 'hdf5':
        outFile = sncConv.save_parameters(uinfo['casename'],uinfo['output_directory'])

    for frame_number in frame_list:
        print(f' ** Reading frame {frame_number}')
        for surface_name in uinfo['face_names'].keys():
            sncConv.read_frame_data(surface_name,frame_number)

            if save_format == 'vtk':
                sncConv.save_vtk(f"{uinfo['casename']}_{frame_number:04d}",uinfo['output_directory'])
            elif save_format == 'hdf5':                
                sncConv.save_instant(outFile,surface_name,frame_number)
