#!/usr/bin/python3
from importlib import reload  # If on Python 3

import yaml
import numpy as np
import sys
import os.path
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
    makedirs(uinfo['output_directory'])

sncfile = os.path.join(uinfo['directory'],uinfo['snc_filename'])

if not os.path.exists(sncfile):
    sys.exit('File {0:s} cannot be found'.format(sncfile))
if 'casename' not in uinfo.keys():
    uinfo['casename'] = uinfo['snc_filename'].replace('.snc','')

if 'verbose' in uinfo.keys():
    sncConv = pfsnct.sncConversion(sncfile,uinfo['verbose'])
else:
    sncConv = pfsnct.sncConversion(sncfile)
    
sncConv.read_conversion_parameters()

sncConv.read_surface_names()
list_av_face_names = sncConv.surface_list['names']

if uinfo['face_names'] is not None:
    for name in uinfo['face_names'].keys():
        for ptn in uinfo['face_names'][name]:
            if not ptn in list_av_face_names:
                print('Face {0:s} is not available in {1:s}'.format(ptn,uinfo['snc_filename']))
else:
    uinfo['face_names'] = list_av_face_names

sncConv.read_coordinates()
sncConv.read_connectivity()

if 'surface_min' in uinfo.keys():
    for surface_name in uinfo['face_names'].keys():
        sncConv.triangulate_surface(surface_name,uinfo['face_names'][surface_name],
            uinfo['surface_min'])
else:
    for surface_name in uinfo['face_names'].keys():
        sncConv.triangulate_surface(surface_name,uinfo['face_names'][surface_name])

for surface_name in uinfo['face_names'].keys():
    sncConv.read_frame_data(surface_name,0)

sncConv.save_vtk(uinfo['casename'],uinfo['output_directory'])
