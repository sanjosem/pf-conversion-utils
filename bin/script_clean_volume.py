#!/usr/bin/env python
from importlib import reload  # If on Python 3

import yaml
import numpy as np
import sys
import os.path
from os import makedirs

import pftools.fnc_conv as pft
reload(pft)

if len(sys.argv)>1:
    input_file = sys.argv[1]
else:
    input_file = "../input_post_volume.yaml"

if not os.path.exists(input_file):
    print('Usage: script_clean_volume.py [input_post_volume.yaml]')
    sys.exit('Input file \'{0:s}\' cannot be found'.format(input_file))

fio = open(input_file,'r')
uinfo = yaml.load(fio,Loader=yaml.Loader)
fio.close()

for required_field in ['directory','fnc_filename']:
    if required_field not in uinfo.keys():
        sys.exit('Missing required parameter {0:s} in \'{1:s}\''.format(required_field,input_file))
if 'output_directory' not in uinfo.keys():
    uinfo['output_directory'] = uinfo['directory']
    
if not os.path.isdir(uinfo['directory']):
    sys.exit('Directory {0:s} cannot be found'.format(uinfo['directory']))
if not os.path.isdir(uinfo['output_directory']):
    makedirs(uinfo['output_directory'])

fncfile = os.path.join(uinfo['directory'],uinfo['fnc_filename'])

if not os.path.exists(fncfile):
    sys.exit('File {0:s} cannot be found'.format(fncfile))
if 'casename' not in uinfo.keys():
    uinfo['casename'] = uinfo['fnc_filename'].replace('.fnc','')

if 'verbose' in uinfo.keys():
    fncConv = pft.fncConversion(fncfile,uinfo['verbose'])
else:
    fncConv = pft.fncConversion(fncfile)
    
fncConv.read_conversion_parameters()

fncConv.read_volume_mesh()

fncConv.read_frame_data(0)

fncConv.save_vtk(uinfo['casename'],uinfo['output_directory'])
