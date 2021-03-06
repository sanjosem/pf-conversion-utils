#!/usr/bin/python3
# from importlib import reload  # If on Python 3

import yaml
import numpy as np
import sys
import os.path
from os import makedirs
from glob import glob

import pftools.pnc_conv

if len(sys.argv)>1:
    input_file = sys.argv[1]
else:
    input_file = "../input_post_probes.yaml"

if not os.path.exists(input_file):
    print('Usage: read_temporal.py [input_post_probes.yaml]')
    sys.exit('Input file \'{0:s}\' cannot be found'.format(input_file))

fio = open(input_file,'r')
uinfo = yaml.load(fio,Loader=yaml.Loader)
fio.close()

for required_field in ['directory','probe_name','file_type']:
    if required_field not in uinfo.keys():
        sys.exit('Missing required parameter {0:s} in \'{1:s}\''.format(required_field,input_file))
if 'output_directory' not in uinfo.keys():
    uinfo['output_directory'] = uinfo['directory']
if 'delimiter' not in uinfo.keys():
    uinfo['delimiter'] = ' '
        
if not os.path.isdir(uinfo['directory']):
    sys.exit('Directory {0:s} cannot be found'.format(uinfo['directory']))
if not os.path.isdir(uinfo['output_directory']):
    makedirs(uinfo['output_directory'])
    
liste = glob('{0:s}/{1:s}*.{2:s}'.format(uinfo['directory'],uinfo['probe_name'],uinfo['file_type']))

if len(liste) == 0:
    sys.exit('No probe file \'{1:s}*.{2:s}\' found in directory {0:s}'.format(uinfo['directory'],uinfo['probe_name'],uinfo['file_type']))

for fl in liste:
    probe_file = os.path.basename(fl)

    print('*** Processing {0:s} ***'.format(probe_file))
    
    if 'verbose' in uinfo.keys():
        pncConv = pftools.pnc_conv.pncConversion(fl,uinfo['verbose'])
    else:
        pncConv = pftools.pnc_conv.pncConversion(fl)
                
    pncConv.export_temporal_data(probe_file.replace('.psnc','').replace('.pfnc',''),
                                    uinfo['output_directory'],delimiter=uinfo['delimiter'])
    
