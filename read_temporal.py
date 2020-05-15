#!/usr/bin/python3
from importlib import reload  # If on Python 3

import yaml
import numpy as np
import sys
import os.path
from os import makedirs
from glob import glob

import pftools # as pfpnct
reload(pftools)
import pftools # as pfpnct
reload(pftools)

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
        
    pncConv.read_conversion_parameters()
    pncConv.define_measurement_variables()
    
    toto

    f1 = netcdf.netcdf_file('{0:s}/{1:s}'.format(uinfo['directory'],probe_file), 'r')
    
    meas1 = f1.variables['measurements']
    st1 = f1.variables['start_time'][()]
    et1 = f1.variables['end_time'][()]

    lx_scales = f1.variables['lx_scales']
    lx_offsets = f1.variables['lx_offsets']
    
    if uinfo['file_type'] == 'pfnc':
        weight = f1.variables['fluid_volumes']
    elif uinfo['file_type'] == 'psnc':
        weight = f1.variables['surfel_area']

    if sys.version_info[0] < 3:
        vnames = f1.variables['variable_short_names']
        variable_names = deepcopy(vnames[0: -1]) # copy without last 0
        variable_names[variable_names == ''] = ' ' # replace 0 by space
        var_list = ''.join(variable_names).split(' ') # join and split
    else:
        x = f1.variables['variable_short_names'][()]
        var_list = x.tostring().decode('utf-8').split('\x00')[:-1]

    print("Variables found:",var_list)

    # Parse var name to find index
    ipress=-1
    irho=-1
    ivx=-1
    for iv,var in enumerate(var_list):
        if var=='x_velocity':
            ivx=iv
        elif var=='y_velocity':
            ivy=iv
        elif var=='z_velocity':
            ivz=iv
        elif var=='static_pressure':
            ipress=iv
        elif var=='density':
            irho=iv
    var_list=[]
    var_idx=[]
    if ipress==-1:
        var_list.append('density')
        var_idx.append(irho)
    else:
        var_list.append('static_pressure')
        var_idx.append(ipress)
    if ivx>=0:
        var_list.append('x_velocity')
        var_list.append('y_velocity')
        var_list.append('z_velocity')
        var_idx.append(ivx)
        var_idx.append(ivy)
        var_idx.append(ivz)

    coeff_dx = lx_scales[0]
    CFL_number = lx_scales[-1]
    coeff_dt = lx_scales[0] / lx_scales[3]
    dt = CFL_number * coeff_dt
    coeff_press = lx_scales[1] * lx_scales[3]**2
    coeff_vel = lx_scales[3]
    coeff_density = lx_scales[1]
    offset_pressure = lx_offsets[4]

    t_lat = 1.0/3.0 # lattive temperature
    r_lat = 1.0 # lattice r_gas constant
    rTlat = r_lat*t_lat # weight for density to pressure
    irTlat = 1.0/rTlat # weight for pressure to density

    # Average point
    iscale = 1.0/float(np.sum(weight[:]))

    mean_meas=np.sum(meas1[()]*weight[:],axis=2)*iscale

    data = OrderedDict()

    data['time'] = 0.5*(st1+et1)*dt

    if ipress==-1:
        data[var_list[0]]=mean_meas[:,var_idx[0]] * coeff_density
        data['static_pressure']=(mean_meas[:,var_idx[0]] * rTlat + offset_pressure) * coeff_press
    else:
        data[var_list[0]]=(mean_meas[:,var_idx[0]]+offset_pressure)*coeff_press
        data['density']=mean_meas[:,var_idx[0]] * irTlat * coeff_density

    if ivx>=0:
        data[var_list[1]]=mean_meas[:,var_idx[1]] * coeff_vel
        data[var_list[2]]=mean_meas[:,var_idx[2]] * coeff_vel
        data[var_list[3]]=mean_meas[:,var_idx[3]] * coeff_vel

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f1.close()

    export_file = '{1:s}/{0:s}.txt'.format(probe_file.replace('.psnc','').replace('.pfnc',''),uinfo['output_directory'])
    delim = uinfo['delimiter']
    print("Exporting data into {0:s}".format(export_file))
    header = delim.join(data.keys())
    nf = len(data.keys())
    nt = data['time'].size
    X = np.zeros((nt,nf))
    for iff,fn in enumerate(data.keys()):
        X[:,iff] = data[fn]

    np.savetxt(export_file,X,header=header,comments='',delimiter=delim)
