'''Class to store parameters and functions for the convertion of PowerFLOW files

Author: Marlene Sanjose
'''

class PFConversion:
    """Initiate the convertion class associated to a given PowerFLOW file.

    Parameters
    ----------
    pfFile : string
        File name (and path) of the PowerFLOW file to convert
    verbose : bool
        Activate or desactivate informative prints (default: True)

    Methods:
    ----------
    read_conversion_parameters
        Read and compute conversion parameters
    read_rotation_information
        Read rotating domain parameters (if existing)
    define_measurement_variables
        Read available variables
    extract_time_info
        Read time information such as sampling and averaging time
    save_parameters
        Save conversion parameters
    """
    def __init__(self,pfFile,verbose=True):

        from os.path import exists as fexists
        if not fexists(pfFile):
            raise RuntimeError('File {0:s} could not be found'.format(pfFile))

        self.pfFile = pfFile
        self.format = None
        self.params = None
        self.rotation = False
        self.vars = None
        self.time = None
        self.verbose = verbose
        self.probe = None
        self.fapi = None

        if self.verbose:
            print('PFConversion class instance created for file:\n  {0:s}'.format(pfFile))


    def read_conversion_parameters(self):
        """Function to read conversion factor, parameters are stored in params dictionnary in the class

        """

        import netCDF4 as netcdf

        f = netcdf.Dataset(self.pfFile, 'r')

        lx_offsets = f.variables['lx_offsets'][()]
        lx_scales = f.variables['lx_scales'][()]
        param_values = f.variables['param_values'][()]

        self.params = dict()
        self.params['coeff_dx'] = lx_scales[0]
        self.params['CFL_number'] = lx_scales[-1]
        self.params['coeff_dt'] = lx_scales[0] / lx_scales[3]
        self.params['dt'] = self.params['CFL_number'] * self.params['coeff_dt']
        self.params['coeff_press'] = lx_scales[1] * lx_scales[3]**2
        self.params['coeff_force'] = self.params['coeff_press']
        self.params['coeff_vel'] = lx_scales[3]
        self.params['coeff_massflux'] = lx_scales[1] * lx_scales[3]
        self.params['coeff_density'] = lx_scales[1]
        self.params['offset_pressure'] = lx_offsets[4]

        t_lat = 1.0/3.0 # lattive temperature
        r_lat = 1.0 # lattice r_gas constant
        self.params['weight_rho_to_pressure'] = r_lat*t_lat # weight for density to pressure
        self.params['weight_pressure_to_rho'] = 1.0/(r_lat*t_lat) # weight for pressure to density

        self.params['p_ref_lat'] = param_values[-1]
        self.params['rho_ref_lat'] = param_values[3]
        self.params['U_ref_lat'] = param_values[1]

        if self.verbose:
            print('  CFL number {0:.3f}'.format(self.params['CFL_number']))
            print('  Time step {0:.4e} s'.format(self.params['dt']))
            print('  Minimum cell size {0:.4e} m'.format(self.params['coeff_dx']))

        # Spatial offset
        self.params['offset_coords'] = f.variables['csys'][0, 0:3, 3].astype('float')

        f.close()

        self.read_rotation_information()

    def read_rotation_information(self):
        """Function to read rotating domain information, if present

        """

        import netCDF4 as netcdf
        from numpy import where,pi

        f = netcdf.Dataset(self.pfFile, 'r')

        if 'ref_frame_indices' in f.variables.keys():
            nlrfs = f.dimensions['nlrfs'].size
            if nlrfs == 1:
                self.rotation=True
                if self.verbose:
                    print('Rotating domain (assuming 1) present in the computational domain')
            else:
                raise NotImplementedError('Class only support single rotating domain')

            tmp = f.variables['lrf_constant_angular_vel_mag'][0]
            self.params['omega'] = tmp / self.params['coeff_dt']

            if 'lrf_initial_n_revolutions' in f.variables.keys():
                n_rot0 = f.variables['lrf_initial_n_revolutions'][0] * pi * 2.0
            else:
                n_rot0 = 0.0
            alpha0 = f.variables['lrf_initial_angular_rotation'][0] + n_rot0

            rotation_axis_origin = f.variables['lrf_axis_origin'][0,:]  # axe de rotation
            rotation_axis = f.variables['lrf_axis_direction'][0,:]  # axe de rotation
            axis_index = where(rotation_axis)[0][0]
            sign_rotation = rotation_axis[axis_index]
            if abs(sign_rotation)<1:
                print('Rotation axis',rotation_axis)
                raise NotImplementedError('Class only support principal axis for rotating')
            self.params['axis'] = 'xyz'[axis_index]
            self.params['iaxis'] = axis_index

            self.params['sign_rotation'] = sign_rotation
            self.params['init_angle'] = alpha0*sign_rotation

            if self.verbose:
                print("Rotating speed {0:.3f} rad/s".format(self.params['omega']))
                print("Initial mesh position {0:.3f} rad".format(self.params['init_angle']))

        else:
            self.rotation=False
            if self.verbose:
                print('No rotating domain in the computational domain')

        f.close()


    def define_measurement_variables(self):
        """Extract the available variable names and store in vars dictionary

        """


        import netCDF4 as netcdf
        from numpy import where

        f = netcdf.Dataset(self.pfFile, 'r')

        if self.verbose:
            print("Extracting Variable names")

        x = f.variables['variable_short_names'][()]
        nvars = f.dimensions['n_variables'].size
        var_list = x.tostring().decode('utf-8').split('\x00')[:nvars]
        if self.verbose:
            print("  -> Found variables:",var_list)

        self.vars = dict()

        # Parse var name to find index
        for iv,var in enumerate(var_list):
            if var in ['x_velocity','y_velocity','z_velocity','static_pressure',
                       'density','surface_x_force','surface_y_force','surface_z_force',
                        'mass_flux']:
                self.vars[var] = iv
                if var == 'static_pressure':
                    self.vars['density'] = -1
                if var == 'density':
                    self.vars['static_pressure'] = -1
            else:
                print('  ! Conversion for variable {0:s} is not yet implemented'.format(var))

        f.close()

    def extract_time_info(self):
        """Extract the time information and store the parameters in the time dictionnary in the class instance.

        """

        import netCDF4 as netcdf
        from numpy import diff
        if self.fapi:
            if self.format == 'volume':
                import pftools.fextend.fnc_reader as Fapi
            elif self.format == 'surface':
                import pftools.fextend.snc_reader as Fapi
            else:
                self.fapi=False

        if self.params is None:
            self.read_conversion_parameters()

        if self.verbose:
            print("Extracting Time information")

        self.time = dict()

        f = netcdf.Dataset(self.pfFile, 'r')
        nsets = f.dimensions['nsets'].size
        self.time['nsets'] = nsets
        f.close()

        if self.verbose:
            print('  -> Number of time frames: {0:d}'.format(self.time['nsets']))

        if self.fapi and nsets>10:
            self.time['Tavg'],self.time['Tsampling'], \
            self.time['time_center'] = Fapi.pf_params.read_time(self.pfFile,
                                                            self.params['dt'],
                                                            self.time['nsets'])
            self.time['Fsampling'] = 1.0/self.time['Tsampling']
            print('  -> Frequency sampling: {0:,.0f} Hz'.format(self.time['Fsampling']))
        else:

            f = netcdf.Dataset(self.pfFile, 'r')

            st = f.variables['start_time'][()]
            et = f.variables['end_time'][()]
            self.time['time_center'] = 0.5*(st+et)*self.params['dt']

            if self.verbose:
                print('  -> First frame ({0:d}): {1:g} s'.format(0,self.time['time_center'][0]))
                if self.time['nsets']>1:
                    print('  -> Last frame ({0:d}): {1:g} s'.format(self.time['nsets'],self.time['time_center'][-1]))

            avg_period = et - st

            if not avg_period.min() == avg_period.max():
                print('  ! Averaging period is not constant')
                print('  ! Min/max sampling [timestep]: {0:.0f} {1:.0f}'.format(avg_period.min(),
                                    avg_period.max()))
            else:
                self.time['Tavg'] = avg_period.mean()*self.params['dt']
                if self.verbose:
                    print('  -> Time averaging: {0:.6e} s'.format(self.time['Tavg']))

            if nsets>1:
                sampling_period = diff(self.time['time_center'])

                if not (sampling_period.min()-sampling_period.max())<self.params['dt']:
                    print('  ! Sampling period is not constant')
                    print('  ! Min/max sampling [sec]: {0:.6e} {1:.6e}'.format(sampling_period.min(),
                            sampling_period.max()))
                else:
                    self.time['Tsampling'] = sampling_period.mean()
                    self.time['Fsampling'] = 1.0/self.time['Tsampling']
                    if self.verbose:
                        print('  -> Time sampling: {0:.6e} s'.format(self.time['Tsampling']))
                        print('  -> Frequency sampling: {0:,.0f} Hz'.format(self.time['Fsampling']))

            f.close()

    def save_parameters(self,casename,dirout):
        """Function to export convertion parameters in a separated hdf5 file.

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        import h5py
        import os.path
        from numpy import array


        outFile = os.path.join(dirout,'param_pf_{0:s}.hdf5'.format(casename))
        print("Creating a new convertion parameter file:\n  ->  {0:s}".format(outFile))

        if self.params is not None:
          param_list = self.params.keys()

        fparams = h5py.File(outFile,'w')

        fparams.create_dataset('filename',data=self.pfFile)
        fparams.create_dataset('format',data=self.format)

        grp = fparams.create_group("constants")

        for par in param_list:
            if par not in ['omega','init_angle']:
                grp.create_dataset(par, data=self.params[par])

        if self.rotation:
            grot = grp.create_group("rotation")
            for par in ['omega','init_angle']:
                grot.create_dataset(par, data=self.params[par])

        if self.vars is not None:
            gvar = fparams.create_group("variables")
            varnames = []
            varidx = []
            for var in self.vars.keys():
                varnames.append(var)
                varidx.append(self.vars[var])
            gvar.create_dataset('names',data=array(varnames,dtype='S'))
            gvar.create_dataset('indices',data=array(varidx))

        if self.time is not None:
            gtime = fparams.create_group("time")
            for par in self.time.keys():
                gtime.create_dataset(par,data=self.time[par])

        fparams.close()
        return outFile

    def load_parameters(self,h5file):

        import h5py

        fparams = h5py.File(h5file,'r')

        if not self.pfFile == fparams['filename'][()]:
            RuntimeError('Parameters do not correspond to present analysed file {0:s}'.format(self.pfFile))
        if not self.format == fparams['format'][()]:
            RuntimeError('Parameters do not correspond to present analysed file format {0:s}'.format(self.format))

        print("Loading the convertion parameter file:\n  ->  {0:s}".format(h5file))

        if fparams.get('constants',getclass=True) is None:
            RuntimeError('No constants group file, {0:s} is an invalid convertion parameter file'.format(self.format))
        else:
            if not fparams.get('constants',getclass=True) == h5py.Group:
                RuntimeError('No constants group file, {0:s} is an invalid convertion parameter file'.format(self.format))

        grp = fparams["constants"]
        self.params = dict()
        for par in list(grp.keys()):
            if grp.get(par,getclass=True) == h5py.Dataset:
                self.params[par] = grp[par][()]
            elif par == 'rotation':
                self.rotation = True
                grot = fparams["constants/rotation"]
                for key in list(grot.keys()):
                    if grot.get(key,getclass=True) == h5py.Dataset:
                        self.params[key] = grot[key][()]

        if not fparams.get('variables',getclass=True) is None:
            if fparams.get('variables',getclass=True) == h5py.Group:

                grp = fparams["variables"]
                varnames = grp['names'][()]
                varidx = grp['indices'][()]
                self.vars = dict()

                for name, idx in zip(varnames,varidx):
                    self.vars[name.decode('utf-8')] = idx

        if not fparams.get('time',getclass=True) is None:
            if fparams.get('time',getclass=True) == h5py.Group:
                grp = fparams["time"]
                self.time = dict()
                for par in list(grp.keys()):
                    self.time[par] = grp[par][()]

        fparams.close()


    def export_temporal_data(self,casename,dirout,delimiter=' ',index=False,
                                 extension='txt'):
        """Function to export probe temporal data to a text file.
        All quantities will be written in SI units

        Parameters
        ----------
        casename : string
            Name assiciated to the present probe conversion.
            The casename is used to build the output file name
            temporal_<casename>.<extension>
        dirout : string
            Absolute path/relative path where the converted file will be written
        delimiter : char
            Field delimiter (default is space).
            If comma ',' is specified the file extension will be 'csv'
        index : bool
            Append index of rows as the first column in the text file
        extension : string
            Extension of the text file (by default txt)


        """
        import os.path

        if delimiter == ',':
            ext = 'csv'
        else:
            ext = extension

        if self.probe is None:
            raise RuntimeError('No probe to export, please use extract_probe')

        for probe_name in self.probe.keys():
            outFile = os.path.join(dirout,'temporal_{0:s}_{2:s}.{1:s}'.format(
                        casename,ext,probe_name))
            if extension == 'pkl':
                print("Exporting in pickle format:\n  ->  {0:s}".format(outFile))
                self.probe[probe_name].to_pickle(outFile)
            else:
                print("Exporting in ascii column format:\n  ->  {0:s}".format(outFile))
                self.probe[probe_name].to_csv(outFile,sep=delimiter,index=index)
