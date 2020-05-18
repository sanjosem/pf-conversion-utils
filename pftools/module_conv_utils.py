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
        
        if self.verbose:
            print('PFConversion class instance created for file:\n  {0:s}'.format(pfFile))
        

    def read_conversion_parameters(self):
        """Function to read conversion factor, parameters are stored in params dictionnary in the class

        """
        
        from scipy.io import netcdf
        
        f = netcdf.netcdf_file(self.pfFile, 'r', mmap=False)
        
        lx_offsets = f.variables['lx_offsets'][()]
        lx_scales = f.variables['lx_scales'][()]
        param_values = f.variables['param_values'][()]

        self.params = dict()
        self.params['coeff_dx'] = lx_scales[0]
        self.params['CFL_number'] = lx_scales[-1]
        self.params['coeff_dt'] = lx_scales[0] / lx_scales[3]
        self.params['dt'] = self.params['CFL_number'] * self.params['coeff_dt']
        self.params['coeff_press'] = lx_scales[1] * lx_scales[3]**2
        self.params['coeff_vel'] = lx_scales[3]
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
        
        from scipy.io import netcdf
        from numpy import where,pi
        
        f = netcdf.netcdf_file(self.pfFile, 'r', mmap=False)
            
        if 'ref_frame_indices' in f.variables.keys():
            nlrfs = f.dimensions['nlrfs']
            if nlrfs == 1:
                self.rotation=True
                if self.verbose:
                    print('Rotating domain (assuming 1) present in the computational domain')
            else:
                raise NotImplementedError('Class only support single rotating domain')
                
            tmp = f.variables['lrf_constant_angular_vel_mag'][0]
            self.params['omega'] = tmp / self.params['coeff_dt']
            
            n_rot0 = f.variables['lrf_initial_n_revolutions'][0] * pi * 2.0
            alpha0 = f.variables['lrf_initial_angular_rotation'][0] + n_rot0

            rotation_axis_origin = f.variables['lrf_axis_origin'][0,:]  # axe de rotation
            print(rotation_axis_origin - self.params['offset_coords'])
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
                print("Rotation {0:.3f} rad".format(self.params['init_angle']))
                print("Initial mesh position {0:.3f} rad".format(self.params['init_angle']))
            
        else:
            self.rotation=False
            if self.verbose:
                print('No rotating domain in the computational domain')
                
        f.close()
        
            
    def define_measurement_variables(self):
        """Extract the available variable names and store in vars dictionary

        """
        
        
        from scipy.io import netcdf
        from numpy import where
        
        f = netcdf.netcdf_file(self.pfFile, 'r', mmap=False)
        
        if self.verbose:
            print("Extracting Variable names")
                        
        x = f.variables['variable_short_names'][()]
        var_list = x.tostring().decode('utf-8').split('\x00')[:f.dimensions['n_variables']]
        if self.verbose:
            print("  -> Found variables:",var_list)
            
        self.vars = dict()

        # Parse var name to find index
        for iv,var in enumerate(var_list):
            if var in ['x_velocity','y_velocity','z_velocity','static_pressure','density']:
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
        
        from scipy.io import netcdf
        from numpy import diff

        if self.params is None:
            self.read_conversion_parameters()

        f = netcdf.netcdf_file(self.pfFile, 'r', mmap=False)
        
        if self.verbose:
            print("Extracting Time information")
            
        self.time = dict()
        
        nsets = f.variables['start_time'].shape[0]
        self.time['nsets'] = nsets
        
        if self.verbose:
            print('  -> Number of time frames: {0:d}'.format(self.time['nsets']))
            
        st = f.variables['start_time'][()]
        et = f.variables['end_time'][()]
        self.time['time_center'] = 0.5*(st+et)*self.params['dt']
        
        if self.verbose:
            print('  -> First frame ({0:d}): {1:g} s'.format(0,self.time['time_center'][0]))
            if self.time['nsets']>1:
                print('  -> Last frame ({0:d}): {1:g} s'.format(self.time['nsets'],self.time['time_center'][-1]))
                
        avg_period = et - st
        sampling_period = diff(self.time['time_center'])
        
        if not avg_period.min() == avg_period.max():
            print('  ! Averaging period is not constant')
            print('  ! Min/max sampling [timestep]: {0:.0f} {1:.0f}'.format(avg_period.min(),
                                avg_period.max()))
        else:
            self.time['Tavg'] = avg_period.mean()*self.params['dt']
            if self.verbose:
                print('  -> Time averaging: {0:.6e} s'.format(self.time['Tavg']))

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
        
        
        outFile = os.path.join(dirout,'param_pf_{0:s}.hdf5'.format(casename))
        print("Creating a new convertion parameter file:\n  ->  {0:s}".format(outFile))
        
        if self.params is not None:
          param_list = self.params.keys()

        fparams = h5py.File(outFile,'w')

        grp = fparams.create_group("constants")
        
        for par in param_list:
            if par not in ['omega','init_angle']:
                grp.create_dataset(par, data=self.params[par])

        if self.rotation:
            grot = grp.create_group("rotation")
            for par in ['omega','init_angle']:
                grot.create_dataset(par, data=self.params[par])
                                
        fparams.close()
        
