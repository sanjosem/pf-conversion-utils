from pftools.module_conv_utils import PFConversion

class pncConversion(PFConversion):
    """Class for probe data conversion
    inherits from PFConversion

    Parameters
    ----------
    pfFile : string 
        probe file (surface or volume)
    verbose : bool
        Activate or desactivate informative prints (default: True)
        
    Methods
    -------
    compute_probe_weight
        Collect volume/surface scaling
    read_measurement
        Convert probe measurement to pandas data frame
    export_temporal_data
        Write column text file 

    """
    def __init__(self,pfFile,verbose=True):
        super().__init__(pfFile,verbose)
                    
        ext = pfFile.split('.')[-1]
        if ext == 'psnc':
            self.format = 'surface-probe'
        elif ext == 'pfnc':
            self.format = 'volume-probe'
        else:
            raise RuntimeError('This is not a probe file')
        self.iscale = None
        self.weight = None
        
        if self.verbose:
            print('PF file format is: {0:s}'.format(self.format))
        
    def get_probe_location(self):
        """Collect probe location"""

        if self.params is None:
            self.read_conversion_parameters()

        if self.format == 'volume-probe':
            f = netcdf.Dataset(self.pfFile, 'r')

            self.params['ncells'] = f.dimensions['npoints'].size
            self.params['ndims'] = f.dimensions['ndims'].size

            element_coords = f.variables['coords'][()].astype('float')

            f.close()

            for idim in range(self.params['ndims']):
                element_coords[:,idim]+=self.params['offset_coords'][idim]
            element_coords *=  self.params['coeff_dx']
            
            # if self.verbose:
            print(f'Probe elements: {self.params['ncells']}')
            print('Probe location:')
            print(element_coords.mean(axis=0))



    def compute_probe_weight(self):
        """Collect volume/surface scaling and store it in the class instance
        The results is stored in the class instance 
        and there is no input except from the class instance.


        """
        
        import netCDF4 as netcdf
        from numpy import pi
        
        if self.params is None:
            self.read_conversion_parameters()

        f = netcdf.Dataset(self.pfFile, 'r')
        
        if self.format == 'volume-probe':
            self.weight = f.variables['fluid_volumes'][()] * self.params['coeff_dx']**3
        elif self.format == 'surface-probe':
            self.weight = f.variables['surfel_area'][()] * self.params['coeff_dx']**2

        # Average point
        intv = self.weight.sum()
        self.iscale = 1.0/float(intv)
        
        f.close()

        if self.verbose:
            print('Probe size')
            if self.format == 'volume-probe':
                print('  -> Volume: {0:e} m3'.format(intv))
                rad = (3*intv/(4*pi))**(1./3.)
                print('  -> Radius: {0:e} m'.format(rad))
            if self.format == 'surface-probe':
                print('  -> Area: {0:e} m2'.format(intv))
                rad = (intv/pi)**0.5
                print('  -> Radius: {0:e} m'.format(rad))

    def extract_probe(self):
        """Function that read and convert data as probe data in SI units
        The results is stored in the class instance 
        and there is no input except from the class instance.

        """
        
        import netCDF4 as netcdf
        from pandas import DataFrame
                
        if self.iscale is None:
            self.compute_probe_weight()
        
        if self.vars is None:
            self.define_measurement_variables()
            
        if self.time is None:
            self.extract_time_info()
        
        data = dict()
        data['time'] = self.time['time_center']
        
        f = netcdf.Dataset(self.pfFile, 'r')
        meas = f.variables['measurements'][()] * self.weight
        mean_meas = meas.sum(axis=-1) * self.iscale
        f.close()
        
        for var in self.vars.keys():
            
            idx = self.vars[var]
            if var == 'static_pressure':
                if idx>=0 :
                    data[var] = ( ( mean_meas[:,idx] + self.params['offset_pressure'] ) 
                                * self.params['coeff_press'] )
                else:
                    idx = self.vars['density']
                    data[var] =  ( mean_meas[:,idx] * self.params['weight_rho_to_pressure']
                                + self.params['offset_pressure'] ) * self.params['coeff_press']
            if var == 'density':
                if idx>=0:
                    data[var] = mean_meas[:,idx] * self.params['coeff_density'] 
                else:
                    idx = self.vars['static_pressure']
                    data[var] =  ( mean_meas[:,idx] * self.params['weight_pressure_to_rho']
                                * self.params['coeff_press'] )
            if var in ['x_velocity','y_velocity','z_velocity']:
                data[var] =  mean_meas[:,idx] * self.params['coeff_vel'] 
        
        if self.probe is None:
            self.probe = dict()
            
        self.probe[self.format] = DataFrame(data=data)
        
        
