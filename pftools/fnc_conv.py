from pftools.module_conv_utils import PFConversion

class fncConversion(PFConversion):
    """Class for volume mesh conversion
    inherits from PFConversion

    Parameters
    ----------
    pfFile : string
        powerFLOW surface file
    verbose : bool
        Activate or desactivate informative prints (default: True)

    Methods:
    ----------
    read_domains
        Identify rotating and stationnary domain
    read_volume_mesh
        Generate the volume mesh and connectivity
    read_frame_data
        Convert data in SI and at mesh nodes
    save_vtk
        Create multiblock VTK format (UnstructuredGrid)


    """
    def __init__(self,pfFile,verbose=True,use_fapi=None):
        super().__init__(pfFile,verbose)
        self.format = 'volume'
        self.cell_coords = None
        self.node_coords = None
        self.cell_conn = None
        self.vertex_to_node = None
        self.volume_cell = None
        self.data = None
        self.domain = None
        self.vtk_object = None


        if use_fapi is None:
            try:
                import pftools.fextend.fnc_reader as Ffnc
                self.fapi = True
            except:
                self.fapi = False
        else:
            self.fapi=use_fapi


        if self.verbose:
            print('PF file format is: {0:s}'.format(self.format))

    def read_domains(self):
        """Method to read the domains contained in the powerflow file

        """
        import netCDF4 as netcdf
        import numpy as np

        if self.params is None:
            self.read_conversion_parameters()


        self.domain = dict()
        if self.rotation:

            f = netcdf.Dataset(self.pfFile,'r')

            ref_frame_indices = f.variables['ref_frame_indices'] # -1 stationnary 0 LRF
            self.domain['rotor'] = np.where(ref_frame_indices[:] == 0)[0]
            self.domain['stator'] = np.where(ref_frame_indices[:] == -1)[0]

            f.close()


    def read_volume_mesh(self,store=True):
        """Method to read and scale node coordinates. The array is stored in the class

        Parameters
        ----------
        store : bool
            if True (default), coordinates are stored in the class

        """

        import netCDF4 as netcdf
        import numpy as np
        from copy import deepcopy
        if self.fapi:
            import pftools.fextend.fnc_reader as Ffnc

        if self.domain is None:
            self.read_domains()

        f = netcdf.Dataset(self.pfFile, 'r')

        self.params['ncells'] = f.dimensions['npoints'].size
        self.params['ndims'] = f.dimensions['ndims'].size

        if not self.fapi:
            element_coords = f.variables['coords'][()].astype('float')
            voxel_scales = f.variables['voxel_scales'][()]

        f.close()

        assert self.params['ndims']==3, "Wrong coordinate dimensions"

        # Compute cell connectivity
        if self.verbose:
            print("Compute center coordinates")

        nb_vertex = 8

        if not self.rotation:
            self.domain['stator'] = np.arange(self.params['ncells'])

        # Offset coordinates
        if not self.fapi:
            for idim in range(self.params['ndims']):
                element_coords[:,idim]+=self.params['offset_coords'][idim]

        self.cell_coords = dict()
        self.node_coords = dict()
        self.cell_conn = dict()
        self.vertex_to_node = dict()
        self.volume_cell = dict()

        for dom in self.domain.keys():
            lst = self.domain[dom]
            nelm = lst.size

            if self.verbose:
                print("  -> {0:d} elements".format(nelm))

            if self.verbose:
                print('Computing vertices coordinates for {0:s}'.format(dom))

            if nelm>0:
                if self.fapi:
                    cell_volumes,cell_coords,vertices_coords = Ffnc.read_fnc_mesh(self.pfFile,
                                                self.params['coeff_dx'],self.params['offset_coords'],
                                                self.params['ncells'],lst)

                    # scale coordinates
                    self.cell_coords[dom] = cell_coords
                    self.volume_cell[dom] = cell_volumes

                else:

                    # scale coordinates
                    self.cell_coords[dom] = element_coords[lst,:] * self.params['coeff_dx']

                    dx = 2**(1.0 + voxel_scales[lst])
                    self.volume_cell[dom] = (self.params['coeff_dx']*dx)**3

                    vertices_coords = np.zeros((nelm * 8, 3))
                    basis = np.eye(3)
                    vx = dx[:,np.newaxis]*basis[:,0][np.newaxis,:]
                    vy = dx[:,np.newaxis]*basis[:,1][np.newaxis,:]
                    vz = dx[:,np.newaxis]*basis[:,2][np.newaxis,:]

                    vertices_coords[0::8, :] = element_coords[lst,:]
                    vertices_coords[1::8, :] = element_coords[lst,:] + vx
                    vertices_coords[2::8, :] = element_coords[lst,:] + vx + vy
                    vertices_coords[3::8, :] = element_coords[lst,:]      + vy
                    vertices_coords[4::8, :] = element_coords[lst,:]           + vz
                    vertices_coords[5::8, :] = element_coords[lst,:] + vx      + vz
                    vertices_coords[6::8, :] = element_coords[lst,:] + vx + vy + vz
                    vertices_coords[7::8, :] = element_coords[lst,:]      + vy + vz

                if self.verbose:
                    print("Compute connectivity")

                # this part is done for both py and fapi
                coords_view = np.ascontiguousarray(vertices_coords.round(5)).view(
                                np.dtype((np.void, vertices_coords.dtype.itemsize * vertices_coords.shape[1])))
                _, idx, inv = np.unique(coords_view, return_index=True, return_inverse=True)

                coords_unique = self.params['coeff_dx']*vertices_coords[idx]

                connectivity = np.arange(nelm * 8).reshape((nelm, 8))
                new_connectivity = inv[connectivity]
                nnodes = coords_unique.shape[0]
                if self.verbose:
                    print("  -> {0:d} nodes".format(nnodes))

                self.node_coords[dom] = coords_unique
                self.cell_conn[dom] = new_connectivity
                self.vertex_to_node[dom] = deepcopy(idx)

                if self.verbose:
                    print('Bounding box in meters for domain {0:s}:'.format(dom))
                    for idim,coor in enumerate(['x','y','z']):
                        print('  {2:s}: [{0:.3e},{1:.3e}]'.format(
                                self.cell_coords[dom][:,idim].min(),self.cell_coords[dom][:,idim].max(),coor))

            else:
                print('No element for domain {0:s}'.format(dom))

        return vertices_coords

    def read_frame_data(self,frame):

        import netCDF4 as netcdf
        import numpy as np
        import pandas as pd
        if self.fapi:
            import pftools.fextend.fnc_reader as Ffnc

        if self.vars is None:
            self.define_measurement_variables()

        if self.cell_conn is None:
            self.read_volume_mesh()

        # before py2f wrapper
        if not type(frame) == int:
            raise RuntimeError('frame input must be an integer')

        if not self.fapi:
            slicing_instant = slice(frame,frame+1)

            nvars = len(self.vars.keys())

            f = netcdf.Dataset(self.pfFile, 'r')

            for dom in self.domain.keys():

                print('Converting data for domain \'{0:s}\''.format(dom))

                lst = self.domain[dom]

                # Get data (selected faces, selected variables)
                tmp = f.variables['measurements'][slicing_instant,:,lst]

                # Conversion in SI
                data_cell = np.zeros((tmp.shape[0],nvars,tmp.shape[-1]))
                for ivar,var in enumerate(self.vars.keys()):
                    idx = self.vars[var]

                    if var == 'static_pressure':
                        if idx>=0 :
                            data_cell[:,ivar,:] = ( ( tmp[:,idx,:] + self.params['offset_pressure'] )
                                                  * self.params['coeff_press'] )
                        else:
                            idx = self.vars['density']
                            data_cell[:,ivar,:] = ( tmp[:,idx,:] * self.params['weight_rho_to_pressure']
                                             + self.params['offset_pressure'] ) * self.params['coeff_press']

                    if var == 'density':
                        if idx>=0:
                            data_cell[:,ivar,:] = tmp[:,idx,:] * self.params['coeff_density']
                        else:
                            idx = self.vars['static_pressure']
                            data_cell[:,ivar,:]  =  ( tmp[:,idx,:] * self.params['weight_pressure_to_rho']
                                               * self.params['coeff_density'] )
                    if var in ['x_velocity','y_velocity','z_velocity']:
                        data_cell[:,ivar,:] =  tmp[:,idx,:] * self.params['coeff_vel']

                if self.verbose:
                    stats = pd.DataFrame(data=data_cell.mean(axis=-1),columns=self.vars.keys())
                    print('  -> Stats (cell)')
                    print(stats)

                ninst,nvars,ncells = data_cell.shape

                # Cell to node
                new_connectivity = self.cell_conn[dom]
                nnodes = self.node_coords[dom].shape[0]

                data_node = np.zeros((ninst,nvars,nnodes))

                scale_nvert = 1./8.
                cell_weight = self.volume_cell[dom]*scale_nvert

                volume_node = np.zeros((nnodes,))

                for node_idx in range(8):
                    if self.verbose:
                        print('  -> vertex {0:d}/8'.format(node_idx+1))

                    data_node[:,:,new_connectivity[:,node_idx]] += data_cell*cell_weight
                    volume_node[new_connectivity[:, node_idx]] += cell_weight

                # To avoid divide by 0 error
                eps = self.volume_cell[dom].min()/100.
                volume_node[volume_node<eps] = eps
                # Scale
                data_node = data_node/volume_node[np.newaxis,np.newaxis,:]

                # Storage
                if self.data is None:
                    self.data = dict()

                self.data[dom] = data_node

                if self.verbose:
                    print(self.data[dom].shape)
                    print(self.data[dom].mean(axis=-1).shape)
                    stats = pd.DataFrame(data=self.data[dom].mean(axis=-1),columns=self.vars.keys())
                    print('  -> Stats (nodes)')
                    print(stats)


            f.close()
        else:
            # Store in fortran module
            self.params_to_fapi()

            # Prepare var extraction
            retrieve_index = []
            scale_type = []

            for ivar,var in enumerate(self.vars.keys()):
                idx = self.vars[var]
                if var == 'static_pressure':
                    if idx>=0 :
                        retrieve_index.append(idx)
                        scale_type.append(int(Ffnc.pf_params.type_pressure))
                    else:
                        retrieve_index.append(self.vars['density'])
                        scale_type.append(int(Ffnc.pf_params.type_pressure_from_density))
                if var == 'density':
                    if idx>=0 :
                        retrieve_index.append(idx)
                        scale_type.append(int(Ffnc.pf_params.type_density))
                    else:
                        retrieve_index.append(self.vars['density'])
                        scale_type.append(int(Ffnc.pf_params.type_density_from_pressure))
                if var in ['x_velocity','y_velocity','z_velocity']:
                    retrieve_index.append(idx)
                    scale_type.append(int(Ffnc.pf_params.type_velocity))

            for dom in self.domain.keys():

                if self.domain[dom].size>0:

                    nnodes = self.node_coords[dom].shape[0]
                    ncells = self.params['ncells']
                    
                    print('Converting data for domain \'{0:s}\''.format(dom))

                    data_node = Ffnc.read_fnc_frame(self.pfFile,frame,
                                    nnodes,ncells,self.domain[dom],self.volume_cell[dom],
                                    self.cell_conn[dom],
                                    retrieve_index,scale_type)


                    # Storage
                    if self.data is None:
                        self.data = dict()

                    nvars = len(retrieve_index)
                    self.data[dom] = data_node.reshape(1,nvars,nnodes)

                    if self.verbose:
                        stats = pd.DataFrame(data=self.data[dom].mean(axis=-1),columns=self.vars.keys())
                        print('  -> Stats (nodes)')
                        print(stats)

                else: 
                    print(f'Nothing to do for {dom}')


    def params_to_fapi(self):
        if self.fapi:
            import pftools.fextend.fnc_reader as Ffnc

        if self.params is None:
            self.read_conversion_parameters()

        if self.fapi:
            Ffnc.pf_params.coeff_dx     = self.params['coeff_dx']
            Ffnc.pf_params.timestep     = self.params['dt']
            Ffnc.pf_params.coeff_vel    = self.params['coeff_vel']
            Ffnc.pf_params.coeff_rho    = self.params['coeff_density']
            Ffnc.pf_params.coeff_press  = self.params['coeff_press']
            Ffnc.pf_params.offset_press = self.params['offset_pressure']
            Ffnc.pf_params.weight_p2r   = self.params['weight_pressure_to_rho']
            Ffnc.pf_params.weight_r2p   = self.params['weight_rho_to_pressure']

    def save_parameters(self,casename,dirout):
        """Method to export convertion parameters in a separated hdf5 file.

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        import h5py

        outFile = super().save_parameters(casename,dirout)
        print("Adding volume mesh conversion arrays into:\n  ->  {0:s}".format(outFile))

        fparams = h5py.File(outFile,'a')
        if self.data is not None and self.vertex_to_node is not None:
            gdom = fparams.create_group("domains")
            for dom in self.domain.keys():
                gcur = gdom.create_group(dom)
                gcur.create_dataset('glo_cell_indices', data=self.domain[dom])
                gcur.create_dataset('cell_coords', data=self.cell_coords[dom])
                gcur.create_dataset('volume_cell', data=self.volume_cell[dom])
                gcur.create_dataset('node_coords', data=self.node_coords[dom])
                gcur.create_dataset('connectivity', data=self.cell_conn[dom])
                gcur.create_dataset('vertex_to_node', data=self.vertex_to_node[dom])

                fparams.flush()

        fparams.close()
        return outFile


    def save_instant(self,outFile,frame):
        """Method to export convertion parameters in a separated hdf5 file.

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        import h5py

        print("Adding variable data into:\n  ->  {0:s}".format(outFile))

        fparams = h5py.File(outFile,'a')
        if self.data is not None:
            if not 'data' in list(fparams.keys()):
                gdata = fparams.create_group("data")
            else:
                gdata = fparams['data']

            instant = f'{frame:04d}'

            if instant in list(gdata.keys()):
                gframe = gdata[instant]
            else:
                gframe = gdata.create_group(instant)
            
            if self.vars is None:
                self.define_measurement_variables()
                
            if self.time is None:
                self.extract_time_info()

            print(self.time)

            gframe.create_dataset('time',data=self.time[frame])

            for dom in self.domain.keys():

                if self.domain[dom].size>0:

                    gcur = gframe.create_group(dom)

                    for ivar,var in enumerate(self.vars.keys()):
                        gcur.create_dataset(var, data=self.data[dom][0,ivar,:])
                
                        fparams.flush()

        fparams.close()


    def load_parameters(self,h5file):

        import h5py

        pfFile = self.pfFile
        verbosity = self.verbose
        fapi = self.fapi

        self.__init__(pfFile,verbose=verbosity,use_fapi=fapi)

        super().load_parameters(h5file)

        fparams = h5py.File(h5file,'r')

        print("Loading the fnc connectivity and geometry file:\n  ->  {0:s}".format(h5file))

        if fparams.get('domains',getclass=True) is None:
            RuntimeError('No domains group file, {0:s} is is not valid for fnc'.format(self.format))

        if not fparams.get('domains',getclass=True) is None:
            if fparams.get('domains',getclass=True) == h5py.Group:
                self.domain = dict()
                self.cell_coords = dict()
                self.volume_cell = dict()
                self.node_coords = dict()
                self.cell_conn = dict()
                self.vertex_to_node = dict()
                for dom in ['stator','rotor']:
                    if not fparams.get('domains/{0:s}'.format(dom),getclass=True) is None:
                        if fparams.get('domains/{0:s}'.format(dom),getclass=True) == h5py.Group:

                            grp = fparams['domains/{0:s}'.format(dom)]
                            self.domain[dom] = grp['glo_cell_indices'][()]
                            self.cell_coords[dom] = grp['cell_coords'][()]
                            self.volume_cell[dom] = grp['volume_cell'][()]
                            self.node_coords[dom] = grp['node_coords'][()]
                            self.cell_conn[dom] = grp['connectivity'][()]
                            self.vertex_to_node[dom] = grp['vertex_to_node'][()]


        fparams.close()

    def create_vtk(self):
        """Method to create VTK volume mesh

        Parameters
        ----------
        """
        from pftools.module_vtk_utils import cells_to_vtkUnstruct,data_to_vtkBlock

        if self.domain is None:
            self.read_domains()

        if self.cell_conn is None:
            self.read_volume_mesh()

        list_unStruct_Blocks = dict()

        for dom in self.domain.keys():

            if self.domain[dom].size>0:
                loc_nodes = self.node_coords[dom]
                cell_elm = self.cell_conn[dom]

                list_unStruct_Blocks[dom] = cells_to_vtkUnstruct(loc_nodes,cell_elm)

                if self.data is not None:
                    list_unStruct_Blocks[dom] = data_to_vtkBlock(
                        list_unStruct_Blocks[dom],self.data[dom][0,:,:],self.vars.keys())

        self.vtk_object = list_unStruct_Blocks

    def save_vtk(self,casename,dirout):
        """Method to export volume mesh for paraview

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        from pftools.module_vtk_utils import save_MultiBlock
        import os.path

        if self.vtk_object is None:
            self.create_vtk()

        outFile = os.path.join(dirout,'surface_mesh_{0:s}.vtm'.format(casename))
        print("Exporting in VTK format:\n  ->  {0:s}".format(outFile))

        save_MultiBlock(self.vtk_object,outFile)


    def extract_probe(self,probe_name,probe_coords):

        if self.fapi:
            import pftools.fextend.fnc_reader as Ffnc

        import numpy as np
        from pandas import DataFrame
        import netCDF4 as netcdf

        if self.time is None:
            self.extract_time_info()

        if self.vars is None:
            self.define_measurement_variables()

        dloc = 1.0e8
        for dom in self.domain.keys():
            if self.fapi:
                icell_closest,dist = Ffnc.find_closest_cell(probe_coords,self.cell_coords[dom])
            else:
                dom_dist = np.linalg.norm(np.asarray(probe_coords)[np.newaxis,:]
                                            - self.cell_coords[dom])
                icell_closest = dom_dist.argmin()
                dist = dom_dist[icell_closest]

            if dist<dloc:
                dloc = dist
                pid = icell_closest
                domn = dom

        print('Probe {0:s} found in {1:s} at distance {2:e} m'.format(
                probe_name,domn,dloc))

        glo_id = self.domain[domn][pid]

        data = dict()
        data['time'] = self.time['time_center']

        if self.fapi:

            # Store in fortran module
            self.params_to_fapi()

            # Prepare var extraction
            retrieve_index = []
            scale_type = []

            for ivar,var in enumerate(self.vars.keys()):
                idx = self.vars[var]
                if var == 'static_pressure':
                    if idx>=0 :
                        retrieve_index.append(idx)
                        scale_type.append(int(Ffnc.pf_params.type_pressure))
                    else:
                        retrieve_index.append(self.vars['density'])
                        scale_type.append(int(Ffnc.pf_params.type_pressure_from_density))
                if var == 'density':
                    if idx>=0 :
                        retrieve_index.append(idx)
                        scale_type.append(int(Ffnc.pf_params.type_density))
                    else:
                        retrieve_index.append(self.vars['density'])
                        scale_type.append(int(Ffnc.pf_params.type_density_from_pressure))
                if var in ['x_velocity','y_velocity','z_velocity']:
                    retrieve_index.append(idx)
                    scale_type.append(int(Ffnc.pf_params.type_velocity))

            temporal_data = Ffnc.get_cell_data(glo_id,self.time['nsets'],self.pfFile,retrieve_index,scale_type)

            for ivar,var in enumerate(self.vars.keys()):
                data[var] = temporal_data[:,ivar]

        else:
            f = netcdf.Dataset(self.pfFile, 'r')
            tmp = f.variables['measurements'][:,:,glo_id]
            f.close()

            nvars = len(self.vars.keys())

            # Conversion in SI
            data_cell = np.zeros((tmp.shape[0],nvars))
            for ivar,var in enumerate(self.vars.keys()):
                idx = self.vars[var]

                if var == 'static_pressure':
                    if idx>=0 :
                        data_cell[:,ivar] = ( ( tmp[:,idx] + self.params['offset_pressure'] )
                                              * self.params['coeff_press'] )
                    else:
                        idx = self.vars['density']
                        data_cell[:,ivar] = ( tmp[:,idx] * self.params['weight_rho_to_pressure']
                                         + self.params['offset_pressure'] ) * self.params['coeff_press']

                if var == 'density':
                    if idx>=0:
                        data_cell[:,ivar] = tmp[:,idx] * self.params['coeff_density']
                    else:
                        idx = self.vars['static_pressure']
                        data_cell[:,ivar]  =  ( tmp[:,idx] * self.params['weight_pressure_to_rho']
                                           * self.params['coeff_density'] )
                if var in ['x_velocity','y_velocity','z_velocity']:
                    data_cell[:,ivar] =  tmp[:,idx] * self.params['coeff_vel']

            for ivar,var in enumerate(self.vars.keys()):
                data[var] = data_cell[:,ivar]

        if self.probe is None:
            self.probe = dict()

        self.probe[probe_name] = DataFrame(data=data)
