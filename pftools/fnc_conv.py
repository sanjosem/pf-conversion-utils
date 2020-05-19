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
    def __init__(self,pfFile,verbose=True):
        super().__init__(pfFile,verbose)
        self.format = 'volume'
        self.cell_coords = None
        self.node_coords = None
        self.cell_conn = None
        self.vertex_to_node = None
        self.volume_cell = None
        self.data = None
        self.domain = None
        
        if self.verbose:
            print('PF file format is: {0:s}'.format(self.format))

    def read_domains(self):
        """Method to read the domains contained in the powerflow file

        """
        from scipy.io import netcdf
        from numpy import unique
        
        if self.params is None:
            self.read_conversion_parameters()
            
        
        self.domain = dict()
        if self.rotation:
        
            f = netcdf.Dataset(self.pfFile,'r',mmap=False)
            
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
        
        from scipy.io import netcdf
        import numpy as np 
        from copy import deepcopy
        
        if self.domain is None:
            self.read_domains()

        f = netcdf.Dataset(self.pfFile, 'r')

        element_coords = f.variables['coords'][()].astype('float')
        voxel_scales = f.variables['voxel_scales'][()]
            
        f.close()
        
        # Compute cell connectivity
        if self.verbose:
            print("Compute center coordinates")
            
        nb_vertex = 8

        self.params['ncells'] = element_coords.shape[0]
        self.params['ndims'] = element_coords.shape[1]
        assert self.params['ndims']==3, "Wrong coordinate dimensions"
        
        if not self.rotation:
            self.domain['stator'] = np.arange(self.params['ncells'])
        
        # Offset coordinates
        for idim in range(self.params['ndims']):
            element_coords[:,idim]+=self.params['offset_coords'][idim]
                
        self.cell_coords = dict()
        self.node_coords = dict()
        self.cell_conn = dict()
        self.vertex_to_node = dict()
        self.volume_cell = dict()
        
        for dom in self.domain.keys():
            lst = self.domain[dom]
            
            # scale coordinates
            self.cell_coords[dom] = element_coords[lst,:] * self.params['coeff_dx']
            
            if self.verbose:
                print('Computing vertices coordinates for {0:s}'.format(dom))
            dx = 2**(1.0 + voxel_scales[lst])
            self.volume_cell[dom] = (self.params['coeff_dx']*dx)**3
            nelm = lst.size
            if self.verbose:
                print("  -> {0:d} elements".format(nelm))

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
            
            coords_view = np.ascontiguousarray(vertices_coords.round(4)).view(
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
                            
        return vertices_coords
            
    def read_frame_data(self,frame):
        
        from scipy.io import netcdf
        import numpy as np
        import pandas as pd
        
        if self.vars is None:
            self.define_measurement_variables()
            
        if self.cell_conn is None:
            self.read_volume_mesh()
            
        # before py2f wrapper
        if not type(frame) == int:
            raise RuntimeError('frame input must be an integer')
        
        
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
                stats = pd.DataFrame(data=self.data[dom].mean(axis=-1),columns=self.vars.keys())
                print('  -> Stats (nodes)')
                print(stats)
        
        
        f.close()
        

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
        import os.path
        
        super().save_parameters(casename,dirout)

        outFile = os.path.join(dirout,'param_pf_{0:s}.hdf5'.format(casename))
        print("Adding surface mesh convertion arrays:\n  ->  {0:s}".format(outFile))
        
        fparams = h5py.File(outFile,'a')
        if self.face_conn is not None:
            gcon = fparams.create_group("connectivity")
            gcon.create_dataset('glo_vertex_number', data=self.face_conn['vert_per_face'])
            gcon.create_dataset('glo_face_vertex_list', data=self.face_conn['face_vertex_list'])
            gcon.create_dataset('glo_node_surface', data=self.face_conn['node_weight'])
            gcon.create_dataset('glo_face_weight', data=self.face_conn['face_weight'])
            fparams.flush()

        if self.node_coords is not None:
            geo = fparams.create_group("coordinates")
            geo.create_dataset('node_coords', data=self.node_coords)
            fparams.flush()

        if self.mesh is not None:
            for surface_name in self.mesh.keys():
                res = self.mesh[surface_name]
                subgrp = gcon.create_group('surface_name')
                subgrp.create_dataset('lst_faces', data=res['glo_faces']) #, dtype='i8')
                subgrp.create_dataset('glo_node_list', data=res['glo_nodes'])
            
        fparams.close()


    def save_vtk(self,casename,dirout):
        """Method to export volume mesh for paraview

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        from pftools.module_vtk_utils import cells_to_vtkUnstruct,save_MultiBlock,data_to_vtkBlock
        import os.path

        outFile = os.path.join(dirout,'surface_mesh_{0:s}.vtm'.format(casename))
        print("Exporting in VTK format:\n  ->  {0:s}".format(outFile))

        list_unStruct_Blocks = dict()
        for dom in self.domain.keys():
            loc_nodes = self.node_coords[dom]
            cell_elm = self.cell_conn[dom]

            list_unStruct_Blocks[dom] = cells_to_vtkUnstruct(loc_nodes,cell_elm)
            
            if self.data is not None:
                list_unStruct_Blocks[dom] = data_to_vtkBlock(
                    list_unStruct_Blocks[dom],self.data[dom][0,:,:],self.vars.keys())

        save_MultiBlock(list_unStruct_Blocks,outFile)
