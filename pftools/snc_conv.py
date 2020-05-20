from pftools.module_conv_utils import PFConversion

class sncConversion(PFConversion):
    """Class for surface mesh conversion
    inherits from PFConversion

    Parameters
    ----------
    pfFile : string
        powerFLOW surface file
    verbose : bool
        Activate or desactivate informative prints (default: True)

    Methods:
    ----------
    read_surface_names
        Read surfaces contained in file
    read_coordinates
        Read and store node coordinates
    read_connectivity
        Read, prepare and store connectivity
    triangulate_surface
        Triangulate surface elements that are neither triangles or quads
    save_vtk
        Create multiblock VTK format (PolyData)

    """
    def __init__(self,pfFile,verbose=True):
        super().__init__(pfFile,verbose)
        self.format = 'surface'
        self.surface_list = None
        self.node_coords = None
        self.face_conn = None
        self.mesh = None
        self.data = None
        self.vtk_object = None
        
        if self.verbose:
            print('PF file format is: {0:s}'.format(self.format))

    def read_surface_names(self):
        """Method to read the surface patchs (called faces) contained in the powerflow file

        Returns
        -------
        list_av_face_names: list
            List of the available faces
        list_av_faces_ids: list
            List of the corresponding ids

        """
        import netCDF4 as netcdf
        from numpy import unique
        
        f = netcdf.Dataset(self.pfFile,'r')
        
        x = f.variables['face_names'][()]
        nb_fnames = f.dimensions['nfaces'].size
        face_names = x.tostring().decode('utf-8').split('\x00')[:nb_fnames]
        face_id=f.variables['face'][()]
        
        f.close()

        list_av_faces_ids, face_counts = unique(face_id,return_counts=True)
        list_av_face_names = []
        
        if self.verbose:
            print('Available faces for extraction:')
        for idx,nf in zip(list_av_faces_ids,face_counts):
            list_av_face_names.append(face_names[idx])
            if self.verbose:
                print('  ({2:d}) {0:s} nb facets {1:d}'.format(face_names[idx],nf,idx))
                
        self.surface_list = dict()
        self.surface_list['names'] = list_av_face_names
        self.surface_list['ids'] = list_av_faces_ids
        self.surface_list['face_counts'] = face_counts

    def read_coordinates(self,store=True):
        """Method to read and scale node coordinates. The array is stored in the class

        Parameters
        ----------
        store : bool
            if True (default), coordinates are stored in the class

        """
        
        import netCDF4 as netcdf

        f = netcdf.Dataset(self.pfFile, 'r')

        coords=f.variables['vertex_coords'][()]
        
        f.close()
        
        self.params['nnodes'] = coords.shape[0]
        self.params['ndims'] = coords.shape[1]
        assert self.params['ndims']==3, "Wrong coordinate dimensions"
        
        # Offset coordinates
        for idim in range(self.params['ndims']):
            coords[:,idim]+=self.params['offset_coords'][idim]
        
        # scale coordinates    
        coords *= self.params['coeff_dx']
        
        if self.verbose:
            print('Bounding box in meters:')
            for idim,coor in enumerate(['x','y','z']):
                print('  {2:s}: [{0:.3e},{1:.3e}]'.format(coords[:,idim].min(),coords[:,idim].max(),coor))
                
        self.node_coords = coords.astype('float')

    def read_connectivity(self):
        """Method to read and prepare connectivity.
        Arrays are stored in the class

        """
        
        import netCDF4 as netcdf
        import numpy as np
        
        if self.node_coords is None:
            self.read_coordinates()
        
        self.face_conn = dict()

        f = netcdf.Dataset(self.pfFile, 'r')
        
        first_vertex = f.variables['first_vertex_refs'][()]
        vertex_list = f.variables['vertex_refs']
        self.params['nfaces'] = first_vertex.size
        self.params['nvertices'] = vertex_list[:].size
        
        vert_per_face = np.zeros((self.params['nfaces'],),dtype='int')
        vert_per_face[:-1] = first_vertex[1:] - first_vertex[:-1]
        vert_per_face[-1] = self.params['nvertices'] - first_vertex[-1]
        max_vertex_per_face = vert_per_face.max()
        min_vertex_per_face = vert_per_face.min()
        self.params['min_vertex_per_face'] = min_vertex_per_face
        self.params['max_vertex_per_face'] = min_vertex_per_face

        if self.verbose:
            print('Reading connectivity')
            print('  Total number of faces: {0:d}'.format(self.params['nfaces']))
            print('  Min/Max vertices per faces: {0:d} {1:d}'.format(min_vertex_per_face,max_vertex_per_face))

        # Compute global connectivity
        face_vertex_list = np.zeros((self.params['nfaces'],max_vertex_per_face),dtype=np.dtype('int64'))
        for iface in range(self.params['nfaces']-1):
            face_vertex_list[iface,0:vert_per_face[iface]]=vertex_list[first_vertex[iface]:first_vertex[iface+1]]
        face_vertex_list[-1,0:vert_per_face[-1]]=vertex_list[first_vertex[-1]:]
        
        self.face_conn['vert_per_face'] = vert_per_face
        self.face_conn['face_vertex_list'] = face_vertex_list
        
        surfel_area = f.variables['surfel_area'][()] * self.params['coeff_dx']**2
        face_weight = surfel_area / vert_per_face.astype('float')

        node_weight = np.zeros((self.params['nnodes'],))
        for nvert in range(min_vertex_per_face,max_vertex_per_face+1):
            lst_face = np.where(vert_per_face == nvert)[0]
            for iv in range(nvert):
                node_weight[face_vertex_list[lst_face,iv]] += face_weight[lst_face]
                
        self.face_conn['face_weight'] = face_weight
        self.face_conn['node_weight'] = node_weight
        self.face_conn['face_area'] = surfel_area
        self.face_conn['face_norm'] = f.variables['surfel_normal'][()].astype('float')
        
        f.close()
        
    def triangulate_surface(self,surface_name,face_list,min_threshold=1.0e-15):
        """Method to generate mesh of the given surface.
        Any element with more than 4 nodes are triangulated
        Results are stored in the class

        Parameters
        ----------
        surface_name : string
            name of the surface assembly
        face_list : list of strings
            list of the surface names to be gather in the assembly
            the surface names must exists in the snc file to be converted
        min_threshold : float
            Threshold for triangulation. Elements below this thimit will be droped

        """
        
        import netCDF4 as netcdf
        import numpy as np
        from scipy.spatial.distance import pdist, squareform
        from copy import deepcopy
        from scipy.spatial.qhull import Delaunay,QhullError 

        
        if self.surface_list is None:
            self.read_face_names()
        
        if self.node_coords is None:
            self.read_coordinates()
        if self.face_conn is None:
            self.read_connectivity()
        
        f = netcdf.Dataset(self.pfFile, 'r')
        face_id=f.variables['face'][()]
        f.close()
        
        print('Triangulation of \'{0:s}\' surface'.format(surface_name))
        # Get facet list
        lst_face = np.empty((0,),dtype='int')
        for face_name in face_list:
            id_face = self.surface_list['ids'][self.surface_list['names'].index(face_name)]
            print('  + \'{0:s}\' face'.format(face_name))
            lst_face_part = np.where((face_id==id_face))[0]
            lst_face = np.concatenate((lst_face,lst_face_part))

        surf_nfaces=lst_face.size            
        
        
        vert_per_face = self.face_conn['vert_per_face']
        face_vertex_list = self.face_conn['face_vertex_list']
        

        if surf_nfaces>0:

            tri_elm=np.zeros((4*surf_nfaces,3),dtype=np.dtype('int64')) # Larger array to accept triangulation

            tri_faces  = np.where(vert_per_face[lst_face]==3)[0] # Existing triangles
            quad_faces = np.where(vert_per_face[lst_face]==4)[0] # Existing quads
            bad_faces  = np.where(vert_per_face[lst_face]>4)[0] # Existing bad faces

            ntri  = len(tri_faces)
            nquad = len(quad_faces)
            nbad  = len(bad_faces)

            tri_elm[0:ntri,:] = face_vertex_list[lst_face[tri_faces],:3].astype('int64')
            quad_elm = face_vertex_list[lst_face[quad_faces],:4].astype('int64')

            if self.verbose:
                print("  Initial facets:")
                print("   -> Tri: {0:d}".format(ntri))
                print("   -> Quad: {0:d}".format(nquad))
                print("   -> Poly: {0:d}".format(nbad))
                print("   -> Total: {0:d}".format(ntri+nquad+nbad))

            count_null_surface = 0
            smax_ignore = 0.
            rel_err_area = 0.

            for glo_face_num in lst_face[bad_faces]:
                
                nvert = vert_per_face[glo_face_num]
                glo_vertex_list=face_vertex_list[glo_face_num,:nvert]
                point_list=self.node_coords[glo_vertex_list,:]
                npoints = point_list.shape[0]

                n = self.face_conn['face_norm'][glo_face_num,:] # normalized vector
                area = self.face_conn['face_area'][glo_face_num] # normalized vector
                
                # Find the most distant vertices of the face
                idx_1, idx_2 = np.unravel_index(np.argmax(squareform(pdist(point_list))), (npoints,npoints))
                
                # Create a (u,v) 2D basis describing the face
                u=point_list[idx_1,:]-point_list[idx_2,:]
                u=u/np.linalg.norm(u)
                
                v=np.cross(n,u) # u,v,n is an orthonormal base
                
                UV_pts=np.zeros((npoints,2))
                UV_pts[:,0]=np.dot(point_list,u)
                UV_pts[:,1]=np.dot(point_list,v)
                
                check_area = 0.
                try:
                    tri=Delaunay(UV_pts)
                    for simplex in tri.simplices:
                        simplex_list=point_list[simplex,:]
                        a=simplex_list[1,:]-simplex_list[0,:]
                        b=simplex_list[2,:]-simplex_list[0,:]
                        S=0.5*np.linalg.norm(np.cross(a,b))
                        if S > min_threshold:
                            tri_elm[ntri,:]=glo_vertex_list[simplex]
                            ntri=ntri+1
                            check_area += S
                            
                    rel_err_area = max(rel_err_area,abs(check_area - area)/area)
                    
                except QhullError:
                    pass
                    count_null_surface+=1
                    print('No connectivity found for face #{0:d}'.format(glo_face_num))

                    
            if rel_err_area > 5.0e-3:
                print('  Max relative difference in facet area after triangulation: {0:.2f} %'.format(100.*rel_err_area))
            if count_null_surface>0:
                print("  Number of facet that have been ignored: {0:d}".format(count_null_surface))

            # remove the additional lines
            tri_elm = np.delete(tri_elm,slice(ntri,None),0)
            
            if self.verbose:
                print("  After triangulation facets:")
                print("   -> Tri: {0:d}".format(ntri))
                print("   -> Quad: {0:d}".format(nquad))
                print("   -> Total: {0:d}".format(ntri+nquad))
            
            results = dict()
            glo_lst_node=np.unique(np.hstack((tri_elm.ravel(),quad_elm.ravel())))

            tri_connectivity = glo_lst_node.searchsorted(tri_elm)
            quad_connectivity = glo_lst_node.searchsorted(quad_elm)
            
            results['glo_faces'] = lst_face
            results['glo_nodes'] = glo_lst_node
            results['glo_tri'] = tri_elm
            results['glo_qua'] = quad_elm
            results['loc_tri'] = tri_connectivity
            results['loc_qua'] = quad_connectivity
            results['loc_nodes'] = glo_lst_node.size
            
            if self.mesh is None:
                self.mesh = dict()
            self.mesh[surface_name] = results
            
    def read_frame_data(self,surface_name,frame):
        
        import netCDF4 as netcdf
        import numpy as np
        import pandas as pd
        
        if self.vars is None:
            self.define_measurement_variables()
            
        if self.face_conn is None:
            self.read_connectivity()
        
        if self.mesh is None:
            raise RuntimeError('Mesh must be generated first using method triangulate_surface')
            
        if not surface_name in self.mesh.keys():
            raise RuntimeError('Mesh must be generated first using method triangulate_surface for {0:s}'.format(surface_name))
            
        # before py2f wrapper
        if not type(frame) == int:
            raise RuntimeError('frame input must be an integer')
        
        print('Converting data of \'{0:s}\' surface'.format(surface_name))
        
        slicing_instant = slice(frame,frame+1)
        
        nvars = len(self.vars.keys())
        
        lst_face = self.mesh[surface_name]['glo_faces']
        vert_per_face = self.face_conn['vert_per_face'][lst_face]
        
        # Get data (selected faces, selected variables)
        f = netcdf.Dataset(self.pfFile, 'r')
        tmp = f.variables['measurements'][slicing_instant,:,lst_face]
        f.close()

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
        nnodes = self.params['nnodes']
        
        # Scatter / Gather operation
        data_node = np.zeros((ninst,nvars,nnodes))
        surface_node = np.zeros((nnodes,))
        
        for nf,iface in enumerate(lst_face):
            nvert = self.face_conn['vert_per_face'][iface]
            face_weight = self.face_conn['face_weight'][iface]
            glo_nodes = self.face_conn['face_vertex_list'][iface,:nvert]
            data_node[:,:,glo_nodes] += data_cell[:,:,nf][:,:,np.newaxis]*face_weight
            surface_node[glo_nodes] += face_weight

        # To avoid divide by 0 error
        eps = self.face_conn['face_area'].min()/100.
        surface_node[surface_node<eps] = eps
        # Scale
        data_node = data_node/surface_node[np.newaxis,np.newaxis,:]
        
        # Storage
        if self.data is None:
            self.data = dict()
            
        lst_node = self.mesh[surface_name]['glo_nodes']
        self.data[surface_name] = data_node[:,:,lst_node]

        if self.verbose:
            stats = pd.DataFrame(data=self.data[surface_name].mean(axis=-1),columns=self.vars.keys())
            print('  -> Stats (nodes)')
            print(stats)
        

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


    def create_vtk(self):
        """Method to create vtk object and store it in the class instance

        """
        from pftools.module_vtk_utils import faces_to_vtkPolyData,data_to_vtkBlock
        import os.path
        
        list_polyData_Blocks = dict()
        for surface_name in self.mesh.keys():
            res = self.mesh[surface_name]
            loc_nodes = self.node_coords[res['glo_nodes']]
            tri = res['loc_tri']
            qua = res['loc_qua']

            list_polyData_Blocks[surface_name] = faces_to_vtkPolyData(loc_nodes,tri,qua)
            
            if surface_name in self.data.keys():
                list_polyData_Blocks[surface_name] = data_to_vtkBlock(
                    list_polyData_Blocks[surface_name],self.data[surface_name][0,:,:],self.vars.keys())

        self.vtk_object = list_polyData_Blocks
        
    def save_vtk(self,casename,dirout):
        """Method to export surface mesh for paraview

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
