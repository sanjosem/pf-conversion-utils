'''Class to store parameters and functions for the convertion of PowerFLOW files

    Author: Marlene Sanjose
    '''
    


class PFConvertion:
    """Initiate the convertion class associated to a given PowerFLOW file.

    Parameters
    ----------
    pfFile : string
        File name (and path) of the PowerFLOW file to convert
    verbose : bool
        Activate or desactivate informative prints (default: True)
            
    Methods:
    ----------
    read_surface_names()
        Read surfaces contains in file
    read_conversion_parameters()
        Read and compute convertion parameters
    save_parameters(outFile)
        Save convertion parameters
    """
    def __init__(self,pfFile,verbose=True):

        from os.path import exists as fexists
        if not fexists(pfFile):
            raise RuntimeError('File {0:s} could not be found'.format(pfFile))
            
        ext = pfFile.split('.')[-1]
        if ext == 'snc':
            self.format = 'surface'
            self.sncFile = pfFile
        else:
            raise NotImplementedError('Only snc file are supported')
            
        self.params = None
        self.rotation = False
        self.surface_list = None
        self.verbose = verbose
        self.node_coords = None
        self.face_conn = None
        self.mesh = None
        

    def read_face_names(self):
        """Function to read the surface patchs (called faces) contained in the powerflow file

        Returns
        -------
        list_av_face_names: list
            List of the available faces
        list_av_faces_ids: list
            List of the corresponding ids

        """
        from scipy.io import netcdf
        from numpy import unique
        
        f = netcdf.netcdf_file(self.sncFile,'r',mmap=False)
        
        x = f.variables['face_names'][()]
        nb_fnames = f.dimensions['nfaces']
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
        
    def read_conversion_parameters(self):
        """Function to read conversion factor, parameters are stored in params dictionnary in the class

        """
        
        from scipy.io import netcdf
        
        f = netcdf.netcdf_file(self.sncFile, 'r', mmap=False)
        
        lx_offsets = f.variables['lx_offsets'][()]
        lx_scales = f.variables['lx_scales'][()]
        param_values = f.variables['param_values'][()]

        print("Computing scaling factors")
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
        from numpy import where
        
        f = netcdf.netcdf_file(self.sncFile, 'r', mmap=False)
            
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
            
            n_rot0 = f.variables['lrf_initial_n_revolutions'][0] * np.pi * 2.0
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
        

    def read_coordinates(self,store=True):
        """Function to read and scale node coordinates. The array is stored in the class

        Parameters
        ----------
        store : bool
            if True (default), coordinates are stored in the class

        """
        
        from scipy.io import netcdf

        f = netcdf.netcdf_file(self.sncFile, 'r', mmap=False)

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
        """Function to read and prepare connectivity.
        Arrays are stored in the class

        """
        
        from scipy.io import netcdf
        import numpy as np
        
        if self.node_coords is None:
            self.read_coordinates()
        
        self.face_conn = dict()

        f = netcdf.netcdf_file(self.sncFile, 'r', mmap=False)
        
        first_vertex = f.variables['first_vertex_refs'][()]
        vertex_list = f.variables['vertex_refs']
        self.params['nfaces'] = first_vertex.size
        self.params['nvertices'] = vertex_list[:].size
        
        vert_per_face = np.zeros((self.params['nfaces'],),dtype='int')
        vert_per_face[:-1] = first_vertex[1:] - first_vertex[:-1]
        vert_per_face[-1] = self.params['nvertices'] - first_vertex[-1]
        max_vertex_per_face = vert_per_face.max()
        min_vertex_per_face = vert_per_face.min()

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
        
        face_weight = self.params['coeff_dx']**2 / vert_per_face.astype('float')

        node_weight = np.zeros((self.params['nnodes'],))
        for nvert in range(min_vertex_per_face,max_vertex_per_face):
            lst_face = np.where(vert_per_face == nvert)[0]
            for iv in range(nvert):
                node_weight[face_vertex_list[lst_face,iv]] += face_weight[lst_face]
                
        self.face_conn['face_weight'] = face_weight
        self.face_conn['node_weight'] = node_weight
        self.face_conn['face_area'] = f.variables['surfel_area'][()] * self.params['coeff_dx']**2
        self.face_conn['face_norm'] = f.variables['surfel_normal'][()].astype('float')
        
        f.close()
        
    def triangulate_surface(self,surface_name,face_list,min_threshold=1.0e-15):
        """Function to generate mesh of the given surface.
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
        
        from scipy.io import netcdf
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
        
        f = netcdf.netcdf_file(self.sncFile, 'r', mmap=False)
        face_id=f.variables['face'][()]
        f.close()
        
        # Get facet list
        lst_face = np.empty((0,),dtype='int')
        for face_name in face_list:
            id_face = self.surface_list['ids'][self.surface_list['names'].index(face_name)]
            lst_face_part = np.where((face_id==id_face))[0]
            lst_face = np.concatenate((lst_face,lst_face_part))

        surf_nfaces=lst_face.size            
        
        vert_per_face = self.face_conn['vert_per_face']
        face_vertex_list = self.face_conn['face_vertex_list']
        
        # subgrp = grp.create_group(surface_name)
        # subgrp.create_dataset('lst_faces', data=lst_face, dtype='i8')
        
        if self.verbose:
            print('Triangulation of {0:s}'.format(surface_name,surf_nfaces))

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
        """Function to export surface mesh for paraview

        Parameters
        ----------
        casename : string
            Label of the configuration
        dirout : string
            Directory for file export

        """
        from module_vtk_utils import faces_to_vtkPolyData,save_MBPolyData
        import os.path
        
        outFile = os.path.join(dirout,'surface_mesh_{0:s}.vtm'.format(casename))
        print("Exporting in VTK format:\n  ->  {0:s}".format(outFile))

        list_polyData_Blocks = dict()
        for surface_name in self.mesh.keys():
            res = self.mesh[surface_name]
            loc_nodes = self.node_coords[res['glo_nodes']]
            tri = res['loc_tri']
            qua = res['loc_qua']
        
            list_polyData_Blocks[surface_name] = faces_to_vtkPolyData(loc_nodes,tri,qua)
        
        save_MBPolyData(list_polyData_Blocks,outFile)
