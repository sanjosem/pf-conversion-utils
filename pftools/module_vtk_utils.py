import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

def faces_to_vtkPolyData(coords,tri_elm,quad_elm):
    import numpy as np

    mesh = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(coords.astype(np.float64), deep=1))

    mesh.SetPoints(points)

    cells = vtk.vtkCellArray()
    
    ntri = tri_elm.shape[0]
    nqua = quad_elm.shape[0]
    ncells = ntri+nqua

    tris =  np.c_[np.tile(3, ntri),tri_elm].flatten().astype(np.int64)
    quads =  np.c_[np.tile(4, nqua),quad_elm].flatten().astype(np.int64)

    faces = np.concatenate((tris,quads))
    faceIds = numpy_to_vtkIdTypeArray(faces, deep=1)

    cells.SetCells(ncells,faceIds)
    
    mesh.SetPolys(cells)
    
    return mesh
    

def cells_to_vtkUnstruct(coords,cell_elm):
    import numpy as np

    # import vtk
    # from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    # coords = fncConv.node_coords['stator']
    # cell_elm = fncConv.cell_conn['stator']

    mesh = vtk.vtkUnstructuredGrid()

    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(coords.astype(np.float64), deep=1))

    mesh.SetPoints(points)

    cells = vtk.vtkCellArray()
    
    ncells = cell_elm.shape[0]

    hexs =  np.c_[np.tile(8, ncells),cell_elm].flatten().astype(np.int64)
    
    cell_type = vtk.vtkHexahedron().GetCellType()
    cell_types = np.tile(cell_type, (ncells, 1))

    cellIds = numpy_to_vtkIdTypeArray(hexs, deep=1)

    cells.SetCells(ncells,cellIds)
    
    mesh.SetCells(cell_type,cells)
    
    return mesh    

def data_to_vtkBlock(mesh,data,vars):

    for ivar,var in enumerate(vars):
        parray = numpy_to_vtk(data[ivar,:])
        parray.SetName(var)
        mesh.GetPointData().AddArray(parray)

    return mesh

def save_polydata(polydata, file_name, binary=True, color_array_name=None):
    
    writer = vtk.vtkXMLPolyDataWriter()

    writer.SetFileName(file_name)
    writer = writer.SetInputData(polydata)
    if color_array_name is not None:
        writer.SetArrayName(color_array_name)

    if binary:
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAscii()
            
    writer.Update()
    print('Writing file {0:s}'.format(file_name))
    writer.Write()

def save_MultiBlock(struct_PolyData, outFile):
    
    surface_names = struct_PolyData.keys()
    nb_blk = len(surface_names)

    mb = vtk.vtkMultiBlockDataSet()
    mb.SetNumberOfBlocks(nb_blk)

    for i,surfn in enumerate(surface_names):
        mb.SetBlock(i, struct_PolyData[surfn])
        mb.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), surfn ) 

    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(outFile)
    writer.SetInputData(mb)
    writer.Write()

def cut_by_plane(vtk_obj,origin=(0.,0.,0.),norm=(1.,0.,0.)):
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin[0],origin[1],origin[2])
    plane.SetNormal(norm[0],norm[1],norm[2])
    
    planeCut = vtk.vtkCutter()
    planeCut.SetInputData(vtk_obj)
    planeCut.SetCutFunction(plane)
    planeCut.Update()
    
    return planeCut.GetOutput()
    
def find_closest_point(vtk_obj,point=(0.,0.,0.)):
    
    pointTree = vtk.vtkKdTreePointLocator()
    pointTree.SetDataSet(vtk_obj)
    pointTree.BuildLocator()
    
    radius = None
    if radius is None:
        vtkId = pointTree.FindClosestPoint(point)
        pcoord = vtk_obj.GetPoint(vtkId)
        
    return vtkId,pcoord
