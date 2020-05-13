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

def save_MBPolyData(struct_PolyData, outFile):
    
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
