#! /usr/bin/env python
import random as rand
import numpy as np
import vtk
import molar.pdb

PATH     = "/Users/alinar/Dropbox/The Project/MicroBio/molecule_pool/"
FILES    =  ["lip1.pdb","lip2.pdb","lip3.pdb","cerebroside.pdb","cerns_martini_single_1.pdb"]
DISTANCE =  7.0   # distance between two neighbouring ceremides. 
SURFACE_DENSITY = 2.0 / ((DISTANCE**2)*(3**0.5))
SURFACE_DENSITY = 0.1 # 0.2

class Membrane(molar.pdb.Pdb):
    
    def __init__(self,layer_dist=1.115,number_of_cells_per_side_=1,surface_density = 0.1,unit_cell_size_ = 300):
        molar.pdb.Pdb.__init__(self)
        self.layer_dist=layer_dist
        self.step_num       = int(32)
        self.marching_cubes = vtk.vtkMarchingCubes()
        self.poly_data      = vtk.vtkPolyData()
        self.data           = vtk.vtkDataObject()
        self.unit_cell_size = unit_cell_size_
        self.form           = "G"
        self.number_of_cells_per_side = number_of_cells_per_side_
        self.cerebroside    = False
        self.Update()
        self.units         = []
        self.units_transed = []
        self.points        = vtk.vtkPoints() #for visualization
        self.surf_dens     = surface_density
     
        self.total_area    = 0
        self.total_pair    = 0
        for f in FILES:
            self.units.append(molar.pdb.Pdb())
            self.units[-1].ReadFile(PATH + f)
        for f in FILES:
            self.units_transed.append(molar.pdb.Pdb())
            self.units_transed[-1].ReadFile(PATH + f)
            self.units_transed[-1].RotateX(180)
        self.units[-1].BringToNegativeY()
        #self.units[-1].Show()
    
    def Update(self):
        self.side_length    = float(self.number_of_cells_per_side * self.unit_cell_size)
        self.sampling_size = self.side_length / (self.step_num - 1)
            
    def Func(self,pos):
        _x_ = pos[0]
        _y_ = pos[1]
        _z_ = pos[2]
        if self.form == "G":
            return np.sin(_x_) * np.cos(_y_) + np.sin(_y_) * np.cos(_z_) + np.sin(_z_) * np.cos(_x_)
        elif self.form == "P":
            return np.cos(_x_) + np.cos(_y_) + np.cos(_z_)
        elif self.form == "D":
            return np.cos(_x_-_y_) * np.cos(_z_) + np.sin(_x_+_y_) * np.sin(_z_)
        
    def Gradient(self,pos):
        step      =    self.number_of_cells_per_side * 2 * np.pi / float(self.step_num)
        _x_ = pos[0] * step
        _y_ = pos[1] * step
        _z_ = pos[2] * step
        out = np.zeros([3],dtype="f")
        if self.form == "G":
            out[0] =   np.cos(_x_) * np.cos(_y_) - np.sin(_z_) * np.sin(_x_)
            out[1] =   np.cos(_y_) * np.cos(_z_) - np.sin(_x_) * np.sin(_y_)
            out[2] =   np.cos(_x_) * np.cos(_z_) - np.sin(_z_) * np.sin(_y_)
        elif self.form == "P":
            out[0] =    - np.sin(_x_)
            out[1] =    - np.sin(_y_) 
            out[2] =    - np.sin(_z_) 
        elif self.form == "D":
            out[0] =    -np.sin(_x_-_y_) * np.cos(_z_) + np.sin(_x_+_y_) * np.sin(_z_)
            out[1] =     np.sin(_x_-_y_) * np.cos(_z_) + np.cos(_x_+_y_) * np.sin(_z_)
            out[2] =    -np.cos(_x_-_y_) * np.sin(_z_) + np.sin(_x_+_y_) * np.cos(_z_) 
        return out
        
    def MarchingCubes(self,C=0):
        step      =    self.number_of_cells_per_side * 2 * np.pi / float(self.step_num)
        vol = vtk.vtkImageData()
        vol.SetDimensions(self.step_num , self.step_num , self.step_num )
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        for i in range(self.step_num):
            for j in range(self.step_num):
                for k in range(self.step_num):
                    f = self.Func([i*step,j*step,k*step])
                    vol.SetScalarComponentFromDouble(i,j,k,0,f)
        self.marching_cubes.SetInputData(vol)
        self.marching_cubes.ComputeNormalsOff()
        self.marching_cubes.ComputeGradientsOff()
        self.marching_cubes.SetValue(0, C)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()
        self.poly_data = self.marching_cubes.GetOutput(0)
        self.normals   = self.poly_data.GetPointData().GetNormals()
        self.cell_n    = self.poly_data.GetNumberOfCells()
        #print type(self.normals)
        #print self.normals.GetNumberOfTuples(),",", self.normals.GetNumberOfComponents()
        #print "cell_n:",self.cell_n
        #print "GetComputeNormals():",self.marching_cubes.GetComputeNormals()

    def CalcTransform(self,triangle):
        center = np.zeros([3],dtype="f")
        trans  = vtk.vtkTransform()
        trans.PostMultiply() 
        for i in range(3):
            center = center + triangle[i]
        center = center * self.sampling_size / 3.0
        
        AreaVec = CalcAreaVec(triangle * self.sampling_size)
        normal  = AreaVec / np.linalg.norm(AreaVec) 
        
        trans.RotateY(360 * rand.random())
        trans.Concatenate(molar.pdb.RotateToParallel(normal,[0,1,0]))
        trans.Translate(center)
        trans.Translate([1*rand.random(),1*rand.random(),1*rand.random()])
        return trans
    
    def CalcTransform2(self,triangle,sign):
        center      = TriangleMidPoint(triangle)
        direction   = sign * self.Gradient(center)
        trans  = vtk.vtkTransform()
        trans.PostMultiply() 
        
        #AreaVec = CalcAreaVec(triangle * self.sampling_size)
        #normal  = AreaVec / np.linalg.norm(AreaVec) 
        
        trans.RotateY(360 * rand.random())
        trans.Concatenate(molar.pdb.RotateToParallel(direction,[0,1,0]))
        trans.Translate(center * self.sampling_size)
        #trans.Translate([1*rand.random(),1*rand.random(),1*rand.random()])
        return trans

    def MakeMartini(self):
        """make with Martini coarse-grained molecule
        """
        self.Update()
        for sign in [-1,1]:
            self.MarchingCubes(sign * 0.5*self.layer_dist)
            triangle   = np.zeros([3,3],dtype="f")
            count = 0
            for i in range(self.cell_n):                
                cell  = self.poly_data.GetCell(i)
                for j in range(3):
                    point_id    =    cell.GetPointId(j)
                    self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
                self.points.InsertNextPoint(TriangleMidPoint(triangle))
                self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
                s = self.TriangleAnalysisMartini(triangle,sign)
                count = count + s
            
        print "ignored triangles   : " , count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "total pairs         : ", self.total_pair 
        if (self.total_area != 0 and self.total_pair != 0):
            print "average density     : ", self.total_pair / self.total_area
            print "average pair's area : ", self.total_area / self.total_pair
        
    def TriangleAnalysisMartini(self,t,sign): # recursive
        """ to be used with MakeMartini()
        """
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surf_dens * area 
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            trans=self.CalcTransform2(t,sign)
            self.CatTransformed(self.units[4], trans)
            #self.CatTransformed(self.units_transed[4], trans)
            self.total_pair = self.total_pair + 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini(np.array([t[0],t[2],MidPoint(t[0],t[1])]),sign)
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[1])]),sign)
                else:
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[0],t[2])]),sign)
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[2])]),sign)
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini(np.array([t[0],t[2],MidPoint(t[2],t[1])]),sign)
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[2],t[1])]),sign)
                else:
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[0],t[2])]),sign)
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[2])]),sign)
        else: # rand_point_num == 0
            return 1
        return 0
    
    def AddMol(self,triangle):
        pair = self.RandomPair()
        trans=self.CalcTransform(triangle)
        self.CatTransformed(pair[0], trans)
        self.CatTransformed(pair[1], trans)
        self.total_pair = self.total_pair + 1
        for side in triangle:#for visualization
            self.points.InsertNextPoint(side)
            
    def RandomPair(self):
        r = rand.random()
        a = 1.0 / 3.0
        if r <  a:
            if self.cerebroside == True:
                return [self.units[3],self.units_transed[3]]
            else:
                return [self.units[0],self.units_transed[0]]
        elif a < r < (2*a):
            return [self.units[1],self.units_transed[2]]
        else:
            return [self.units[2],self.units_transed[1]]
    
    def VisVertices(self):
        actor     = vtk.vtkActor()
        glyph3d   = vtk.vtkGlyph3D()
        point_s   = vtk.vtkPointSource()
        mapper    = vtk.vtkPolyDataMapper()
        window    = vtk.vtkRenderWindow()
        renderer  = vtk.vtkRenderer()
        interactor= vtk.vtkRenderWindowInteractor()
        poly_data = vtk.vtkPolyData()                
        point_s.SetNumberOfPoints(1)
        interactor.SetRenderWindow(window)
        poly_data.SetPoints(self.points)
        glyph3d.SetSourceConnection(point_s.GetOutputPort())
        glyph3d.SetInputData(poly_data)
        mapper.SetInputConnection(glyph3d.GetOutputPort())
        actor.SetMapper(mapper)
        window.AddRenderer(renderer)
        renderer.AddActor(actor)
        renderer.SetBackground(0.1,0.2,0.3)
        renderer.ResetCamera()
        window.Render()
        interactor.Start()

#     def TriangleAnalysis(self,t): # recursive
#         area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
#         l1   = np.linalg.norm(t[0] - t[1])
#         l2   = np.linalg.norm(t[1] - t[2])
#         l3   = np.linalg.norm(t[0] - t[2])
#         point_num      = self.surf_dens * area 
#         rand_point_num =  np.random.poisson(point_num)
#         if (rand_point_num == 1):
#             self.AddMol(t)
#         elif (rand_point_num > 1):
#             if l1 > l2:
#                 if l1 > l3:
#                     self.TriangleAnalysis(np.array([t[0],t[2],MidPoint(t[0],t[1])]))
#                     self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[1])]))
#                 else:
#                     self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
#                     self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
#             else:
#                 if l2 > l3:
#                     self.TriangleAnalysis(np.array([t[0],t[2],MidPoint(t[2],t[1])]))
#                     self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[2],t[1])]))
#                 else:
#                     self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
#                     self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
#         else: # rand_point_num == 0
#             return 1
#         return 0

#     def Make(self):
#         self.Update()
#         self.MarchingCubes()
#         triangle   = np.zeros([3,3],dtype="f")
#         count = 0
#         for i in range(self.cell_n):                
#             cell        = self.poly_data.GetCell(i)
#             for j in range(3):
#                 point_id    =    cell.GetPointId(j)
#                 self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
#             s = self.TriangleAnalysis(triangle)
#             count = count + s
#         #print "ignored triangles : " , count
        
        
##########################################################

def MidPoint(v1,v2):
    return 0.5 * (v1 + v2)

def TriangleMidPoint(triangle):
    mid   = np.zeros(3,dtype="f")
    for point in triangle:
        mid = mid + point
    return mid / 3.0

def CalcAreaVec(in_array):
    """returning normal vector of a triangle with sides of v1 and v2
       that has the length of the area of the triangle. 
    """
    coor = np.array(in_array)
    v1 = coor[1] - coor[0]
    v2 = coor[2] - coor[0]
    c  = 0.5*np.cross(v1,v2)
    return c #np.linalg.norm(c)

def Distance(v1,v2):
    return np.linalg.norm(v1-v2)

