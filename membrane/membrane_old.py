#! /usr/bin/env python
import random as rand
import numpy as np
import vtk
import molar.pdb
#from membrane.martini_membrane import SURFACE_DENSITY

PATH     = "/Users/alinar/Dropbox/The Project/MicroBio/molecule_pool/"
FILES    =  ["lip1.pdb","lip2.pdb","lip3.pdb","cerebroside.pdb","cerns_martini_single_1.pdb",
             "fah24_martini.pdb","chol_martiniv2.0.pdb","cerns_martini_2218.pdb","cerns_martini_2618.pdb","cerns_martini_3018.pdb",
             "fah22_martini.pdb","fah26_martini.pdb","fah30_martini.pdb","dpgs_martini_new.pdb"]
CER_NS   =  ["cerns_martini_2218.pdb","cerns_martini_2418.pdb","cerns_martini_2618.pdb","cerns_martini_3018.pdb"]
CER_NP   =  ["cernp_martini_2218.pdb","cernp_martini_2418.pdb","cernp_martini_2618.pdb","cernp_martini_3018.pdb"]
GCER_NS  =  ["glycerns_martini_2218.pdb","glycerns_martini_2418.pdb","glycerns_martini_2618.pdb","glycerns_martini_3018.pdb"]
GCER_NP  =  ["glycernp_martini_2218.pdb","glycernp_martini_2418.pdb","glycernp_martini_2618.pdb","glycernp_martini_3018.pdb"]
EOS      =  "cereos_martini_hairpin.pdb"
GEOS     =  "glycereos_martini_hairpin.pdb" 
FAH      =  ["fah22_martini.pdb","fah24_martini.pdb","fah26_martini.pdb","fah30_martini.pdb"]
CHOL     =  "chol_martiniv2.0.pdb"

DISTANCE = 5.0  #7.0 # distance between two neighbouring ceremides. 
#SURFACE_DENSITY   = 2.0 / ((DISTANCE**2)*(3**0.5))
SURFACE_DENSITY    = 0.0638   # 0.053
UNIT_CELL_SIZE     = 320      # used in the function
WATER_MARGIN_SURF_VALUE = 1.19# 1.15  #This value should be changed for different unit cell sizes (1.15 for 200)

class Membrane(molar.pdb.Pdb):
    
    def __init__(self,number_of_cells_per_side_=1,unit_cell_size_ = 200):

        molar.pdb.Pdb.__init__(self)
        self.surface_density= SURFACE_DENSITY
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
        self.attempt       = 0
        self.collisioncounter = 0
        self.total_area    = 0
        self.total_pair    = 0
        for f in FILES:
            self.units.append(molar.pdb.Pdb())
            self.units[-1].ReadFile(PATH + f)
        for f in FILES:
            self.units_transed.append(molar.pdb.Pdb())
            self.units_transed[-1].ReadFile(PATH + f)
            self.units_transed[-1].RotateX(180)
        self.pointlocator   = vtk.vtkPointLocator()
        self.pointlocator.Initialize()
        self.pointlocator.SetDataSet(vtk.vtkPolyData())
        self.pointlocator.InitPointInsertion(vtk.vtkPoints(), [0,unit_cell_size_,0,unit_cell_size_,0,unit_cell_size_])
    
    def Update(self):
        self.side_length    = float(self.number_of_cells_per_side * self.unit_cell_size)
        self.sampling_size = self.side_length / (self.step_num - 1)
            
    def Func(self,_x_,_y_,_z_):
        if self.form == "G":
            return np.sin(_x_) * np.cos(_y_) + np.sin(_y_) * np.cos(_z_) + np.sin(_z_) * np.cos(_x_)
        elif self.form == "P":
            return np.cos(_x_) + np.cos(_y_) + np.cos(_z_)
        elif self.form == "D":
            return np.cos(_x_-_y_) * np.cos(_z_) + np.sin(_x_+_y_) * np.sin(_z_)
        
    def MarchingCubes(self):
        step    =    self.number_of_cells_per_side * 2 * np.pi / float(self.step_num)
        vol = vtk.vtkImageData()
        vol.SetDimensions(self.step_num , self.step_num , self.step_num )
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        for i in range(self.step_num):
            for j in range(self.step_num):
                for k in range(self.step_num):
                    f = self.Func(i*step,j*step,k*step)
                    vol.SetScalarComponentFromDouble(i,j,k,0,f)
        self.marching_cubes.SetInputData(vol)
        self.marching_cubes.ComputeNormalsOn()
        self.marching_cubes.SetValue(0, 0)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()

        self.poly_data = self.marching_cubes.GetOutput(0)
        self.normals   = self.poly_data.GetPointData().GetNormals()
        self.cell_n    = self.poly_data.GetNumberOfCells()
        print self.cell_n

    def CalcTransform(self,triangle,rotation_angle = False):
        center = np.zeros([3],dtype="f")
        trans  = vtk.vtkTransform()
        trans.PostMultiply() 
        for i in range(3):
            center = center + triangle[i]
        center = center * self.sampling_size / 3.0
        
        AreaVec = CalcAreaVec(triangle * self.sampling_size)
        normal  = AreaVec / np.linalg.norm(AreaVec) 
        if rotation_angle:
            trans.RotateY(rotation_angle)
        else:
            trans.RotateY(360 * rand.random())
        trans.Concatenate(molar.pdb.RotateToParallel(normal,[0,1,0]))
        trans.Translate(center)
        trans.Translate([1*rand.random(),1*rand.random(),1*rand.random()])
#         trans  = vtk.vtkTransform()
#         trans.PostMultiply() 
#         phi    = np.degrees(np.arccos(normal[1]))
#         if normal[0]<0.0001: #escape dividing by zero
#             theta = 90
#         else:
#             theta  = np.degrees(np.arctan(normal[2]/normal[0]) )
#         if (normal[0]<0):
#             theta = theta + 180
#         trans.RotateZ(-1*phi)
#         trans.RotateY(-1*theta)
#         # randomization
#         #trans.Translate([10*rand.random(),10*rand.random(),10*rand.random()])
#         trans.Translate(center)
        return trans

    def Make(self):
        self.Update()
        self.MarchingCubes()
        triangle   = np.zeros([3,3],dtype="f")
        count = 0
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)
            cell        = self.poly_data.GetCell(i)
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
            s = self.TriangleAnalysis(triangle)
            count = count + s
        self.MakeMolList()
        print "\ntotal area      : %f angstrom squered" % self.total_area
        #print "total ceramides   : %d " % self.mol_set['CER ']
        #print "area per cer      : %f nm^2 / CER" % (2*self.total_area/(100*self.mol_set['CER ']))
        #print "cer density       : %f CER / nm^2" % (100*self.mol_set['CER ']/(2*self.total_area))
        #print "ignored triangles : " , count
        
    def MakeMartini(self):
        """make with Martini coarse-grained molecule with cholesterol and fatty acid.
        """
        self.Update()
        self.MarchingCubes()
        triangle   = np.zeros([3,3],dtype="f")
        count = 0
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
            s = self.TriangleAnalysisMartini(triangle)
            count = count + s
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
        self.PrintInfo()
    
    def MakeMartini2(self):
        """make with Martini coarse-grained molecule without cholesterol and fatty acid.
        with varying lengths.
        """
        self.Update()
        self.MarchingCubes()
        triangle   = np.zeros([3,3],dtype="f")
        count = 0
        for i in range(self.cell_n): 
            molar.pdb.update_progress(float(i)/self.cell_n)        
            cell        = self.poly_data.GetCell(i)
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
            s = self.TriangleAnalysisMartini2(triangle)
            count = count + s
        self.PrintInfo()
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
    
    def MakeMartini3(self, glyco=0.0):

        """make with Martini coarse-grained molecule with cholesterol and fatty acid.
        NOT with varying lengths. with some glycoseramides. 
        """
        self.Update()
        self.MarchingCubes()
        triangle   = np.zeros([3,3],dtype="f")
        count = 0
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
            s = self.TriangleAnalysisMartini3(triangle,glyco)
            count = count + s
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
        self.PrintInfo()
        pass
    
    def MakeMartini4(self, glyco=0.0):

        """make with Martini coarse-grained molecule with cholesterol and fatty acid.
        AND with varying lengths. with some glycoseramides. 
        """

#         dist_   = np.array([0.6/4 , 2.47/4 , 0.41/4 , 0.52/4],dtype=float)
#         dist_   = dist_ / dist_.sum()
#         self.dist  = np.add.accumulate(dist_)
        ########### loading component files ###########
        self.cerns      = []
        self.cernp      = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah        = []
        self.chol       = molar.pdb.Pdb()
        self.chol.ReadFileGMX(PATH+CHOL)
        self.chol.BringToCenter()
        self.chol.Translate([0,10,0])
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(PATH+EOS)
        self.eos.BringToCenter()
        self.eos.Translate([0,10,0])
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(PATH+GEOS)
        self.geos.BringToCenter()
        self.geos.Translate([0,10,0])
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(PATH+file_name)
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(PATH+file_name)
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(PATH+file_name)
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(PATH+file_name)
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(PATH+file_name)
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05        # overall ceramide concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute distributions."""
        self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,  self.cerns[1]: h*(1-glyco)*d[1] * cer_conc, self.cerns[2]: h*(1-glyco)*d[2] * cer_conc, self.cerns[3]: h*(1-glyco)*d[3] * cer_conc, 
                       self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,  self.cernp[1]: h*(1-glyco)*d[1] * cer_conc, self.cernp[2]: h*(1-glyco)*d[2] * cer_conc, self.cernp[3]: h*(1-glyco)*d[3] * cer_conc, 
                       self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,    self.gcerns[2]:h*glyco*d[2] * cer_conc,     self.gcerns[3]:h*glyco*d[3] * cer_conc,
                       self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,    self.gcernp[2]:h*glyco*d[2] * cer_conc,     self.gcernp[3]:h*glyco*d[3] * cer_conc,
                       self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,            self.fah[2]:   d[2] * fah_conc,             self.fah[3]:   d[3] * fah_conc + extra_fah30,
                       self.chol:   chol_conc,
                       self.eos:  (1-glyco)*eos_conc,
                       self.geos:    glyco *eos_conc
                   }
        self.componants = []
        
        dist_           = np.array([])
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            dist_=np.append(dist_,value)
        dist_ = dist_ / dist_.sum()
        self.dist = np.add.accumulate(dist_)
        
        ################################################
        self.Update()
        self.MarchingCubes()
        triangle   = np.zeros([3,3],dtype="f")
        count = 0
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
             
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle[j])
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle * self.sampling_size))
            s = self.TriangleAnalysisMartini4(triangle,glyco)
            count = count + s
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
        self.PrintInfo()
        pass
    
    def TriangleAnalysis(self,t): # recursive
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area 
        rand_point_num =  np.random.poisson(point_num)
        if (rand_point_num == 1):
            self.AddMol(t)
        elif (rand_point_num > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysis(np.array([t[0],t[2],MidPoint(t[0],t[1])]))
                    self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[1])]))
                else:
                    self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
            else:
                if l2 > l3:
                    self.TriangleAnalysis(np.array([t[0],t[2],MidPoint(t[2],t[1])]))
                    self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[2],t[1])]))
                else:
                    self.TriangleAnalysis(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysis(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
        else: # rand_point_num == 0
            return 1
        return 0
        
    def TriangleAnalysisMartini(self,t): # recursive
        """ to be used with MakeMartini()
        """
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area 
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            trans=self.CalcTransform(t)
            random_number_1 = rand.random()
            random_number_2 = rand.random()
            if random_number_1 < 0.33: # add two ceramides
                for i in range(11):
                    success = self.CatTransformedCautious(self.units[8], trans, self.pointlocator, cutoff = 0.85) 
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                trans=self.CalcTransform(t)
                for i in range(11):
                    success = self.CatTransformedCautious(self.units_transed[8], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
            else: # add chol fatty acid pair  
                for i in range(11):
                    if random_number_2 < 0.5:
                        success = self.CatTransformedCautious(self.units[5], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units[6], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                for i in range(11):
                    if random_number_2 < 0.5:
                        success = self.CatTransformedCautious(self.units_transed[6], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units_transed[5], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
        
            self.total_pair = self.total_pair + 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini(np.array([t[0],t[2],MidPoint(t[0],t[1])]))
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[1])]))
                else:
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini(np.array([t[0],t[2],MidPoint(t[2],t[1])]))
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[2],t[1])]))
                else:
                    self.TriangleAnalysisMartini(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
        else: # rand_point_num == 0
            return 1
        return 0
    
    def TriangleAnalysisMartini2(self,t): # recursive
        """ to be used with MakeMartini2()
        """
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area 
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            trans=self.CalcTransform(t)
            random_number_1 = rand.random()
            random_number_2 = rand.random()
            if random_number_1 < 0.33: # add two ceramides
                random_number_3 = rand.random()
                for i in range(11):
                    if random_number_3 < 0.3/2:
                        success = self.CatTransformedCautious(self.units[7], trans, self.pointlocator, cutoff = 0.85)
                    elif 0.3/2 < random_number_3 < 1.44/2:
                        success = self.CatTransformedCautious(self.units[8], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units[9], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                trans=self.CalcTransform(t)
                random_number_3 = rand.random()
                for i in range(11):
                    if random_number_3 < 0.3/2:
                        success = self.CatTransformedCautious(self.units_transed[7], trans, self.pointlocator, cutoff = 0.85)
                    elif 0.3/2 < random_number_3 < 1.44/2:
                        success = self.CatTransformedCautious(self.units_transed[8], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units_transed[9], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
            else: # add chol fatty acid pair
                random_number_3 = rand.random()
                for i in range(11):
                    if random_number_2 < 0.5:
                        ##  insert fah  ##
                        if random_number_3 < 0.3/2:
                            success = self.CatTransformedCautious(self.units_transed[10], trans, self.pointlocator, cutoff = 0.85)
                        elif 0.3/2 < random_number_3 < 1.44/2:
                            success = self.CatTransformedCautious(self.units_transed[11], trans, self.pointlocator, cutoff = 0.85)
                        else:
                            success = self.CatTransformedCautious(self.units_transed[12], trans, self.pointlocator, cutoff = 0.85)
                        ## finished ##
                    else:
                            success = self.CatTransformedCautious(self.units[6], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                trans=self.CalcTransform(t)
                random_number_3 = rand.random()
                for i in range(11):
                    if random_number_2 < 0.5:
                        success = self.CatTransformedCautious(self.units_transed[6], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        ##  insert fah  ##
                        if random_number_3 < 0.3/2:
                            success = self.CatTransformedCautious(self.units_transed[10], trans, self.pointlocator, cutoff = 0.85)
                        elif 0.3/2 < random_number_3 < 1.44/2:
                            success = self.CatTransformedCautious(self.units_transed[11], trans, self.pointlocator, cutoff = 0.85)
                        else:
                            success = self.CatTransformedCautious(self.units_transed[12], trans, self.pointlocator, cutoff = 0.85)
                        ## finished ##
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
        
            self.total_pair = self.total_pair + 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini2(np.array([t[0],t[2],MidPoint(t[0],t[1])]))
                    self.TriangleAnalysisMartini2(np.array([t[1],t[2],MidPoint(t[0],t[1])]))
                else:
                    self.TriangleAnalysisMartini2(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini2(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini2(np.array([t[0],t[2],MidPoint(t[2],t[1])]))
                    self.TriangleAnalysisMartini2(np.array([t[0],t[1],MidPoint(t[2],t[1])]))
                else:
                    self.TriangleAnalysisMartini2(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini2(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
        else: # rand_point_num == 0
            return 1
        return 0
            
    def TriangleAnalysisMartini3(self,t,glyco=0.0): # recursive
        """ to be used with MakeMartini()
        """
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area 
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            trans=self.CalcTransform(t)
            random_number_1 = rand.random()
            random_number_2 = rand.random()
            if random_number_1 < 0.33: # add two ceramides
                for i in range(11):
                    if rand.random() > glyco: # add ceramide
                        success = self.CatTransformedCautious(self.units[8], trans, self.pointlocator, cutoff = 0.85) 
                    else:   #add glycoceramide
                        success = self.CatTransformedCautious(self.units[13], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                trans=self.CalcTransform(t)
                for i in range(11):
                    if rand.random() > glyco:
                        success = self.CatTransformedCautious(self.units_transed[8], trans, self.pointlocator, cutoff = 0.85) 
                    else:   #add glycoceramide
                        success = self.CatTransformedCautious(self.units[13], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
            else: # add chol fatty acid pair  
                for i in range(11):
                    if random_number_2 < 0.5:
                        success = self.CatTransformedCautious(self.units[5], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units[6], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
                ## insert the opposite molecule.
                for i in range(11):
                    if random_number_2 < 0.5:
                        success = self.CatTransformedCautious(self.units_transed[6], trans, self.pointlocator, cutoff = 0.85)
                    else:
                        success = self.CatTransformedCautious(self.units_transed[5], trans, self.pointlocator, cutoff = 0.85)
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
                else:
                    self.collisioncounter += 1
                self.attempt += 1
        
            self.total_pair = self.total_pair + 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[2],MidPoint(t[0],t[1])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[1])]),glyco)
                else:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[0],t[2])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[2])]),glyco)
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[2],MidPoint(t[2],t[1])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[2],t[1])]),glyco)
                else:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[0],t[2])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[2])]),glyco)
        else: # rand_point_num == 0
            return 1
        return 0
    
    def TriangleAnalysisMartini4(self,t,glyco=0.0): # recursive
        """ to be used with MakeMartini()
        """
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area 
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            trans=self.CalcTransform(t)
            for i in range(11):
                    j       = self.ProduceIndex()
                    new_mol = self.componants[j]
                    success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 0.85) 
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
            else:
                    self.collisioncounter += 1
            self.attempt += 1
            ## insert the opposite molecule.
            trans=vtk.vtkTransform()
            trans.PostMultiply()
            trans_aux=self.CalcTransform(t)
            trans.RotateX(180)
            trans.Concatenate(trans_aux)
            for i in range(11):
                    j       = self.ProduceIndex()
                    new_mol = self.componants[j]
                    success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 0.85) 
                    if success:
                        break
                    trans=self.CalcTransform(t,i*30)
            else:
                    self.collisioncounter += 1
            self.attempt += 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[2],MidPoint(t[0],t[1])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[1])]),glyco)
                else:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[0],t[2])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[2])]),glyco)
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[2],MidPoint(t[2],t[1])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[2],t[1])]),glyco)
                else:
                    self.TriangleAnalysisMartini3(np.array([t[0],t[1],MidPoint(t[0],t[2])]),glyco)
                    self.TriangleAnalysisMartini3(np.array([t[1],t[2],MidPoint(t[0],t[2])]),glyco)
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
    
    def SurfMarchingCubes(self):
        """ To make two surfaces on the headgoups.
        """
        step    =    self.number_of_cells_per_side * 2 * np.pi / float(self.step_num)
        self.s_marching_cubes = vtk.vtkMarchingCubes()
        self.s_marching_cubes.ComputeScalarsOff()
        vol = vtk.vtkImageData()
        vol.SetDimensions(self.step_num , self.step_num , self.step_num ) 
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        for i in range(self.step_num):
            for j in range(self.step_num):
                for k in range(self.step_num):
                    f = self.Func(i*step,j*step,k*step)
                    vol.SetScalarComponentFromDouble(i,j,k,0,f)
        self.s_marching_cubes.SetInputData(vol)
        self.s_marching_cubes.SetValue(0, WATER_MARGIN_SURF_VALUE)
        self.s_marching_cubes.SetValue(1, -1*WATER_MARGIN_SURF_VALUE)
        self.s_marching_cubes.Update()
        self.s_marching_cubes.UpdateInformation()
        self.surf_actor = vtk.vtkActor()
        scaleTrans=vtk.vtkTransform()
        scaleTrans.Scale(self.sampling_size,self.sampling_size,self.sampling_size)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.s_marching_cubes.GetOutputPort())
        self.surf_actor.SetMapper(mapper)
        self.surf_actor.GetProperty().SetColor(0.0,0.95,0.0)
        self.surf_actor.SetUserTransform(scaleTrans)
        pass
    
    def AddWater(self):
        water_box  = molar.pdb.Pdb()
        water_box.ReadFileGMX(PATH + "water_antifreeze.pdb")
        box_length = water_box.Bounds()[0,1] - water_box.Bounds()[0,0]
        trans = vtk.vtkTransform()
        trans.PostMultiply()
        for x in np.arange(0,self.unit_cell_size,box_length):
            for y in np.arange(0,self.unit_cell_size,box_length):
                for z in np.arange(0,self.unit_cell_size,box_length):
                    trans.Identity()
                    trans.Translate([x,y,z])
                    self.CatTransformed(water_box, trans)
                pass
            pass
        ## remove redundant
        global UNIT_CELL_SIZE
        UNIT_CELL_SIZE = self.unit_cell_size #change the global value to match.
        self.RemoveAtomByFunc(DelFunc)
        print "Pruning after removing extra water..."
        ## Prune:  [pdb.Prune() takes a long time!]
        mol_num = len(self.molecules)
        i   = 0 
        new_mol_list=[]
        for mol in self.molecules:
            molar.pdb.update_progress(float(i)/mol_num) 
            if len(mol.chains[0].residues[0].atoms) != 0:
                new_mol_list.append(mol)
                pass
            i+=1
        self.molecules=new_mol_list
            
    def Show(self):
        ## override
        self.SurfMarchingCubes()
        molar.pdb.Pdb.Show(self,"sphere",ext_actors=[self.surf_actor])
        pass
    
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

    def ProduceIndex(self):
        r        = rand.random()
        for i in range(len(self.dist)):
            if r < self.dist[i]:
                break
        return i
        
##########################################################

def Func_aux(pos):
        return np.sin(pos[0]) * np.cos(pos[1]) + np.sin(pos[1]) * np.cos(pos[2]) + np.sin(pos[2]) * np.cos(pos[0])

def DelFunc(atom):
        # works only for the box size = 200
        global UNIT_CELL_SIZE
        if ((atom.name.replace(" ","") == "W" or atom.name.replace(" ","") == "WF")
            and 
            (abs(Func_aux( (atom.pos/UNIT_CELL_SIZE) * 2 * np.pi )) <= WATER_MARGIN_SURF_VALUE
            or atom.pos[0]<0 or atom.pos[0]>UNIT_CELL_SIZE
            or atom.pos[1]<0 or atom.pos[1]>UNIT_CELL_SIZE
            or atom.pos[2]<0 or atom.pos[2]>UNIT_CELL_SIZE) ):
            return True
        else:
            return False

def MidPoint(v1,v2):
    return 0.5 * (v1 + v2)

def CalcAreaVec(in_array):
    """returning normal vector of a triangle with sides of v1 and v2
       that has the length equal to the area of the triangle. 
    """
    coor = np.array(in_array)
    v1 = coor[1] - coor[0]
    v2 = coor[2] - coor[0]
    c  = 0.5*np.cross(v1,v2)
    return c #np.linalg.norm(c)

def Distance(v1,v2):
    return np.linalg.norm(v1-v2)


################################################################
#     def CalcTransform_old(self,cell_num):
#         coor   = np.zeros([4,3],dtype="f")
#         center = np.zeros([3],dtype="f")
#         normal = np.zeros([3],dtype="f")
#         cell   = self.poly_data.GetCell(cell_num)
#         N      = cell.GetNumberOfPoints()
#         for i in range(N):
#             point_id    =    cell.GetPointId(i)
#             self.poly_data.GetPoints().GetPoint(point_id,coor[i])
#             center = center + coor[i]
#             normal = normal + np.array(self.normals.GetTuple(point_id))
#             
#         center = center * self.side_length / (N * (self.step_num - 1))
#         normal = normal / np.linalg.norm(normal) 
#         trans  = vtk.vtkTransform()
#         trans.PostMultiply() 
#         phi    = np.degrees(np.arccos(normal[1]))
#         theta  = np.degrees(np.arctan(normal[2]/normal[0]) )
#         if (normal[0]<0):
#             theta = theta + 180
#         trans.RotateZ(-1*phi)
#         trans.RotateY(-1*theta)
#         # randomization
#         trans.Translate([10*rand.random(),10*rand.random(),10*rand.random()])
#         trans.Translate(center)
#         area = CalcArea(coor * self.side_length / (N * (self.step_num - 1)))
#         return [trans,area]

#     def Make_old(self):
#         self.MarchingCubes()
#         count = 0
#         for i in range(self.cell_n):
#             pair         = self.RandomPair()
#             [trans,area] = self.CalcTransform(i)
#             if area < 10:
#                 count = count + 1
#                 continue
#             self.CatTransformed(pair[0], trans)
#         print count