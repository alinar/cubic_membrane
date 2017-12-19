#! /usr/bin/env python
import random as rand
import numpy as np
import vtk
import molar,molar.pdb
import os
import sys
import math

molar.Minium_version_required("2.0.0")

import molecule_pool
PATH     =  molecule_pool.__path__[0] 

FILES    =  ["lip1.pdb","lip2.pdb","lip3.pdb","cerebroside.pdb","cerns_martini_single_1.pdb",
             "fah24_martini.pdb","chol_martiniv2.0.pdb","cerns_martini_2218.pdb","cerns_martini_2618.pdb","cerns_martini_3018.pdb",
             "fah22_martini.pdb","fah26_martini.pdb","fah30_martini.pdb","dpgs_martini_new.pdb"]
CER_NS   =  ["cerns_martini_2218.pdb","cerns_martini_2418.pdb","cerns_martini_2618.pdb","cerns_martini_3018.pdb"]
CER_NP   =  ["cernp_martini_2218.pdb","cernp_martini_2418.pdb","cernp_martini_2618.pdb","cernp_martini_3018.pdb"]
GCER_NS  =  ["glycerns_martini_2218.pdb","glycerns_martini_2418.pdb","glycerns_martini_2618.pdb","glycerns_martini_3018.pdb"]
GCER_NP  =  ["glycernp_martini_2218.pdb","glycernp_martini_2418.pdb","glycernp_martini_2618.pdb","glycernp_martini_3018.pdb"]
EOS        =  "cereos_martini_hairpin.pdb"
GEOS     =  "glycereos_martini_hairpin.pdb" 
FAH        =  ["fah22_martini.pdb","fah24_martini.pdb","fah26_martini.pdb","fah30_martini.pdb"]
CHOL     =  "chol_martiniv2.0.pdb"

sqrt3          =  np.sqrt(3.0)
sqrt3o2      =  0.5 * np.sqrt(3.0)

#DISTANCE          = 5.0  #7.0 # distance between two neighbouring ceremides. 
#SURFACE_DENSITY   = 2.0 / ((DISTANCE**2)*(3**0.5))
SURFACE_DENSITY    = 0.0317 # 0.0638        # 0.053 # chris:0.057
UNIT_CELL_SIZE     = False                # used in the function
WATER_MARGIN_SURF_VALUE = 1.16#1.05 #1.19# 1.15  #This value should be changed for different unit cell sizes (1.15 for 200, 1.05 for 250)

class Membrane(molar.pdb.Pdb):
    
    def __init__(self,number_of_cells_per_side_=1,unit_cell_size_ = 200):
        
        molar.pdb.Pdb.__init__(self)
        self.surface_density= SURFACE_DENSITY
        self.step_num            = int(16) #32
        self.marching_cubes = vtk.vtkMarchingCubes()
        self.poly_data            = vtk.vtkPolyData()
        self.data                    = vtk.vtkDataObject()
        self.unit_cell_size      = unit_cell_size_
        self.form                    = "G"
        self.number_of_cells_per_side = number_of_cells_per_side_
        self.cerebroside    = False
        self.units               = []
        self.units_transed  = []
        self.points             = vtk.vtkPoints() #for visualization
        self.attempt               = 0
        self.collisioncounter  = 0
        self.total_area           = 0
        self.total_pair            = 0
        for f in FILES:
            self.units.append(molar.pdb.Pdb())
            self.units[-1].ReadFile(os.path.join(PATH , f))
        for f in FILES:
            self.units_transed.append(molar.pdb.Pdb())
            self.units_transed[-1].ReadFile(os.path.join(PATH , f))
            self.units_transed[-1].RotateX(180)
        self.pointlocator   = vtk.vtkPointLocator()
        self.pointlocator.Initialize()
        self.pointlocator.SetDataSet(vtk.vtkPolyData())
        self.pointlocator.InitPointInsertion(vtk.vtkPoints(), [0,unit_cell_size_,0,unit_cell_size_,0,unit_cell_size_])
        self.dist_dic                        = False
        self.micelle_radius             = 25.0
        self.tube_radius                 = 15.0
        self.tube_length                 = 100.0
        self.tube_hex_side             = self.unit_cell_size # also equal to the small diagonal of the rhombus of the unit-cell
        self.disc_height                  = 1.5 ## 
        self.disc_cell_height           =  110
        self.disc_radius                  = False
        self.test_dic = dict()
        self.test_dic1 = dict()
        
    def Update(self): #override
        self.side_length    = float(self.number_of_cells_per_side * self.unit_cell_size)
        self.sampling_size = self.side_length / (self.step_num) #.../(self.step_num - 1)
        self.micelle_scaled_radius  = self.micelle_radius * (2*np.pi) / self.unit_cell_size # radius of micelles scaled (normalized)
        if self.form == "tubular":
            self.sampling_size = np.sqrt(3) *  self.side_length / (self.step_num - 1)
            self.step_num_l = int (np.floor(self.tube_length / self.sampling_size) )# step number along the tube.
            # hexagonal crystal
            self.DefCryst(a=self.tube_hex_side+2,b=self.tube_hex_side+2,c=self.tube_length+5, alpha=90,beta=90,gamma=120,sGroup="P 1", zvalue=1)
        elif self.form == "disc" :
            if not self.disc_radius:
                self.disc_radius  = 0.5 * self.unit_cell_size - 3.0
            self.sampling_size = np.sqrt(3) *  self.side_length / (self.step_num - 1)
            self.sampling_size_z = self.sampling_size / 30.0
            self.step_num_l = int (np.floor(self.disc_height / self.sampling_size_z) )
            ## calculate the pbc box angle parameters
            a_=self.tube_hex_side+55.0
            b_=self.tube_hex_side+55.0
            h=self.disc_cell_height+20.0
            idea_r = 0.5*a_ ## used for unit cell parameters 
            d=idea_r / np.cos(np.radians(30))
            c_=np.sqrt( h**2 + d**2 )
            alpha_ = np.degrees ( np.arccos( idea_r /c_) )
            self.DefCryst(a=a_,b=b_,c=c_, alpha=alpha_,beta=90,gamma=120,sGroup="P 1", zvalue=1)
            
    def Func(self,_x_,_y_,_z_):
        if self.form == "G":
            return np.sin(_x_) * np.cos(_y_) + np.sin(_y_) * np.cos(_z_) + np.sin(_z_) * np.cos(_x_)
        elif self.form == "P": 
            return np.cos(_x_) + np.cos(_y_) + np.cos(_z_)
        elif self.form == "D":
            #return np.cos(_x_-_y_) *np.sin(_z_)  + np.sin(_x_+_y_) * np.cos(_z_)
            return np.cos(_x_-_y_) *np.sin(_z_)  + np.sin(_x_+_y_) * np.cos(_z_)
        elif self.form == "micelle" or self.form == "reverse-micelle":
            return np.sqrt (  min ( 
                                            (_x_-np.pi)**2    + (_y_-np.pi)**2    + (_z_-np.pi)**2, # center of the box
                                            (_x_-0.)**2         + (_y_-0.)**2         + (_z_-0.)**2,      # vertices of the box
                                            (_x_-2*np.pi)**2 + (_y_-0.)**2         + (_z_-0.)**2,
                                            (_x_-0.)**2         + (_y_-2*np.pi)**2 + (_z_-0.)**2,
                                            (_x_-0.)**2         + (_y_-0.)**2         + (_z_-2*np.pi)**2,
                                            (_x_-2*np.pi)**2 + (_y_-2*np.pi)**2 + (_z_-0.)**2,
                                            (_x_-0.)**2         + (_y_-2*np.pi)**2 + (_z_-2*np.pi)**2,
                                            (_x_-2*np.pi)**2 + (_y_-0.)**2         + (_z_-2*np.pi)**2,
                                            (_x_-2*np.pi)**2 + (_y_-2*np.pi)**2 + (_z_-2*np.pi)**2) 
                                            ) - self.micelle_scaled_radius
        elif self.form == "tubular":
            # real scale is used.
            return np.sqrt ( (_x_-0.5*self.tube_hex_side)**2         + (_y_-sqrt3o2*self.tube_hex_side)**2 ) - self.tube_radius 
        elif self.form == "disc":
            # real scale is used.
            value = (np.sqrt ( (_x_-0.5*self.unit_cell_size)**2         + (_y_-sqrt3o2*self.unit_cell_size)**2 ) - self.disc_radius  )
            if value<0 and 0.0<_z_<self.disc_height :
                return value
            elif value>0:
                return value
            else:
                 return 1.0
             
    def MarchingCubes(self):
        self.Update()
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
        self.marching_cubes.SetValue(0, 0.0)#(0,0)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()

        self.poly_data = self.marching_cubes.GetOutput(0)
        self.normals    = self.poly_data.GetPointData().GetNormals()
        self.cell_n       = self.poly_data.GetNumberOfCells()
        
    def MarchingCubes2(self):
        """ to be used for tubular form.
        """
        self.Update()
        step    =   self.sampling_size
        vol = vtk.vtkImageData()
        print "step_num: ",self.step_num,self.step_num_l 
        vol.SetDimensions(self.step_num , self.step_num , self.step_num_l )
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        for i in range(self.step_num):
            for j in range(self.step_num):
                for k in range(self.step_num_l):
                    f = self.Func(i*step,j*step,k*step)
                    vol.SetScalarComponentFromDouble(i,j,k,0,f)
        self.marching_cubes.SetInputData(vol)
        self.marching_cubes.ComputeNormalsOn()
        self.marching_cubes.SetValue(0, 0)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()

        self.poly_data = self.marching_cubes.GetOutput(0)
        self.normals    = self.poly_data.GetPointData().GetNormals()
        self.cell_n       = self.poly_data.GetNumberOfCells()
        
    def MarchingCubes3(self):
        """ to be used for disc form.
        """
        self.Update()
        step    =   self.sampling_size
        step_z =  self.sampling_size_z
        vol = vtk.vtkImageData()
        print "step_num: ",self.step_num,self.step_num_l 
        vol.SetDimensions(self.step_num , self.step_num , self.step_num_l +3) # if the disc diameter is larger than the unit_cell_size then it is cleaved off in X dimension,
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        for i in range(self.step_num):
            for j in range(self.step_num):
                for k in range(-1,self.step_num_l+2):
                    f = self.Func(i*step,j*step,k*step_z)
                    vol.SetScalarComponentFromDouble(i,j,k+1,0,f)
        self.marching_cubes.SetInputData(vol)
        self.marching_cubes.ComputeNormalsOn()
        self.marching_cubes.SetValue(0, 0)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()

        self.poly_data = self.marching_cubes.GetOutput(0)
        self.normals    = self.poly_data.GetPointData().GetNormals()
        self.cell_n       = self.poly_data.GetNumberOfCells()

    def CalcTransform(self,triangle,rotation_angle = False , x_displacement = False , y_displacement = False):
        center = np.zeros([3],dtype="f")
        trans  = vtk.vtkTransform()
        trans.PostMultiply() 
        for i in range(3):
            center = center + triangle[i]
        center = center * self.sampling_size / 3.0
        
        AreaVec = CalcAreaVec(triangle * self.sampling_size)
        AreaVal  = np.linalg.norm(AreaVec)
        if  AreaVal == 0.0:  # Avoiding NaN
            return False
        normal  = AreaVec / AreaVal
        if y_displacement:
            trans.Translate([0,y_displacement,0])
        if x_displacement:
            trans.Translate([x_displacement,0,0])
        if rotation_angle:
            trans.RotateY(rotation_angle)
        else:
            trans.RotateY(360 * rand.random())
        trans.Concatenate(molar.pdb.RotateToParallel(normal,[0,1,0]))
        
        ## test ##
#         test_point = [0,0,0]
#         trans.TransformPoint([0,1,0],test_point)
        #print test_point
#         if np.dot(test_point,normal) < 0:
#             pass
#             print np.dot(test_point,normal)
        ## ##
        
        trans.Translate(center)
        #trans.Translate([1*rand.random(),1*rand.random(),1*rand.random()])

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
        print "ignored triangles    : ", count
        print "total area                : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions                : ", self.collisioncounter
        print "average density      : ", self.total_pair / self.total_area
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
        AND with varying lengths. with some glycoseramides. glyco is between 0 and 1.
        """

#         dist_   = np.array([0.6/4 , 2.47/4 , 0.41/4 , 0.52/4],dtype=float)
#         dist_   = dist_ / dist_.sum()
#         self.dist  = np.add.accumulate(dist_)
        ########### loading component files ###########
        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.Translate([0,10,0])
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.Translate([0,10,0])
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.Translate([0,10,0])
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        #extra_fah30 = eos_conc    # extra fah30 is added to balance [Magnus compensation]
        extra_fah30 = 0.0
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - extra_fah30
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5  # ns & np 
        """self.dist is a dic() of absolute accumulative distributions."""
        self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
                       self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
                       self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
                       self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,          self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
                       self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
                       self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
                       self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
                       self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
                       self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
                       self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
                       self.chol: chol_conc,\
                       self.eos:  (1-glyco)*eos_conc,\
                       self.geos: glyco *eos_conc}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
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
            s = self.TriangleAnalysisMartini4(triangle)
            count = count + s
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
        self.PrintInfo()
        pass
    
    def MakeMartini4_1(self, glyco=1.0):#TODO

        """make with Martini coarse-grained molecule with cholesterol and fatty acid. G D P.
        WHITHOUT varying lengths. with some glycoseramides. glyco is between 0 and 1.
        """

#         dist_   = np.array([0.6/4 , 2.47/4 , 0.41/4 , 0.52/4],dtype=float)
#         dist_   = dist_ / dist_.sum()
#         self.dist  = np.add.accumulate(dist_)
        ########### loading component files ###########
        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.Translate([0,10,0])
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.Translate([0,10,0])
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.Translate([0,10,0])
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute accumulative distributions."""
#         self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,  self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
#                        self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
#                        self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
#                        self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
#                        self.chol: chol_conc,\
#                        self.eos:  (1-glyco)*eos_conc,\
#                        self.geos: glyco *eos_conc}
        #self.dist_dic={self.cerns[1]: (1-glyco)*0.38 , self.gcerns[1]: glyco*0.38 , self.fah[1]: 0.32 , self.chol: 0.30}#{self.cerns[1]: (1-glyco)*0.49 , self.gcerns[1]: glyco*0.49 , self.fah[1]: 0.28 , self.chol: 0.23}
        self.dist_dic={self.cerns[1]: (1-glyco)*0.38 , self.gcerns[1]: glyco*0.38 , self.fah[1]: 0.33 , self.chol: 0.29}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
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
            s = self.TriangleAnalysisMartini4_1(triangle)
            count = count + s
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        print "insertion attempts  : ", self.attempt 
        print "collisions          : ", self.collisioncounter
        print "average density     : ", self.total_pair / self.total_area
        print "average pair's area : ", self.total_area / self.total_pair
        self.PrintInfo()
        pass
    
    def MakeMartini5(self, glyco=1.0):

        """[Micelles].suitable for structures with lipids on one side of zero-value iso-surface.
        make with Martini coarse-grained molecule with cholesterol and fatty acid.
        AND with varying lengths. with some glycoseramides. glyco is between 0 and 1.
        suitable for structures with lipids on one side of zero-value iso-surface.
        """

#         dist_   = np.array([0.6/4 , 2.47/4 , 0.41/4 , 0.52/4],dtype=float)
#         dist_   = dist_ / dist_.sum()
#         self.dist  = np.add.accumulate(dist_)
        ########### loading component files ###########
        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.Translate([0,10,0])
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.Translate([0,10,0])
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.Translate([0,10,0])
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
            #self.cerns[-1].Show("sphere")

        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
            #self.cernp[-1].Show("sphere")

        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
            #self.gcerns[-1].Show("sphere")

        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
            #self.gcernp[-1].Show("sphere")
            
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
            #self.fah[-1].Show("sphere")
            
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute accumulative distributions."""
#         self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,  self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
#                        self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
#                        self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
#                        self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
#                        self.chol: chol_conc,\
#                        self.eos:  (1-glyco)*eos_conc,\
#                        self.geos: glyco *eos_conc}
        self.dist_dic={self.cerns[1]: (1-glyco)*0.4 , self.gcerns[1]: glyco*0.4 , self.fah[1]: 0.34 , self.chol: 0.26}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
            dist_=np.append(dist_,value)
        dist_ = dist_ / dist_.sum()
        self.dist = np.add.accumulate(dist_)
        ################################################
        self.Update()
        self.MarchingCubes()

        micelle_centers = [] 
        # center
        micelle_centers.append([0.5*self.side_length,0.5*self.side_length,0.5*self.side_length])
        # vertices
        
        for x in [0.0,self.side_length]:
            for y in [0.0,self.side_length]:
                for z in [0.0,self.side_length]:
                    point = [x,y,z] 
                    micelle_centers.append( point )
        self.micelle_centers = micelle_centers
        ##########################
        triangle_bs   = np.zeros([3,3],dtype="f") #triangle before scaling!! 
        count = 0

        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle_bs[j])
            triangle = self.sampling_size * triangle_bs  # fix the scale 
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle))
            triangle_center = TriangleCenter(triangle)
            #print triangle_center
            vertex_id = molar.pdb.FindClosestPoint( triangle_center , micelle_centers )
            micelle_center = np.array(micelle_centers[vertex_id])
            #print "\ntriangle: ",triangle , "\nmicelle_center",micelle_center
            micelle_center = micelle_center / self.sampling_size       # micelle_center and  triangle_bs should have the same scale
            triangle_bs = np.array(triangle) / self.sampling_size       # fix the facing in triangle_bs as well.
            s = self.TriangleAnalysisMartini5(triangle_bs,micelle_center)  # TriangleAnalysisMartini5 takes triangle before scaling to the right dimensions.
            count = count + 0

        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        self.PrintInfo()
        pass
    
    def MakeMartini6(self, glyco=1.0):

        """[reverse Micelles].suitable for structures with lipids on one side of zero-value iso-surface.
        make with Martini coarse-grained molecule with cholesterol and fatty acid.
        AND with varying lengths. with some glycoseramides. glyco is between 0 and 1.
        suitable for structures with lipids on one side of zero-value iso-surface.
        """

#         dist_   = np.array([0.6/4 , 2.47/4 , 0.41/4 , 0.52/4],dtype=float)
#         dist_   = dist_ / dist_.sum()
#         self.dist  = np.add.accumulate(dist_)
        ########### loading component files ###########
        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.Translate([0,10,0])
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.Translate([0,10,0])
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.Translate([0,10,0])
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
            #self.cerns[-1].Show("sphere")

        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
            #self.cernp[-1].Show("sphere")

        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
            #self.gcerns[-1].Show("sphere")

        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
            #self.gcernp[-1].Show("sphere")
            
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
            #self.fah[-1].Show("sphere")
            
        ################  fix distributions ################
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute accumulative distributions."""
#         self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,  self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
#                        self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
#                        self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
#                        self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
#                        self.chol: chol_conc,\
#                        self.eos:  (1-glyco)*eos_conc,\
#                        self.geos: glyco *eos_conc}
        self.dist_dic={self.cerns[1]: (1-glyco)*0.40 , self.gcerns[1]: glyco*0.40 , self.fah[1]: 0.3 , self.chol: 0.3}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
            dist_=np.append(dist_,value)
        dist_ = dist_ / dist_.sum()
        self.dist = np.add.accumulate(dist_)
        ################################################
        self.Update()
        self.MarchingCubes()

        micelle_centers = [] 
        # center
        micelle_centers.append([0.5*self.side_length,0.5*self.side_length,0.5*self.side_length])
        # vertices
        
        for x in [0.0,self.side_length]:
            for y in [0.0,self.side_length]:
                for z in [0.0,self.side_length]:
                    point = [x,y,z] 
                    micelle_centers.append( point )
        self.micelle_centers = micelle_centers
        #################################################
        triangle_bs   = np.zeros([3,3],dtype="f") #triangle before scaling!! 
        count = 0

        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle_bs[j])
            triangle = self.sampling_size * triangle_bs  # fix the scale 
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle))
            triangle_center = TriangleCenter(triangle)
            #print triangle_center
            vertex_id = molar.pdb.FindClosestPoint( triangle_center , micelle_centers )
            micelle_center = np.array(micelle_centers[vertex_id])
            #print "\ntriangle: ",triangle , "\nmicelle_center",micelle_center
            micelle_center = micelle_center / self.sampling_size       # micelle_center and  triangle_bs should have the same scale
            triangle_bs = np.array(triangle) / self.sampling_size       # fix the facing in triangle_bs as well.
            s = self.TriangleAnalysisMartini6(triangle_bs,micelle_center,side = -1)  # TriangleAnalysisMartini6 takes triangle before scaling to the right dimensions.
            count = count + 0

        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        self.PrintInfo()
        pass
    
    def MakeMartini7(self, glyco=1.0):
        """[tubular].suitable for structures with lipids on one side of zero-value iso-surface.
        make with Martini coarse-grained molecule with cholesterol and fatty acid.
        AND with varying lengths. with some glycoseramides. glyco is between 0 and 1.
        suitable for structures with lipids on one side of zero-value iso-surface.
        """

        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.BringToNegativeY()
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.BringToNegativeY()
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.BringToNegativeY()
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].BringToNegativeY()
            #self.cerns[-1].Show("sphere")

        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].BringToNegativeY()
            #self.cernp[-1].Show("sphere")

        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].BringToNegativeY()
            #self.gcerns[-1].Show("sphere")

        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].BringToNegativeY()
            #self.gcernp[-1].Show("sphere")
            
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].BringToNegativeY()
            #self.fah[-1].Translate([0,10,0])
            #self.fah[-1].Show("sphere")
            
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute accumulative distributions."""
#         self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,  self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
#                        self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
#                        self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
#                        self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
#                        self.chol: chol_conc,\
#                        self.eos:  (1-glyco)*eos_conc,\
#                        self.geos: glyco *eos_conc}
        self.dist_dic={self.cerns[1]: (1-glyco)*0.46 , self.gcerns[1]: glyco*0.46 , self.fah[1]: 0.33 , self.chol: 0.21}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
            dist_=np.append(dist_,value)
        dist_ = dist_ / dist_.sum()
        self.dist = np.add.accumulate(dist_)
        ################################################
        self.Update()
        self.MarchingCubes2()

        ##########################
        triangle_bs   = np.zeros([3,3],dtype="f") #triangle before scaling!! 
        count = 0
        tube_center_bs = [0.5*self.tube_hex_side / self.sampling_size   , sqrt3o2*self.tube_hex_side / self.sampling_size ,0.0]
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle_bs[j])
            triangle = self.sampling_size * triangle_bs  # fix the scale 
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle))
            s = self.TriangleAnalysisMartini5(triangle_bs,tube_center_bs )  # TriangleAnalysisMartini5 takes triangle before scaling to the right dimensions.
            count = count + 0
        
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
        self.PrintInfo()
        pass
    
    def MakeMartini8(self, glyco=1.0):
        """ Vesicles. self.form is ineffective. 
        """

        self.form = "reverse-micelle"
        self.MakeMartini6(glyco)
        #self.Show()
        self.form = "micelle"
        self.micelle_radius = self.micelle_radius + 0
        self.MakeMartini5(glyco)
        self.form = "vesicle"  # for adding the water

    def MakeMartini9(self, glyco=1.0):
        """[disc of vesicle].suitable for structures with lipids on one side of zero-value iso-surface.
        make with Martini coarse-grained molecule with cholesterol and fatty acid.
        AND with varying lengths. with some glycoseramides. glyco is between 0 and 1.
        suitable for structures with lipids on one side of zero-value iso-surface.
        """

        self.cerns       = []
        self.cernp       = []
        self.gcerns     = []
        self.gcernp     = []
        self.fah           = []
        self.chol         = molar.pdb.Pdb()
        self.chol.ReadFileGMX(os.path.join(PATH , CHOL))
        self.chol.BringToCenter()
        self.chol.BringToNegativeY()
        self.eos        = molar.pdb.Pdb()
        self.eos.ReadFileGMX(os.path.join(PATH , EOS))
        self.eos.BringToCenter()
        self.eos.BringToNegativeY()
        self.geos        = molar.pdb.Pdb()
        self.geos.ReadFileGMX(os.path.join(PATH , GEOS))
        self.geos.BringToCenter()
        self.geos.BringToNegativeY()
        for file_name in CER_NS:
            self.cerns.append(molar.pdb.Pdb())
            self.cerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cerns[-1].BringToCenter()
            self.cerns[-1].Translate([0,10,0])
            #self.cerns[-1].Show("sphere")

        for file_name in CER_NP:
            self.cernp.append(molar.pdb.Pdb())
            self.cernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.cernp[-1].BringToCenter()
            self.cernp[-1].Translate([0,10,0])
            #self.cernp[-1].Show("sphere")

        for file_name in GCER_NS:
            self.gcerns.append(molar.pdb.Pdb())     
            self.gcerns[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcerns[-1].BringToCenter() 
            self.gcerns[-1].Translate([0,10,0])
            #self.gcerns[-1].Show("sphere")

        for file_name in GCER_NP:
            self.gcernp.append(molar.pdb.Pdb())
            self.gcernp[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.gcernp[-1].BringToCenter()
            self.gcernp[-1].Translate([0,10,0])
            #self.gcernp[-1].Show("sphere")
            
        for file_name in FAH:
            self.fah.append(molar.pdb.Pdb())
            self.fah[-1].ReadFileGMX(os.path.join(PATH , file_name))
            self.fah[-1].BringToCenter()
            self.fah[-1].Translate([0,10,0])
            #self.fah[-1].Show("sphere")
            
        ################  fix distributions ############
        ### overall concentrations: ###
        eos_conc    = 0.05             # overall ceramide eos concentration  
        extra_fah30 = eos_conc    # extra fah30 is added to balance
        cer_conc = 0.33 - eos_conc
        fah_conc = 0.33 - eos_conc
        chol_conc= 0.33
        d        = [0.6/4 , 2.47/4 , 0.41/4 , 0.52/4]
        h        = 0.5
        """self.dist is a dic() of absolute accumulative distributions."""
#         self.dist_dic={self.cerns[0]: h*(1-glyco)*d[0] * cer_conc,   self.cerns[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cerns[2]: h*(1-glyco)*d[2] * cer_conc,          self.cerns[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.cernp[0]: h*(1-glyco)*d[0] * cer_conc,          self.cernp[1]:  h*(1-glyco)*d[1] * cer_conc,\
#                        self.cernp[2]: h*(1-glyco)*d[2] * cer_conc,  self.cernp[3]:  h*(1-glyco)*d[3] * cer_conc,\
#                        self.gcerns[0]:h*glyco*d[0] * cer_conc,      self.gcerns[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcerns[2]:h*glyco*d[2] * cer_conc,      self.gcerns[3]: h*glyco*d[3] * cer_conc,\
#                        self.gcernp[0]:h*glyco*d[0] * cer_conc,      self.gcernp[1]: h*glyco*d[1] * cer_conc,\
#                        self.gcernp[2]:h*glyco*d[2] * cer_conc,      self.gcernp[3]: h*glyco*d[3] * cer_conc,\
#                        self.fah[0]:   d[0] * fah_conc,              self.fah[1]:    d[1] * fah_conc,\
#                        self.fah[2]:   d[2] * fah_conc,              self.fah[3]:    d[3] * fah_conc + extra_fah30,\
#                        self.chol: chol_conc,\
#                        self.eos:  (1-glyco)*eos_conc,\
#                        self.geos: glyco *eos_conc}
        self.dist_dic={self.cerns[1]: (1-glyco)*0.36 , self.gcerns[1]: glyco*0.36, self.fah[1]: 0.31 , self.chol: 0.33}
        self.componants = []
        self.componants_inv = []
        dist_           = np.array([])
        print "The aiming relative concentrations are:"
        for key, value in self.dist_dic.iteritems():
            self.componants.append(key)
            self.componants_inv.append(key.MakeCopy())
            self.componants_inv[-1].RotateX(180)
            print "%10s :%10f" % (key.molecules[-1].name,value)
            dist_=np.append(dist_,value)
        dist_ = dist_ / dist_.sum()
        self.dist = np.add.accumulate(dist_)
        ################################################
        self.Update()
        self.MarchingCubes3()

        ##########################
        triangle_bs   = np.zeros([3,3],dtype="f") #triangle before scaling!! 
        count = 0
        #disc_center_bs = [0.5*self.unit_cell_size ,  sqrt3o2*self.unit_cell_size,0.5*self.disc_height]
        for i in range(self.cell_n):
            molar.pdb.update_progress(float(i)/self.cell_n)                
            cell        = self.poly_data.GetCell(i)
            
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                self.poly_data.GetPoints().GetPoint(point_id,triangle_bs[j])
            triangle = self.sampling_size * triangle_bs  # fix the scale 
            self.total_area = self.total_area + np.linalg.norm(CalcAreaVec(triangle))
            s = self.TriangleAnalysisMartini4_1(triangle_bs )  # TriangleAnalysisMartini5 takes triangle before scaling to the right dimensions.
            count = count + 0
        
        self.Update()        # to modify the crystal parameters with the real disc radius
        print "ignored triangles   : ", count
        print "total area          : ", self.total_area ,"  angstrom squared"
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
    
    def TriangleAnalysisMartini4(self,t): # recursive
        """ to be used with MakeMartini4()
        """
        
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size))
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num      = self.surface_density * area
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            j       = self.ProduceIndex()
            new_mol = self.componants[j]
            for y in [0,2,4]:
                for d in [0.0,2.0,4.0]:
                    for i in range(0,360,90):
                        trans   = self.CalcTransform(t,i,d,y)
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 2.0) #cutoff = 1.2
                        if success:
                            break
                    if success:
                        break
                if success:
                    break
            else:
                self.collisioncounter += 1
            self.attempt += 1
            ## insert the opposite molecule.
            j       = self.ProduceIndex()
            new_mol = self.componants_inv[j]
            for y in range(3):
                for d in range(10):
                    for i in range(0,360,10):
                        trans   = self.CalcTransform(t,i,d,-1*y)
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 2.0) 
                        if success:
                            break
                    if success:
                        break
                if success:
                    break
            else:
                self.collisioncounter += 1
            self.attempt += 1
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini4(np.array([t[0],t[2],MidPoint(t[0],t[1])]))
                    self.TriangleAnalysisMartini4(np.array([t[1],t[2],MidPoint(t[0],t[1])]))
                else:
                    self.TriangleAnalysisMartini4(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini4(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini4(np.array([t[0],t[2],MidPoint(t[2],t[1])]))
                    self.TriangleAnalysisMartini4(np.array([t[0],t[1],MidPoint(t[2],t[1])]))
                else:
                    self.TriangleAnalysisMartini4(np.array([t[0],t[1],MidPoint(t[0],t[2])]))
                    self.TriangleAnalysisMartini4(np.array([t[1],t[2],MidPoint(t[0],t[2])]))
        else: # rand_point_num == 0
            return 1
        return 0             

    def TriangleAnalysisMartini4_1(self,t): # recursive
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size)) 
        point_num      = 1
        #rand_point_num =  np.random.poisson(point_num)
        if (np.floor(point_num) == 1):
            j       = self.ProduceIndex()
            new_mol = self.componants[j]
            for y in [0,2]:
                for d in [0.0,2.0,4.0]:
                    for i in range(0,360,45):
                        trans   = self.CalcTransform(triangle=t , rotation_angle=i , x_displacement=d , y_displacement=y)
                        ###
                        if trans == False:
                            success = False
                            continue
                        ###
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 3.0) #5.0
                        if success:
                            break
                    if success:
                        break
                if success:
                    break  
            ## insert the opposite molecule.
            j       = self.ProduceIndex()
            new_mol = self.componants_inv[j]
            for y in [0,2]:
                for d in [0.0,2.0,4.0]:
                    for i in range(0,360,45):
                        trans   = self.CalcTransform(triangle=t , rotation_angle=i , x_displacement=d , y_displacement=-1*y)
                        ###
                        if trans == False:
                            success = False
                            continue
                        ###
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 3.0) #5.0
                        if success:
                            break
                    if success:
                        break
                if success:
                    break   
            else:
                self.collisioncounter += 1
            self.attempt += 1
            
        ## test ##
        if  success:
                name = new_mol.molecules[0].name
                if name in self.test_dic:
                    self.test_dic[name]+=1
                else:
                    self.test_dic[name]=1
        else:
                name = new_mol.molecules[0].name
                if name in self.test_dic1:
                    self.test_dic1[name]+=1
                else:
                    self.test_dic1[name]=1
        ##    ##
        return success

    def TriangleAnalysisMartini5(self,t,micelle_center,side = +1): # recursive
        
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size)) 
        l1   = np.linalg.norm(t[0] - t[1])
        l2   = np.linalg.norm(t[1] - t[2])
        l3   = np.linalg.norm(t[0] - t[2])
        point_num = 1
        #point_num      = self.surface_density * area 
        #rand_point_num =  np.random.poisson(point_num)
        t = TriangleFaceOrderFix(t, micelle_center,side)
        if (np.floor(point_num) == 1):
            j       = self.ProduceIndex()
            new_mol = self.componants[j]
            for y in [0]:
                for d in [0.0,2.0,4.0]:
                    for i in range(0,360,90):
                        trans   = self.CalcTransform(triangle=t , rotation_angle=i , x_displacement=d , y_displacement=y)
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 5.0)
                        if success:
                            break
                    if success:
                        break
                if success:
                    break
            else:
                self.collisioncounter += 1
            self.attempt += 1
            
        elif (np.floor(point_num) > 1):
            if l1 > l2:
                if l1 > l3:
                    self.TriangleAnalysisMartini5(np.array([t[0],t[2],MidPoint(t[0],t[1])]),micelle_center )
                    self.TriangleAnalysisMartini5(np.array([t[1],t[2],MidPoint(t[0],t[1])]),micelle_center )
                else:
                    self.TriangleAnalysisMartini5(np.array([t[0],t[1],MidPoint(t[0],t[2])]),micelle_center )
                    self.TriangleAnalysisMartini5(np.array([t[1],t[2],MidPoint(t[0],t[2])]),micelle_center )
            else:
                if l2 > l3:
                    self.TriangleAnalysisMartini5(np.array([t[0],t[2],MidPoint(t[2],t[1])]),micelle_center )
                    self.TriangleAnalysisMartini5(np.array([t[0],t[1],MidPoint(t[2],t[1])]),micelle_center )
                else:
                    self.TriangleAnalysisMartini5(np.array([t[0],t[1],MidPoint(t[0],t[2])]),micelle_center )
                    self.TriangleAnalysisMartini5(np.array([t[1],t[2],MidPoint(t[0],t[2])]),micelle_center )
        else: # rand_point_num == 0
            return 1
        return 0             
        
    def TriangleAnalysisMartini6(self,t,micelle_center,side = +1): # recursive
        area = np.linalg.norm(CalcAreaVec(t * self.sampling_size)) 
        point_num      = 1
        #rand_point_num =  np.random.poisson(point_num)
        t = TriangleFaceOrderFix(t, micelle_center,side)
        if (np.floor(point_num) == 1):
            j       = self.ProduceIndex()
            new_mol = self.componants[j]
            for y in [0]:
                for d in [0.0,2.0,4.0]:
                    for i in range(0,360,90):
                        trans   = self.CalcTransform(triangle=t , rotation_angle=i , x_displacement=d , y_displacement=y)
                        success = self.CatTransformedCautious(new_mol, trans, self.pointlocator, cutoff = 5.0) #5.0
                        if success:
                            break
                    if success:
                        break
                if success:
                    break  
            else:
                self.collisioncounter += 1
            self.attempt += 1
            
        ## test ##
        if  success:
                name = new_mol.molecules[0].name
                if name in self.test_dic:
                    self.test_dic[name]+=1
                else:
                    self.test_dic[name]=1
        else:
                name = new_mol.molecules[0].name
                if name in self.test_dic1:
                    self.test_dic1[name]+=1
                else:
                    self.test_dic1[name]=1
        ##    ##
        return success

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
    
    def AddWater(self,margin=False):
        water_box  = molar.pdb.Pdb()
        water_box.ReadFileGMX(os.path.join(PATH , "water_antifreeze.pdb"))
        box_length = water_box.Bounds()[0,1] - water_box.Bounds()[0,0]
        disc_max_z  = self.Bounds()[2][1] ## *important* before adding any water
        trans = vtk.vtkTransform()
        trans.PostMultiply()
        if self.form in ["G","micelle","reverse-micelle","P","D","vesicle"]:
            for x in np.arange(0,self.unit_cell_size,box_length):
                for y in np.arange(0,self.unit_cell_size,box_length):
                    for z in np.arange(0,self.unit_cell_size,box_length):
                        trans.Identity()
                        trans.Translate([x,y,z])
                        self.CatTransformed(water_box, trans)
                pass
            pass
        elif self.form in ["tubular"]:
            for x in np.arange(0,self.unit_cell_size,box_length):
                for y in np.arange(0,np.sqrt(3)*self.unit_cell_size,box_length):
                    for z in np.arange(0,self.tube_length,box_length):
                        trans.Identity()
                        trans.Translate([x,y,z])
                        self.CatTransformed(water_box, trans)
                pass
            pass
        elif self.form in ["disc"]:
            for x in np.arange(0,self.unit_cell_size,box_length):
                for y in np.arange(0,np.sqrt(3)*self.unit_cell_size,box_length):
                    for z in np.arange(0 , self.disc_cell_height+box_length  , box_length):
                        trans.Identity()
                        trans.Translate([x,y,z])
                        self.CatTransformed(water_box, trans)
                pass
            pass
        ## ############  remove redundant #####################
        ###################################################
        if self.form in ["G","P","D"]:
            print "Water removal ... "
            new_mol_list=[]
            t   = len(self.molecules)
            i=0
            for mol in self.molecules:
                mol_keep = True
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                if mol.name.replace(" ","") in ["W","WF"]:
                    pos = mol.Center()
                    pos_scaled = ( pos / self.unit_cell_size ) * 2 * np.pi 
                    value = abs(self.Func( pos_scaled[0],pos_scaled[1],pos_scaled[2] ) )
                    if ( value <= WATER_MARGIN_SURF_VALUE )\
                        or pos[0]<0 or pos[0] > self.unit_cell_size\
                        or pos[1]<0 or pos[1] > self.unit_cell_size\
                        or pos[2]<0 or pos[2] > self.unit_cell_size :
                            mol_keep = False
                            pass
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ##########################
        elif self.form == "micelle":
            self.Prune()
            if not margin:
                margin=30
            print "making backward links ..."
            self.MakeBackwardLinks()
            t = len(self.molecules)
            i=0
            print "Water removal ... "
            new_mol_list=[]
            for mol in self.molecules:
                mol_keep = True
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                for micelle_center in self.micelle_centers:
                    if mol.name.replace(" ","") in ["W","WF"]:
                        pos = mol.Center()
                        d = np.sqrt( (pos[0]-micelle_center[0])**2 
                                 + (pos[1]-micelle_center[1])**2    
                                 + (pos[2]-micelle_center[2])**2 )
                        if ( d <= ( self.micelle_radius  + margin )
                             or pos[0]<0 or pos[0] > self.unit_cell_size
                             or pos[1]<0 or pos[1] > self.unit_cell_size
                             or pos[2]<0 or pos[2] > self.unit_cell_size ):
                            mol_keep = False
                            pass
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ##############################   
        elif self.form == "reverse-micelle": 
            self.Prune()
            if not margin:
                margin=25
            print "making backward links ..."
            self.MakeBackwardLinks()
            t = len(self.molecules)
            i=0
            print "Water removal ... "
            new_mol_list=[]
            for mol in self.molecules:
                mol_keep = False
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                for micelle_center in self.micelle_centers:
                    if mol.name.replace(" ","") in ["W","WF"]:
                        pos = mol.Center()
                        d = np.sqrt( (pos[0]-micelle_center[0])**2 
                                 + (pos[1]-micelle_center[1])**2    
                                 + (pos[2]-micelle_center[2])**2 )
                        if ( d < ( self.micelle_radius  - margin )
                             and pos[0]>0 and pos[0] < self.unit_cell_size
                             and pos[1]>0 and pos[1] < self.unit_cell_size
                             and pos[2]>0 and pos[2] < self.unit_cell_size ):
                            mol_keep = True
                    else:
                        mol_keep = True  # keep non water
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ##############################   
        elif self.form == "tubular":                     
            self.Prune()
            a = 0.5 * self.side_length
            b = np.sqrt(3) * 0.5 *  self.side_length
            if not margin:
                margin=5
            print "making backward links ..."
            self.MakeBackwardLinks()
            t = len(self.molecules)
            i=0
            print "Water removal ... "
            new_mol_list=[]
            for mol in self.molecules:
                mol_keep = True
                pos = mol.Center()
                d = np.sqrt( (pos[0] - a)**2 + (pos[1] - b)**2 )
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                if mol.name.replace(" ","") in ["W","WF"]:
                    if pos[1] < ( sqrt3 * pos[0] - sqrt3 * a )\
                        or pos[1] < ( -sqrt3 * pos[0] + b )\
                        or pos[1] > ( sqrt3 * pos[0] + b )\
                        or pos[1] > ( -sqrt3 * pos[0] + sqrt3 * self.side_length + b )\
                        or pos[2] > self.tube_length\
                        or d < self.tube_radius + margin:
                            mol_keep = False
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ##############################   
        elif self.form == "vesicle":
            self.Prune()
            if not margin:
                margin=23
            print "making backward links ..."
            self.MakeBackwardLinks()
            t = len(self.molecules)
            i=0
            print "Water removal ...  "
            new_mol_list=[]
            for mol in self.molecules:
                mol_keep = True
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                for micelle_center in self.micelle_centers:
                    if mol.name.replace(" ","") in ["W","WF"]:
                        pos = mol.Center()
                        d = np.sqrt( (pos[0]-micelle_center[0])**2 
                                 + (pos[1]-micelle_center[1])**2    
                                 + (pos[2]-micelle_center[2])**2 )
                        if ( ( d > ( self.micelle_radius  - margin ) and  d < ( self.micelle_radius  + margin)  )
                             or pos[0]<0 or pos[0] > self.unit_cell_size
                             or pos[1]<0 or pos[1] > self.unit_cell_size
                             or pos[2]<0 or pos[2] > self.unit_cell_size ):
                            mol_keep = False
                    else:
                        mol_keep = True  # keep non water
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ##############################
        elif self.form == "disc":
            self.Prune()
            a = 0.5 * self.side_length
            b = np.sqrt(3) * 0.5 *  self.side_length
            if not margin:
                margin   =25
            margin_z=5
            print "making backward links ..."
            self.MakeBackwardLinks()
            t = len(self.molecules)
            i=0
            print "Water removal ... "
            new_mol_list=[]
            for mol in self.molecules:
                mol_keep = True
                pos = mol.Center()
                d = np.sqrt( (pos[0] - a)**2 + (pos[1] - b)**2 )
                i+=1
                molar.pdb.update_progress(float(i)/(t))
                if mol.name.replace(" ","") in ["W","WF"]:
                    if pos[1] < ( sqrt3 * pos[0] - sqrt3 * a )\
                        or pos[1] < ( -sqrt3 * pos[0] + b )\
                        or pos[1] > ( sqrt3 * pos[0] + b )\
                        or pos[1] > ( -sqrt3 * pos[0] + sqrt3 * self.side_length + b )\
                        or pos[2] > self.disc_cell_height\
                        or (d < self.disc_radius + margin and pos[2] < disc_max_z + margin_z):
                            mol_keep = False
                if mol_keep:
                    new_mol_list.append(mol)
            self.molecules=new_mol_list
        ############################## 
            
    def Show(self):
        ## override
        self.Update()
        if self.form in ["G","P","D"]:
            self.SurfMarchingCubes()
            molar.pdb.Pdb.Show(self,"sphere",ext_actors=[self.surf_actor])
        elif self.form == "micelle":
            micelle_actor = self.VisVertices()
            molar.pdb.Pdb.Show(self,"sphere",ext_actors=[micelle_actor])
        else:
            micelle_actor = self.VisVertices()
            molar.pdb.Pdb.Show(self,"sphere",ext_actors=[micelle_actor])
        pass
    
    def VisVertices(self,show=False,file_name = False):
        actor     = vtk.vtkActor()
        #glyph3d   = vtk.vtkGlyph3D()
        #point_s   = vtk.vtkPointSource()
        scaleTrans=vtk.vtkTransform()
        scaleTrans.Scale(self.sampling_size,self.sampling_size,self.sampling_size)
        mapper    = vtk.vtkPolyDataMapper()
        window    = vtk.vtkRenderWindow()
        renderer  = vtk.vtkRenderer()
        interactor= vtk.vtkRenderWindowInteractor()
        interactor.Initialize()
        #poly_data = vtk.vtkPolyData()
        #point_s.SetNumberOfPoints(1)
        interactor.SetRenderWindow(window)
        #poly_data.SetPoints(self.points)
        #glyph3d.SetSourceConnection(point_s.GetOutputPort())
        #glyph3d.SetInputData(poly_data)
        #mapper.SetInputConnection(glyph3d.GetOutputPort())
        mapper.SetInputConnection(self.marching_cubes.GetOutputPort())
        mapper.ScalarVisibilityOff()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0,1,1)
        print "yes"
        if show:
            window.AddRenderer(renderer)
            renderer.AddActor(actor)
            renderer.SetBackground(0.5,0.5,0.5)
            renderer.ResetCamera()
            window.Render()
            interactor.Start()
        if file_name:
            obj = vtk.vtkOBJExporter()
            obj.SetInput(window)
            obj.SetFilePrefix(file_name)
            obj.Write()
        actor.SetUserTransform(scaleTrans)
        
        return actor
    
    def ProduceIndex(self):
        r        = rand.random()
        for i in range(len(self.dist)):
            if r < self.dist[i]:
                break
        return i

    def WriteOnFileGMX(self,file_name,make_TER=True): #override
        """ To fit the system inside the periodic box for VMD visualiziation.
        """
        if self.form in ["tubular","reverse_tubular","disc"]:
            self.BringToCenter()
            self.Translate( [ 0.25*self.tube_hex_side , np.sqrt(3)*0.25*self.tube_hex_side ,  0.5*self.tube_length ] )
            pass
        molar.pdb.Pdb.WriteOnFileGMX(self,file_name,make_TER)
        pass

##########################################################

def Func_aux(pos):
        return np.sin(pos[0]) * np.cos(pos[1]) + np.sin(pos[1]) * np.cos(pos[2]) + np.sin(pos[2]) * np.cos(pos[0])

def DelFunc(atom):
        """ works only for the box size = 200 because of WATER_MARGIN_SURF_VALUE
            this delete function can only have one argument as atom.
        """
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
    return c

def Distance(v1,v2):
    return np.linalg.norm(v1-v2)

def TriangleFaceOrderFix(triangle_ , point_ , side = +1):
    """makes sure that point_ stands on one side of triangle_ , otherwise, it swaps two vertices of the triangle.
    side determines the direction an is +1 or -1. 
    """
    triangle = np.array( triangle_ )
    point    = np.array( point_ )
    v1        = triangle[1] - triangle[0]
    v2        = triangle[2] - triangle[0]
    v3        = point - triangle[0]
    if side * np.dot ( np.cross(v1,v2) , v3) < 0:
        return triangle
    else: # swap two vertices
        c = np.copy(triangle[1])
        triangle[1] = triangle[2]
        triangle[2] = c
        return triangle
    
def TriangleCenter(np_triangle):
    center = np.zeros(3,dtype=float)
    i=0
    for p in np_triangle:
        center += p
        i+=1
    out = center / float(i)
    return out

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
