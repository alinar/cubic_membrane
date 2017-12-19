import vtk
import numpy as np
import math
import molar.pdb as pdb
import molar.pdb_viewer
import matplotlib.pyplot as plt
import MDAnalysis
import sys

ANG2="( A\\v{0.65}\h{-0.5}\z{0.6}o  -2\\v{}\h{}\z{} )"
ANG1="( A\\v{0.65}\h{-0.5}\z{0.6}o  -1\\v{}\h{}\z{} )"
            
class ExtractSurf():
    """Extracts the Surface mesh from the pdb file and analyze the surface.
    """
    def __init__(self , file_in , traj_file = False , steps_input = 5):
        self.color_series = vtk.vtkColorSeries()
        self.color_series.SetColorScheme( self.color_series.BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_9)   #set the color theme
        
        if not traj_file:
            self.traj_mode = False
            if file_in[-3:] == "pdb":
                self.pdb = pdb.Pdb()
                self.pdb.ReadFileGMX(file_in)
                #self.pdb.BringToPositiveCorner()
                self.marching_cubes=False
                self.univ = False
            elif file_in[-3:] == "gro":
                self.univ = MDAnalysis.Universe(file_in)
        else:
            self.traj_mode = True
            self.univ        = MDAnalysis.Universe(file_in,traj_file)
            self.univ_gro = MDAnalysis.Universe(file_in) #only load gro file to read coordinates and velocities from. not trajectory
            self.steps = steps_input  # number of steps to consider in the analysis in the whole trajectory file.
        
        self.gauss_mean=True # in gaussian curvature case = True / in mean curvature  case = False
        self.test = []
        self.test_lines    = vtk.vtkCellArray()
        self.test_points = vtk.vtkPoints()
        
    def GaugeFunction(self,point3_in):
        """ returning the distance between the point3_in and the closest point in self.points_periodic.
        """
        point_id = self.pointlocator.FindClosestPoint(point3_in)
        pos = self.points_periodic.GetPoint(point_id)
        return  math.sqrt(  vtk.vtkMath.Distance2BetweenPoints(point3_in, pos)  )
   
    def DigestPdb(self):
        """   add atoms too self and 8 neighboring cubes. points except for water ###
        """
        self.pdb_water_head=pdb.Pdb()  ## used for visualizations
        self.pdb_water_head.AddMolecule()
        if self.pdb.cryst:
            self.cryst = self.pdb.cryst
            print "unit cell size: ",self.cryst
        else:
            print "Please provide the unit-cell size in self.cryst"
        num_atoms = self.pdb.NumOfAtoms()
        c = 0
        self.points_periodic = vtk.vtkPoints()
        self.points               =  vtk.vtkPoints()
        print "Digesting the input pdb ..."
        for atom in self.pdb:
            if atom.name.replace(" ","") not in ["W","WF","ROH","H1","H2","H3"] and\
               ( atom.GetMolNameGMX()[0]!="G" or atom.name.replace(" ","") not in ["C1","C2","C3"] ) :
                self.points.InsertNextPoint(atom.pos)
                for q in [-self.cryst[0],0,self.cryst[0]]:
                    for p in [-self.cryst[1],0,self.cryst[1]]:
                        for r in [-self.cryst[2],0,self.cryst[2]]:
                            pos_aux = atom.pos + np.array([q,p,r])
                            self.points_periodic.InsertNextPoint(pos_aux)
                self.pdb_water_head.molecules[-1].AddAtom(atom.line)
            c+=1
            pdb.update_progress(float(c)/num_atoms)
        print "\n"
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.points_periodic)
        ### initialize pointlocator ###
        self.pointlocator = vtk.vtkPointLocator()
        self.pointlocator.SetDataSet( self.polydata )
        self.pointlocator.SetNumberOfPointsPerBucket(10)
        self.pointlocator.BuildLocator()
        
    def DigestGro(self):
        """  same as DigestPdb but for gro files
        """
        self.Atom_Select()
        self.pdb_water_head = False
        self.cryst   = self.univ.dimensions[0:3]
        num_atoms = len( self.selected_atoms )
        c = 0
        self.points_periodic = vtk.vtkPoints()
        self.points               =  vtk.vtkPoints()
        print "Digesting the input gro ..."
        for atom in self.selected_atoms:
            self.points.InsertNextPoint(atom.position)
            for q in [-self.cryst[0],0,self.cryst[0]]:
                for p in [-self.cryst[1],0,self.cryst[1]]:
                    for r in [-self.cryst[2],0,self.cryst[2]]:
                        pos_aux = atom.position + np.array([q,p,r])
                        self.points_periodic.InsertNextPoint(pos_aux)
            c+=1
            pdb.update_progress(float(c)/num_atoms)
        print "\n"
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.points_periodic)
        ### initialize pointlocator ###
        self.pointlocator = vtk.vtkPointLocator()
        self.pointlocator.SetDataSet( self.polydata )
        self.pointlocator.SetNumberOfPointsPerBucket(10)
        self.pointlocator.BuildLocator()
        
    def Atom_Select(self):
        """selecting any atom except for the water and the head-groups. 
        """
        self.selected_atoms = self.univ.select_atoms("(not name W WF ROH H1 H2 H3) and (not resname \"G*\" or not name C1 C2 C3)") 
        #print "number of atoms: ",len(self.selected_atoms)
        
    def DigestTraj(self):
        """the same function as DigestPdb just used with MDAnalysis and Gromacs files within a time step of trajectory"""
        self.points_periodic = vtk.vtkPoints()
        self.points                = vtk.vtkPoints()
        self.cryst   = self.univ.dimensions[0:3]
        for atom in self.selected_atoms:
            self.points.InsertNextPoint(atom.pos)
            for q in [-self.cryst[0],0,self.cryst[0]]:
                for p in [-self.cryst[1],0,self.cryst[1]]:
                     for r in [-self.cryst[2],0,self.cryst[2]]:
                        pos_aux = atom.position + np.array([q,p,r])
                        self.points_periodic.InsertNextPoint(pos_aux)
                
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(self.points_periodic)
        ### initialize pointlocator ###
        self.pointlocator = vtk.vtkPointLocator()
        self.pointlocator.SetDataSet( polydata )
        self.pointlocator.SetNumberOfPointsPerBucket(10)
        self.pointlocator.BuildLocator()

    def MolCurvatures(self,resname,xvg_file=False):
        """to be used ([?]after selection of a time step in self.univ and) after running GaussianCurvature or MeanCurvature 
            resname is the molecule type...
            it writes an xvg histogram file with the distribution of the molecules in different curvature areas at the last frame of 
            self.uniiv  .
        """
        ts = self.univ.trajectory[-1] # go to the last frame.
        whole_atoms = self.univ.select_atoms("resname "+resname) 
        resid_set       = set()
        point_data = self.curvature_polydata.GetPointData()
        data_array = point_data.GetScalars()
        
        for atom in whole_atoms:
            resid_set.add(atom.resid)
        # resid is needed to be able to select the atoms in a single molecule.
        
        curvature_value_list = []
        
        for resid_no in resid_set: ## == each molecule in selection
            molecule_atoms = whole_atoms.select_atoms("resid "+str(resid_no)) 
            t = np.array([0.0,0.0,0.0])
            n = 0
            for atom in molecule_atoms:
                t=t+atom.position
                n += 1
            mol_pos = t / n
            ### find the triangle:
            point_id = self.triangles_centers_pointlocator.FindClosestPoint(mol_pos)
            
            value=data_array.GetTuple1(point_id)
            if abs(value) < 0.1 and abs(value)>1e-6 :
                curvature_value_list.append(value)
            pass
        if not xvg_file:
            xvg_file=resname+".xvg"
        WriteHistXVG( curvature_value_list  ,bins=100, titile=resname , file_name=xvg_file , x_lable=ANG1 )
    
    def InitMarchingCubes(self,verbose=False):
        """Use marching cubes to creat the surface polydata (self.marching_cubes.GetOutput).
        """
        self.marching_cubes = vtk.vtkMarchingCubes()
        self.step_num = int(32) #32
        #bound = np.array(self.pdb.Bounds())
        step3    = [0.0,0.0,0.0]
        for i in range(3):
            step3[i] = self.cryst[i] / self.step_num
            #step3[i] = (bound[i][1] - bound[i][0]) / self.step_num
        vol = vtk.vtkImageData()
        vol.SetDimensions(self.step_num+1 , self.step_num+1 , self.step_num+1 )
        vol.AllocateScalars(vtk.VTK_DOUBLE,int(1));
        vol.SetOrigin([0.0,0.0,0.0])
        vol.SetSpacing([step3[0],step3[1],step3[2]])
        c=0
        if verbose:
            print "Creating the mesh using marching cubes ..."
        total = (self.step_num+1)**3
        for i in range(self.step_num+1):
            for j in range(self.step_num+1):
                for k in range(self.step_num+1):
                    f = self.GaugeFunction( [i*step3[0],j*step3[1],k*step3[2]] )
                    #print f
                    vol.SetScalarComponentFromDouble(i,j,k,0,f)
                    if verbose:
                        pdb.update_progress(float(c)/total)
                    c+=1
        self.marching_cubes.SetInputData(vol)
        self.marching_cubes.ComputeNormalsOn()
        self.marching_cubes.SetValue(0, 4.0)
        self.marching_cubes.Update()
        self.marching_cubes.UpdateInformation()
        self.sampling_size=step3
        
    def SurfaceExtractDelaunay3D(self):
        """[under construction] The same functionality as self.InitMarchingCubes() which is extraction of the surface given the points representing
        atoms. this uses vtkDelaunay3D to form set of tetrahedrons (vtkUnstructuredGrid) and then using vtkGeometryFilter 
        to find the surface mesh as vtkPolyData. TODO
        """
        
        points_vtkPolyData    = vtk.vtkPolyData()
        points_vtkPolyData.SetPoints(self.points)
        
        delaunay3D = vtk.vtkDelaunay3D()
        delaunay3D.SetInputData(points_vtkPolyData)
        delaunay3D.SetAlpha(5)

        geometry_filter = vtk.vtkGeometryFilter()
        geometry_filter.SetInputConnection(delaunay3D.GetOutputPort())
        ShowVtkAlgorithm(geometry_filter)
  
    def SmoothSurf(self):
        """To be run after InitMarchingCubes.
        """
        smoothFilter = vtk.vtkSmoothPolyDataFilter()
        #smoothFilter.BoundarySmoothingOff()
        smoothFilter.SetInputConnection(self.marching_cubes.GetOutputPort())
        smoothFilter.SetNumberOfIterations(40);
        smoothFilter.SetRelaxationFactor(0.1);
        smoothFilter.FeatureEdgeSmoothingOff();
        smoothFilter.BoundarySmoothingOn();
        smoothFilter.Update();
        #print "iterations for smoothing : ", smoothFilter.GetNumberOfIterations() 
 
        ## Update normals on newly smoothed polydata
        normalGenerator = vtk.vtkPolyDataNormals()
        normalGenerator.SetInputConnection(smoothFilter.GetOutputPort());
        normalGenerator.ComputePointNormalsOn();
        normalGenerator.ComputeCellNormalsOn();
        normalGenerator.Update();
        
        self.normal_generator_smooth = normalGenerator
                    
    def Extract(self,verbose=False,gaussian_polydata_file_out=False,mean_polydata_file_out=False):
        """extract the surface. measure the curvature variation in time."""
        if not self.traj_mode: # input was not trajectory 
            if not self.univ: # with pdb file
                print "input: pdb"
                self.DigestPdb()
            else:  #with gro file
                print "input: gro without trajectory"
                self.DigestGro()
            self.InitMarchingCubes()
            self.SmoothSurf()
            self.GaussianCurvature(verbose)
            if gaussian_polydata_file_out:
                SavePolyData(self.curvature_polydata , gaussian_polydata_file_out)
            self.MeanCurvature(verbose)
            if mean_polydata_file_out:
                SavePolyData(self.curvature_polydata , mean_polydata_file_out)
                
        else: 
            print "input: gro with trajectory"
            total_gaussian_array = []
            total_mean_array      = []
            time_array                 = []
            trj_step_no                = len(self.univ.trajectory) 
            ds = int( np.floor( trj_step_no / self.steps ) )
            i=0
            for ts in self.univ.trajectory[0:-1:ds]:
                pdb.update_progress( float(i)/self.steps )
                self.Atom_Select()
                self.DigestTraj()
                self.InitMarchingCubes()
                self.SmoothSurf()
                gc = self.GaussianCurvature()
                mc = self.MeanCurvature()
                if verbose:
                    print  "\n",ts,"\ntotal gaussian curvature:",gc
                    print  "\n",ts,"\ntotal mean curvature     :",mc
                total_gaussian_array.append(gc)
                total_mean_array.append(mc)
                time_array.append(ts.time)
                i+=1
            WriteXVG(x_array=time_array , y_array=total_gaussian_array , titile="Total Gaussian Curvature" , file_name="gaussian_curvature.xvg" , x_lable="Time(ps)" , y_lable=ANG2)
            WriteXVG(x_array=time_array , y_array=total_mean_array       , titile="Total Mean Curvature"      , file_name="mean_curvature.xvg"       , x_lable="Time(ps)" , y_lable=ANG1)
            
    def Extract2(self,frame=0,file_out_gauss="gauss.txt",file_out_mean="mean.txt",file_out_gro=False):
        """calculates the curvature on the given frame and writes the curvature value for every beads
        on the given file to be imported in VMD.
        """
        print "total number of frames in this trajectory is: ",len (self.univ.trajectory)
        self.univ.trajectory[frame]
        self.DigestGro()
        self.InitMarchingCubes()
        self.SmoothSurf()
        ###################
        scalar_range = [-0.002,0.001]
        i=0
        t=float(len(self.univ.atoms))
        self.GaussianCurvature(verbose=False)
        last_resid = -1
        file_string = ""
        print "Molecules gaussian curvature assignment ... "
        for atom in self.univ.atoms:
            pdb.update_progress(float(i)/t)
            i+=1
            if atom.resname in ["W","WF"]:
                file_string += "%.6f "%0.0
            else:
                if last_resid != atom.resid :
                    last_resid = atom.resid
                    molecule_curv = self.MoleculeCurvature(atom.resid)[0]
                    molecule_curv = LimitToRange(molecule_curv,scalar_range)
                else:
                    pass # molecule_curv has already been calculated for this resid.
                file_string += "%.6f "% molecule_curv
        new_file = file(file_out_gauss, 'w')
        new_file.write(file_string)
        new_file.close()        
        #SavePolyData(self.curvature_polydata , "gaussian_polydata.vtp")
        self.test_MoleculeCurvature()
        
        ###  repeat for mean curvature ###
        scalar_range = [-0.002,0.001]
        i=0
        t=float(len(self.univ.atoms))
        self.MeanCurvature()
        last_resid = -1
        file_string = ""
        print "\nMolecules mean curvature assignment ... "
        for atom in self.univ.atoms:
            pdb.update_progress(float(i)/t)
            i+=1
            if atom.resname in ["W","WF"]:
                file_string += "%.6f "%0.0
            else:
                if last_resid != atom.resid :
                    last_resid = atom.resid
                    value_id   = self.MoleculeCurvature(atom.resid)
                    molecule_curv = value_id [0]
                    molecule_curv = LimitToRange(molecule_curv,scalar_range)
                else:
                    pass # molecule_curv has already calculated for this resid.
                file_string += "%.6f " % molecule_curv
                
        new_file = file(file_out_mean, 'w')
        new_file.write(file_string)
        new_file.close()
        #SavePolyData(self.curvature_polydata , "mean_polydata.vtp")
        
        if file_out_gro:
            self.univ.atoms.write(file_out_gro)
            
    def GaussianCurvature(self,verbose=False):
        """ it also fixes self.triangles_centers_pointlocator
        """
        self.gauss_mean=True
        value_list = []
        self.curvaturesFilter_gauss = vtk.vtkCurvatures()
        self.curvaturesFilter_gauss.SetInputConnection( self.normal_generator_smooth.GetOutputPort() )
        self.curvaturesFilter_gauss.SetCurvatureTypeToGaussian()
        self.curvaturesFilter_gauss.Update()
        self.centers_points = vtk.vtkPoints()   # update self.centers_point to contain centers of the triangles.
        
        # Use a color series to create a transfer function
        scalar_range = [-0.002,0.001]
        #self.curvaturesFilter_gauss.GetOutput().GetScalarRange(scalar_range)
        if verbose:
            print "scalar_range:  ",scalar_range
            lut = vtk.vtkColorTransferFunction()
            lut.SetColorSpaceToHSV()
            num_colors = self.color_series.GetNumberOfColors()
            d_color=[0.0,0.0,0.0]
            for i in range(num_colors):
                color = self.color_series.GetColor(i)
                d_color[0] =  color[0] / 255.0
                d_color[1] =  color[1] / 255.0
                d_color[2] =  color[2] / 255.0
                t= scalar_range[0] + (scalar_range[1] - scalar_range[0]) * i / (num_colors - 1) 
                lut.AddRGBPoint(t, d_color[0], d_color[1], d_color[2]);
                # print t, "  :  ",d_color[0], d_color[1], d_color[2]
        ## make the value list:
        poly_data = self.curvaturesFilter_gauss.GetOutput()
        self.curvature_polydata=poly_data
        cell_number = poly_data.GetNumberOfCells()
        #poly_data.GetPointData.SetActiveScalars("Mean_Curvature")
        point_data = poly_data.GetPointData()
        data_array = point_data.GetScalars()
        self.total_GC = 0  # total Gaussian curvature
        triangle   = np.zeros([3,3],dtype="f") # variable to store triangles
        if verbose:
            print "first value: ",data_array.GetTuple1(0)
        t = float(cell_number)
        
        for i in range(cell_number):
            value=data_array.GetTuple1(i) 
            if verbose:
                pdb.update_progress(float(i)/t)
            cell        = poly_data.GetCell(i)
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                poly_data.GetPoints().GetPoint( point_id , triangle[j] )
            if  abs(value) < 1.0: # integrate only if ...
                self.total_GC = self.total_GC + ( CalcArea(triangle)  * value )
                value_list.append(value)
            # store the center of the triangles:
            self.centers_points.InsertNextPoint(TriangleCenter(triangle))
        
        centers_polydata = vtk.vtkPolyData()
        centers_polydata.SetPoints(self.centers_points)
        self.triangles_centers_pointlocator = vtk.vtkPointLocator()
        self.triangles_centers_pointlocator.SetDataSet(centers_polydata)
        self.triangles_centers_pointlocator.SetNumberOfPointsPerBucket(10)
        self.triangles_centers_pointlocator.BuildLocator()
        
        if verbose:
            print "\nmax value: ", max(value_list)
            print "min value: ",min(value_list)
            print "bounds of the mesh triangle centers of the curvature:",self.centers_points.GetBounds( )    
            print "total Gaussian curvature: ", self.total_GC
            plt.hist(value_list,bins=100)
            plt.suptitle('Gaussian curvature distribution')
            plt.show()
    
            ##Create a mapper and actor
            mapper = vtk.vtkPolyDataMapper()
            mapper.UseLookupTableScalarRangeOn()
            mapper.SetInputConnection(self.curvaturesFilter_gauss.GetOutputPort())
            mapper.SetLookupTable(lut)
            mapper.SetScalarRange(scalar_range)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
        
            ##Create a scalar bar actor    
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(mapper.GetLookupTable())
            scalarBar.SetNumberOfLabels(5)
            ##Create text actor
            text_actor = vtk.vtkTextActor()
            text_actor.SetInput ( "Gaussian Curvature" )
            text_actor.GetTextProperty().SetFontSize ( 20 )
            text_actor.GetTextProperty().SetColor ( 0.5, 0.5, 0.5 )
            text_actor.SetPosition(100, 16)
        
            self.gauss_renderer = vtk.vtkRenderer()
            self.gauss_renderer.SetBackground([1.,1.,1.])
            self.gauss_renderer.AddActor(actor)
            self.gauss_renderer.AddActor2D(scalarBar)
            self.gauss_renderer.AddActor2D ( text_actor )
        return self.total_GC
    
    def MeanCurvature(self,verbose=False):
        self.gauss_mean=False
        value_list = []
        self.curvaturesFilter_mean = vtk.vtkCurvatures()
        self.curvaturesFilter_mean.SetInputConnection( self.normal_generator_smooth.GetOutputPort() )
        self.curvaturesFilter_mean.SetCurvatureTypeToMean()
        self.curvaturesFilter_mean.Update();
        # Use a color series to create a transfer function
        scalar_range = [-0.002,0.001]
        #self.curvaturesFilter_gauss.GetOutput().GetScalarRange(scalar_range)
        if verbose:
            print "scalar_range:  ",scalar_range
            lut = vtk.vtkColorTransferFunction()
            lut.SetColorSpaceToHSV()
            num_colors = self.color_series.GetNumberOfColors()
            d_color=[0.0,0.0,0.0]
            for i in range(num_colors):
                color = self.color_series.GetColor(i)
                d_color[0] =  color[0] / 255.0
                d_color[1] =  color[1] / 255.0
                d_color[2] =  color[2] / 255.0
                t= scalar_range[0] + (scalar_range[1] - scalar_range[0]) * i / (num_colors - 1) 
                lut.AddRGBPoint(t, d_color[0], d_color[1], d_color[2]);
                # print t, "  :  ",d_color[0], d_color[1], d_color[2]
        ## make the value list:
        poly_data = self.curvaturesFilter_mean.GetOutput()
        self.curvature_polydata=poly_data
        cell_number = poly_data.GetNumberOfCells()
        #poly_data.GetPointData.SetActiveScalars("Mean_Curvature")
        point_data = poly_data.GetPointData()
        data_array = point_data.GetScalars()
        self.total_MC = 0  # total Gaussian curvature
        triangle   = np.zeros([3,3],dtype="f") # variable to store triangles
        if verbose:
            print "first value: ",data_array.GetTuple1(0)
        t = float(cell_number)
        
        for i in range(cell_number):
            value=data_array.GetTuple1(i)
            if verbose:
                pdb.update_progress(float(i)/t)
            cell        = poly_data.GetCell(i)
            for j in range(3):
                point_id    =    cell.GetPointId(j)
                poly_data.GetPoints().GetPoint( point_id , triangle[j] )
            if  abs(value) < 1.0: # integrate only if ...
                self.total_MC = self.total_MC + ( CalcArea(triangle)  * value )
                value_list.append(float(value))
            
        if verbose:
            print "\nmax value: ", max(value_list)
            print "min value    : ",min(value_list)
            print "total Gaussian curvature: ", self.total_MC
            plt.hist(value_list,bins=100)
            plt.suptitle('Mean curvature distribution')
            plt.show()
    
            ##Create a mapper and actor
            mapper = vtk.vtkPolyDataMapper()
            mapper.UseLookupTableScalarRangeOn()
            mapper.SetInputConnection(self.curvaturesFilter_mean.GetOutputPort())
            mapper.SetLookupTable(lut)
            mapper.SetScalarRange(scalar_range)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
        
            ##Create a scalar bar actor    
            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(mapper.GetLookupTable())
            scalarBar.SetNumberOfLabels(5)
            ##Create text actor
            text_actor = vtk.vtkTextActor()
            text_actor.SetInput ( "Mean Curvature" )
            text_actor.GetTextProperty().SetFontSize ( 20 )
            text_actor.GetTextProperty().SetColor ( 0.5, 0.5, 0.5 )
            text_actor.SetPosition(100, 16)
        
            self.meancurve_renderer = vtk.vtkRenderer()
            self.meancurve_renderer.SetBackground([1.,1.,1.])
            self.meancurve_renderer.AddActor(actor)
            self.meancurve_renderer.AddActor2D(scalarBar)
            self.meancurve_renderer.AddActor2D ( text_actor )

        return self.total_MC

    def MovementCurvature(self,file_name="test.vtp"):
        """forming the movement, curvature dictionary of different molecules (resid). 
        self.curvature_movement_dict = dict()  dictionary {resid:[curvature,movement]} [not used]
        velocities are used as the movement measure.
        """
        ## make the vtkpolydata, .. of the movement profile of the curvature.
        movement_polydata = vtk.vtkPolyData()
        movement_polydata.DeepCopy(self.curvature_polydata)
        no_of_cells = movement_polydata.GetPointData().GetScalars().GetSize()
        
        ##
        resid_set        = set()
        whole_atoms = self.univ.atoms
        for atom in whole_atoms:
            resid_set.add(atom.resid)
        #self.curvature_movement_dict = dict()
        t=len(resid_set)
        i=0
        triangleid_movements_dict = dict() # {triangle_id : list of movements of molecules in the triangle. }
        print "\nfetching movement and curvature for each molecule: "
        for id in resid_set:
            i+=1
            pdb.update_progress(float(i)/t)
            #movement = self.VelocityAlongSurf(id)
            #movement = self.MovementAlongSurfTraj(id)
            movement = self. BeadsMovementTraj(id)
            aux            = self.MoleculeCurvature(id)
            curvature   = aux[0]
            triangle_id = aux[1]
            #self.curvature_movement_dict[id] = [curvature,movement]
            if triangle_id in triangleid_movements_dict:
                triangleid_movements_dict[triangle_id] = np.append(triangleid_movements_dict[triangle_id] ,movement)
            else:
                #initialize np.array
                triangleid_movements_dict[triangle_id] = np.array([movement])
                
        for id in range(no_of_cells):
            if id in triangleid_movements_dict:
                movement_average = np.average( triangleid_movements_dict[id] )
            else:
                movement_average = 0.0
#             if movement_average>1000:
#                 print "\nmovement_average: ",movement_average,"triangle_id : ",id
            movement_polydata.GetPointData().GetScalars().SetTuple1(id,movement_average)
            
        ### save on file ##
        if file_name:
            SavePolyData(movement_polydata,file_name)
        ##### visualizations ###
        ShowPolyData(movement_polydata)
        
    def MovementVisualize(self,m=3):
        """Visualizeing the movement of the beads in the molecules in the last m frames.
        """
        resid_set        = set()
        whole_atoms = self.univ.atoms
        glyph3d = vtk.vtkGlyph3D()
        id_mov  = dict()
        
        for atom in whole_atoms:
            resid_set.add(atom.resid)
        
        scalar_range = [1e6,-1e6]
        t=len(resid_set)
        for id in resid_set:
            i+=1
            pdb.update_progress(float(i)/t)
            movement = self.BeadsMovementTraj(id,m)
            id_mov[id] = movement
            scalar_range[0] = min ( movement , scalar_range[0] )
            scalar_range[1] = max ( movement , scalar_range[1] )
            
        lut = vtk.vtkColorTransferFunction()
        lut.SetColorSpaceToHSV()
        num_colors = self.color_series.GetNumberOfColors()
        d_color=[0.0,0.0,0.0]
        for i in range(num_colors):
            color = self.color_series.GetColor(i)
            d_color[0] =  color[0] / 255.0
            d_color[1] =  color[1] / 255.0
            d_color[2] =  color[2] / 255.0
            t= scalar_range[0] + (scalar_range[1] - scalar_range[0]) * i / (num_colors - 1) 
            lut.AddRGBPoint(t, d_color[0], d_color[1], d_color[2]);
            
        glyph3d = vtk.vtkGlyph3D()
        colors = vtk.vtk
        pass
    
    def MoleculeCurvature(self,in_residue):
        """last calculated curvature which also updates triangles_centers_pointlocator.
            * it takes into account the last curvature extraction that has been run (mean or gaussian).
        """
        point_data = self.curvature_polydata.GetPointData() # the last updated self.curvature_polydata
        data_array = point_data.GetScalars()
        pos = np.array([0.0,0.0,0.0])
        i=0
        molecule_atoms = self.univ.select_atoms("resid "+str(in_residue)+" and name H1 C1" ) # C1 C2 C3 are the head beads for glycoceramides. 
        for atom in molecule_atoms:  # average position if necessary
            pos +=atom.position
            i+=1
        pos /=i
        point_id = self.triangles_centers_pointlocator.FindClosestPoint(pos)
        point_coor=self.triangles_centers_pointlocator.GetDataSet().GetPoints().GetPoint(point_id)
        #distance = vtk.vtkMath.Distance2BetweenPoints(point_coor,pos)
        value=data_array.GetTuple1(point_id)
        
        ## for test ##
        index1 = self.test_points.InsertNextPoint(pos)
        index2 = self.test_points.InsertNextPoint(point_coor)
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, index1)
        line.GetPointIds().SetId(1, index2)
        self.test_lines.InsertNextCell( line )
        ####
        
        return [value,point_id]

    def test_MoleculeCurvature(self):
        """visualize the lines between molecules and their corresponding triangle for testing purposes
        """

        linesPolyData = vtk.vtkPolyData()
        linesPolyData.SetPoints(self.test_points)
        linesPolyData.SetLines(self.test_lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(linesPolyData)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        ##
        # bound curvature polydata:
        scalar_range = [-0.002,0.001]
        print  "\nvtkDataArray max id: ",self.curvature_polydata.GetPointData().GetScalars().GetMaxId(),"number of cell: ", self.curvature_polydata.GetNumberOfCells (),"number of points: ",self.curvature_polydata.GetNumberOfPoints()
        for index in range( self.curvature_polydata.GetNumberOfCells() ) :
            self.curvature_polydata.GetPointData().GetScalars().SetTuple1( 0,LimitToRange(self.curvature_polydata.GetPointData().GetScalars().GetTuple1(index) ,scalar_range) )
        ##
        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(self.curvature_polydata)
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        ##
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.AddActor(actor2)
        window  = vtk.vtkRenderWindow()
        window.AddRenderer(renderer)
        interactor= vtk.vtkRenderWindowInteractor()
        interactor.SetInteractorStyle(molar.pdb_viewer.PdbInteractorStyle())
        interactor.SetRenderWindow(window)
        interactor.Start()
        window.Render()
        pass

    def VelocityAlongSurf(self,resid_):
        """given the resid, it returns the velocity of the molecule with the resid=resid_ along (parallel to/projected on)the surface which should be already calculated.
        this calculation is done using just the position and velocity on the input gro file (self.univ_gro). 
        two other possible options for extracting the velocities: 1. save velocities on the trajectory file and use them here. 2. use positions in
        two consecutive frames to calculate the velocity. 
        """
        atoms=self.univ_gro.select_atoms( "resid "+str(resid_) )
        i=0
        v     = np.array([0.0,0.0,0.0])
        pos = np.array([0.0,0.0,0.0])
        for a in atoms:
            v    +=a.velocity
            pos+=a.position
            i+=1
        v    /=i  # v = velocity vector of the molecule averaged over the atoms
        pos/=i  # pos = position of the molecule averaged over the atoms
        
        point_id = self.triangles_centers_pointlocator.FindClosestPoint(pos)
        normal =  self.normal_generator_smooth.GetOutput().GetPointData().GetNormals().GetTuple(point_id) #self.normal_generator_smooth.GetCellData().GetNormals()
        normal = np.array(normal)
        normal = normal / np.linalg.norm(normal)
        out       = v  -  (np.dot(v,normal)) * normal
        return np.linalg.norm(out)
    
    def MovementAlongSurfTraj(self,resid_,m=5):
        """ movement on the last m frames of trajectory projected on the surface (using self.univ instead of self.univ_gro). similar to VelocityAlongSurf
        """
        atoms  =  self.univ.select_atoms( "resid "+str(resid_) )
        switch  = True
        out       = 0
        for ts in self.univ.trajectory[-m:]:
            i=0
            pos = np.array([0.0,0.0,0.0])
            for a in atoms:
                pos+=a.position
                i+=1
            pos/=i  # pos = position of the molecule averaged over the atoms
            if switch:
                previous_pos = pos
                switch = False
                continue
            dp = pos - previous_pos
            previous_pos = pos
            point_id = self.triangles_centers_pointlocator.FindClosestPoint(pos)
            normal =  self.normal_generator_smooth.GetOutput().GetPointData().GetNormals().GetTuple(point_id) #self.normal_generator_smooth.GetCellData().GetNormals()
            normal = np.array(normal)
            normal = normal / np.linalg.norm(normal)
            out       = out + np.linalg.norm(dp  -  (np.dot(dp,normal)) * normal)  # adding the magnitude of dp to out
        return out
        
    def BeadsMovementTraj(self,resid_,m=5):
        """ movement on the last m frames of trajectory (using self.univ instead of self.univ_gro). similar to MovementAlongSurfTraj but not projected on surface and not averaged position of molecule.
        """
        atoms  =  self.univ.select_atoms( "resid "+str(resid_) )
        switch  = True
        #atoms_position =[np.zeros(3)] * len(atoms) # list of numpy.arraysopen 
        total_disposition = 0.0
        for ts in self.univ.trajectory[-m:]:
            if switch:
                atoms_previous_pos = [ np.copy(a.position) for a in atoms ] # making hard copies of elements of atoms_position
                switch = False
                continue
            i=0
            for a in atoms:
                total_disposition = total_disposition + np.linalg.norm( a.position - atoms_previous_pos[i] )
                atoms_previous_pos[i] = np.copy(a.position)
                i+=1
        return total_disposition
        
    def Show(self):
        """
        """
        window    = vtk.vtkRenderWindow()
        window.SetSize(1800, 1200)
        interactor= vtk.vtkRenderWindowInteractor()
        interactor.SetInteractorStyle(molar.pdb_viewer.PdbInteractorStyle())
        interactor.SetRenderWindow(window)
        top_left_viewport       = [0.0    ,0.5  ,0.33   ,1.0]
        top_center_viewport  = [0.33  ,0.5  ,0.66   ,1.0]
        top_right_viewport     = [0.66  ,0.5  ,1.0     ,1.0]
        down_left_viewport       = [0.0    ,0.0  ,0.33   ,0.5]
        down_center_viewport  = [0.33  ,0.0  ,0.66   ,0.5]
        down_right_viewport     = [0.66  ,0.0  ,1.0     ,0.5]
        # setup three renderers
        if self.pdb_water_head:
            renderer_0 = self.pdb_water_head.RetrunRenderer( background=[.0, .0, .0] )
            renderer_0.SetViewport(top_left_viewport)
        else:
            ## dummy renderer
            renderer_0 = vtk.vtkRenderer()
            renderer_0.SetViewport(down_right_viewport)
        renderer_1 = MakeRenderer(self.marching_cubes, color=[0.0,0.95,0.0], opacity=1.0, backgorund=[.0, .0, .0] ,viewport=top_center_viewport)
        renderer_2 = MakeRenderer(self.normal_generator_smooth, color=[0.0,0.95,0.0], opacity=1.0, backgorund=[.0, .0, .0] ,viewport=top_right_viewport)
        self.gauss_renderer.SetViewport(down_left_viewport)
        self.meancurve_renderer.SetViewport(down_center_viewport)
        print "renderers have been made."
        ## dummy renderer
        dummy_renderer = vtk.vtkRenderer()
        dummy_renderer.SetViewport(down_right_viewport)
        ##
        window.AddRenderer(renderer_0)
        window.AddRenderer(renderer_1)
        window.AddRenderer(renderer_2)
        window.AddRenderer(self.gauss_renderer)
        window.AddRenderer(self.meancurve_renderer)
        window.AddRenderer(dummy_renderer)
        window.Render();
        interactor.Start();
        
##############################
##############################
def LimitToRange(a,scalar_range=[0.0,1.0]):
    if a<scalar_range[0] :
        return scalar_range[0] 
    elif a>scalar_range[1] :
        return scalar_range[1] 
    else:
        return a

def ShowVtkAlgorithm(in_vtkAlgorithm):
    mapper = vtk.vtkDataSetMapper()  #vtkPolyDataMapper() does not work with vtkDelaunay3D
    mapper.SetInputConnection(in_vtkAlgorithm.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground([0.5,0.5,0.5])
    renderer.AddActor(actor)
    window    = vtk.vtkRenderWindow()
    interactor= vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(molar.pdb_viewer.PdbInteractorStyle())
    interactor.SetRenderWindow(window)
    window.SetSize(1800, 1200)
    window.AddRenderer(renderer)
    window.Render()
    interactor.Start()    

def ShowPolyData(in_polydata,caption=""):
    """Visualize the in_polydata surface which has scalar values which will be colour-coded.
    """
    if True: 
        color_series = vtk.vtkColorSeries()
        color_series.SetColorScheme( color_series.BREWER_DIVERGING_SPECTRAL_10 )   #set the color theme
        lut = vtk.vtkColorTransferFunction()
        lut.SetColorSpaceToHSV()
        num_colors = color_series.GetNumberOfColors()
        d_color         =  [0.0,0.0,0.0]
        scalar_range =  [0.0,0.0]  #[-0.002,0.001]
#         in_polydata.GetPointData().GetScalars().GetRange(scalar_range)
#         for i in range(num_colors):
#                 color = color_series.GetColor(i)
#                 d_color[0] =  color[0] / 255.0
#                 d_color[1] =  color[1] / 255.0
#                 d_color[2] =  color[2] / 255.0
#                 t= scalar_range[0] + (scalar_range[1] - scalar_range[0]) * i / (num_colors - 1) 
#                 lut.AddRGBPoint(t, d_color[0], d_color[1], d_color[2]);
                
        ##Create a mapper and actor
        mapper = vtk.vtkPolyDataMapper()
#         mapper.UseLookupTableScalarRangeOn()
        mapper.SetInputData(in_polydata)    #SetInputConnection() is for algorithm not polydata
#         mapper.SetLookupTable(lut)
#         mapper.SetScalarRange(scalar_range)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
         
        ##Create a scalar bar actor    
        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetLookupTable(mapper.GetLookupTable())
        scalarBar.SetNumberOfLabels(5)
        
        ##Create text actor
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput ( caption )
        text_actor.GetTextProperty().SetFontSize ( 20 )
        text_actor.GetTextProperty().SetColor ( 0.5, 0.5, 0.5 )
        text_actor.SetPosition(100, 16)

        renderer = vtk.vtkRenderer()
        renderer.SetBackground([0.5,0.5,0.5])
        renderer.AddActor(actor)
        renderer.AddActor2D(scalarBar)
        renderer.AddActor2D ( text_actor )
        
        ####
        window    = vtk.vtkRenderWindow()
        interactor= vtk.vtkRenderWindowInteractor()
        interactor.SetInteractorStyle(molar.pdb_viewer.PdbInteractorStyle())
        interactor.SetRenderWindow(window)
        window.SetSize(1800, 1200)
        window.AddRenderer(renderer)
        window.Render()
        interactor.Start()    
    
def SavePolyData(vtkpolydata_in,vtp_file_name_in):  #file has vtp extension by convention.#
    ### save on file ##
    if vtp_file_name_in:
        xmlwrilter = vtk.vtkXMLPolyDataWriter()
        xmlwrilter.SetFileName(vtp_file_name_in)
        xmlwrilter.SetInputData(vtkpolydata_in)
        xmlwrilter.Write()
        print "\nmovement polydata is written to the file: ", vtp_file_name_in
    ## ## 

def MakeRenderer(vtkAlgorithm_obj , color=[0.0,0.95,0.0] , opacity = 1.0 , backgorund=[0.1,0.2,0.3] ,viewport=[0.0,0.0,1.0,1.0] ):
    """make a renderer from the vtkAlgorithm
    """
    actor = vtk.vtkActor()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection( vtkAlgorithm_obj.GetOutputPort() )
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetOpacity(opacity)
    renderer  = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetViewport(viewport)
    renderer.SetBackground(backgorund)
    return renderer
    
def CalcArea(in_array):
    """returning normal vector of a triangle with sides of v1 and v2
        that has the length equal to the area of the triangle. 
    """
    coor = np.array(in_array)
    v1 = coor[1] - coor[0]
    v2 = coor[2] - coor[0]
    c  = 0.5*np.cross(v1,v2)
    return np.linalg.norm(c)

def TriangleCenter(triangle):
    return (triangle[0] + triangle[1] + triangle[2]) / 3.0

def WriteXVG(x_array=[] , y_array=[] , titile="" , file_name="out_xvg.xvg" , x_lable="" , y_lable=""):
    # ( A\\v{0.65}\h{-0.5}\z{0.6}o  -2\\v{}\h{}\z{} )
    xvg_header = "@    title \""+titile+"\"\n\
@    xaxis  label \""+x_lable+"\"\n\
@    yaxis  label \""+y_lable+"\"\n\
@    TYPE xy\n"
    if len(x_array) != len(y_array):
        sys.stderr.write("Length of x_array and y_array should match!\n")
    out_str = xvg_header
    for x,y in zip(x_array,y_array):
        out_str = out_str + ( "%15.3f     %15.3f\n"%(x, y ) )
    new_file = file(file_name, 'w')
    new_file.write(out_str)
    new_file.close()
    
def WriteHistXVG(x_array=[]  ,bins=10, titile="" , file_name="out_hist.xvg" , x_lable="" ):
    """Draw histogram of input.
    """
    xvg_header = "@    title \""+titile+"\"\n\
@    xaxis  label \""+x_lable+"\"\n\
@    s0 line type 0\n\
@    s0 symbol fill pattern 4\n\
@    TYPE bar\n"
    hist = np.histogram(x_array, bins, range=None, normed=False, weights=None, density=None)
    mids = []

    for i in range( len(hist[1]) - 1 ):
        mids.append(0.5*(hist[1][i]+hist[1][i+1]) )
    out_str = xvg_header
    for x,y in zip(mids,hist[0]):
            out_str = out_str + ( "%15.4f     %15d\n"%(x, y ) )
    new_file = file(file_name, 'w')
    new_file.write(out_str)
    new_file.close()
