#! /usr/bin/env python

import vtk
import sys
import argparse
###########################

class PdbInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self,parent=None):
        #self.AddObserver("LeftButtonPressEvent",self.middleButtonPressEvent)
        #self.AddObserver("LeftButtonReleaseEvent",self.middleButtonReleaseEvent)
        #self.AddObserver("KeyPressEvent",self.MyOnKeyPress)
        pass
   
    def OnKeyPress(self,a,b):
        print "key pressed"
        super(PdbInteractorStyle,self).OnKeyPress(self,a, b)

    def middleButtonPressEvent(self,obj,event):
        print "Middle Button pressed"
        self.OnMiddleButtonDown()
        return

    def middleButtonReleaseEvent(self,obj,event):
        print "Middle Button released"
        self.OnMiddleButtonUp()
        return

###########################
###########################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This is a tool to 3D-visualize vtp files (vtkPolyData) with scalar values which will be colour-coded.')
    parser.add_argument('-i', '--input',type=str, help = 'Input vtp file name',required = True)
    parser.add_argument('-o', '--output',help = 'Output obj file name',type=str)
    parser.add_argument('-s','--show', dest='show', action='store_true',help = 'Visualize the input')
    parser.add_argument('-r','--range',type=float , help = 'Set the scalar range for visualization to [-RANGE,+RANGE].')
    parser.add_argument('-c', '--color',type=str, help = 'The name of the color scheme.\
                                        The default value is BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_3 .\
                                        For more information visit:\
                                        http://www.vtk.org/doc/nightly/html/classvtkColorSeries.html',\
                                        default='BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_3')
    args = parser.parse_args()
    

vtp_file = args.input
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(vtp_file)
reader.Update()
polydata = reader.GetOutput()

subdivisionFilter = vtk.vtkLoopSubdivisionFilter()
subdivisionFilter.SetNumberOfSubdivisions(2) #set the order
subdivisionFilter.SetInputData(polydata)
subdivisionFilter.Update()

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(subdivisionFilter.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

##### visualizations ###
## color lookup table ##
lut = vtk.vtkColorTransferFunction()
lut.SetColorSpaceToHSV()
color_series = vtk.vtkColorSeries()
#color_series.SetColorScheme(color_series.BREWER_SEQUENTIAL_YELLOW_ORANGE_BROWN_3)
color_series.SetColorScheme( getattr(color_series,args.color) )
print "\ncolor scheme name:  ",color_series.GetColorSchemeName()
num_colors = color_series.GetNumberOfColors()
d_color=[0.0,0.0,0.0]
if args.range:
    scalar_range = [-args.range,args.range]
else:
    scalar_range = [0.0,0.0]
    polydata.GetPointData().GetScalars().GetRange(scalar_range)
print "\ninput scalar range:",scalar_range
for i in range(num_colors):
    color = color_series.GetColor(i)
    d_color[0] =  color[0] / 255.0
    d_color[1] =  color[1] / 255.0
    d_color[2] =  color[2] / 255.0
    t= scalar_range[0] + (scalar_range[1] - scalar_range[0]) * i / (num_colors - 1) 
    lut.AddRGBPoint(t, d_color[0], d_color[1], d_color[2]);
##
mapper.SetLookupTable(lut)

scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(mapper.GetLookupTable())
scalarBar.SetNumberOfLabels(5)
  
renderer = vtk.vtkRenderer()
window  = vtk.vtkRenderWindow()
window.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetInteractorStyle(PdbInteractorStyle())
interactor.SetRenderWindow(window)

renderer.SetBackground(1.0,1.0,1.0)
renderer.AddActor(actor)
renderer.AddActor2D(scalarBar)
if args.show:
    window.Render()
    interactor.Start()
    
if args.output:
    obj = vtk.vtkOBJExporter()
    obj.SetInput(window)
    file_name = args.output
    if file_name[-4:] == ".obj":
        file_prefix = file_name[:-4]
    obj.SetFilePrefix(file_prefix)
    obj.Write()



#############################
#############################

    
