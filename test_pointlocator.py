import vtk
pl = vtk.vtkPointLocator()
pl.Initialize()
pl.SetDataSet(vtk.vtkPolyData())
#pl.InitPointInsertion(vtk.vtkPoints(), [-100,100,-100,100,-100,100] )

#pl.Update()
pl.GetDataSet().SetPoints(vtk.vtkPoints())
print "h1"
pl.InsertNextPoint([0.,0.,0.])
pl.InsertNextPoint([0.,0.,1.])
pl.InsertNextPoint([0.,1.,0.])
pl.InsertNextPoint([1.,0.,0.])
#pl.Update()
print "h2"
#pl.SetNumberOfPointsPerBucket(2)
#pl.BuildLocator()
print pl.FindClosestPoint([0.1,0,0])
