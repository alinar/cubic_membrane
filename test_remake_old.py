import membrane.membrane1 as mem
m = mem.Membrane(1,250)
m.Make()
#print m.Bounds()
#m.Show()
# m.VisVertices()
m.BringToCenter()
m.PrintInfo()
m.WriteOnFile("s25.pdb")
print "FILENAME" , " DONE."
