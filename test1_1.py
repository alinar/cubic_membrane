import membrane.membrane as mem
m = mem.Membrane(2,300)
m.surface_density  = 0.053
m.form = "G"
m.cerebroside = True
m.Make()
print m.Bounds()
#m.Show()
# m.VisVertices()
m.BringToCenter()
m.PrintInfo()
m.WriteOnFile("test.pdb")
