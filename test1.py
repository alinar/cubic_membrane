import membrane.membrane1 as mem
m = mem.Membrane(2,220)
m.form = "G"
if "cereb" == "cereb": 
 m.cerebroside = True
 print "cereb is activated."
m.Make()
#print m.Bounds()
#m.Show()
# m.VisVertices()
m.BringToCenter()
m.PrintInfo()
m.WriteOnFile("cereb_cubic_G_220_-0d.pdb")
print "cereb_cubic_G_220_-0d.pdb" , " DONE."
