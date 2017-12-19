import membrane.membrane1 as mem
m = mem.Membrane(2,CELL_SIZE)
m.form = "MODE"
if "MOLEC" == "cereb": 
 m.cerebroside = True
 print "cereb is activated."
m.Make()
#print m.Bounds()
#m.Show()
# m.VisVertices()
m.BringToCenter()
m.PrintInfo()
m.WriteOnFile("FILENAME")
print "FILENAME" , " DONE."
