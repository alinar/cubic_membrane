import membrane.membrane as mem
import time
m = mem.Membrane(1,220)
m.form = "G"
m.MakeMartini()
#m.MakeMartini2()
m.WriteOnFile("cereb_cubic_G_220_-0d.pdb")
