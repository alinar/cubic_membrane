import membrane.membrane as mem

m=mem.Membrane(unit_cell_size_ = 200)
m.form = "micelle"
m.MakeMartini5(1.0)
m.AddWater()
#m.Show()
m.WriteOnFileGMX("micelle_s20.pdb")
