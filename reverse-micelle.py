import membrane.membrane as mem

m=mem.Membrane(unit_cell_size_ = 200)
m.form = "reverse-micelle"
m.micelle_radius = 80
m.MakeMartini6(1.0)
m.AddWater()
m.Show()
m.WriteOnFileGMX("reverse_micelle_s20_8.pdb")

