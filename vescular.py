# Making vescular systems

import membrane.membrane as mem

m=mem.Membrane(unit_cell_size_ = 150)
m.step_num = 16
m.micelle_radius = 50
m.MakeMartini8(1.0)
m.AddWater()
m.Show()
m.WriteOnFileGMX("vesicle_s15_5.pdb",make_TER=True)

