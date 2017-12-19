import membrane.membrane as mem
import molar.pdb as pdb

m=mem.Membrane()
m.ReadFileGMX("vesicle_s15_5_no-water.pdb")
m.form = "vesicle"
m.side_length = 150


m.AddWater()

m.WriteOnFileGMX("vesicle_s15_5.pdb",make_TER=True)


