# make cubic membrane, visualize and write on file 
import membrane.membrane as mem

m=mem.Membrane(unit_cell_size_=200)
# choose the form between "G", "P", "D", "micelle", "reveerse-micelle"
m.form = "G"

# start making the system with 100% glycoderamides
m.MakeMartini3(glyco=1.0)
m.Show()
m.WriteOnFileGMX("cubic_g_20_gly1.0.pdb")
#m.PrintInfo()
print "number of molecules: ",len(m.molecules)
print "sample_size=",m.sampling_size
