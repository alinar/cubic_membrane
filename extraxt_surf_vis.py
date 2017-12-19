# This will extract the lipid-water surrface from the input.pdb file and visualize it.

import membrane.extract_surf as ex
ex_obj=ex.ExtractSurf("input.pdb")
ex_obj.Extract(verbose=True)
print "Visualizations..."
ex_obj.Show()
