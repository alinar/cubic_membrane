# Extract the water-lipid surface from the input.gro file, do the curvature calculateions (Gaussian and mean)
# and wirte the resaults on two separate files.

import membrane.extract_surf as ex

ex_obj=ex.ExtractSurf("input.gro")
ex_obj.Extract(verbose=True,gaussian_polydata_file_out="gauss.vtp",mean_polydata_file_out="mean.vtp")
# visualize
ex_obj.Show()
