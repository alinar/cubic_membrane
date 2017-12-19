# This will take the coordinate and trajectory. On the specified frame it will extract the surface and calculate the curvature. it will write the curvature value associated with each bead/atom on a txt file. (this is useful when visualizing the curvature values with VMD)

import membrane.extract_surf as ex

ex_obj=ex.ExtractSurf("input_coordinate_file.gro","input_trajectory_file.xtc")

ex_obj.Extract2(frame=500,curve_file_out="gauss.txt",gaussian_mean=True,file_out_gro="frame_gauss.gro")


