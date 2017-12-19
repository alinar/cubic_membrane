#extract surface and calculate curvature and write vtp file.

import membrane.extract_surf as ex
ex_obj=ex.ExtractSurf("input.gro","input_traj.xtc")

ex_obj.DigestGro()
ex_obj.InitMarchingCubes()
ex_obj.SmoothSurf()
ex_obj.GaussianCurvature(verbose=False)
pd = ex_obj.curvature_polydata
ex.SavePolyData(pd,"test_pd.vtp")
