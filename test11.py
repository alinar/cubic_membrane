import membrane.extract_surf as ex
#ex_obj=ex.ExtractSurf("/Users/alinar/tmp/confout.gro","/Users/alinar/tmp/traj_comp.xtc", 2)
ex_obj=ex.ExtractSurf("/Users/alinar/Downloads/cubic/confout.gro","/Users/alinar/Downloads/cubic/traj_comp.xtc", 2)
ex_obj.Extract(verbose=False)
ex_obj.MovementCurvature()

