import membrane.make_box as box

#box.Make_test()
#box=box.Box(x=753,y=753,z=300)
box=box.Box(x=300,y=300,z=100)
box.MaskPdb(pdb_in_file="test.pdb",out_pdb_filename="test_boxed.pdb",specimen_rotation_angles=[45,45,45])

