# This will mask input.pdb file inside a box with a default size and the given specimen_rotation_angles and writes the new file on masked_pdb.pdb.

import membrane.make_box as box

box=box.Box()
box.WriteTclFile()
box.MaskPdb(pdb_in_file="input.pdb",out_pdb_filename="masked_pdb.pdb",specimen_rotation_angles=[50,50,50])

