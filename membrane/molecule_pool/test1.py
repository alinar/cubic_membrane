import molar.pdb as p
import numpy as np
import random
m=p.Pdb("cereos.pdb")

m.WriteOnFile(file_name_str="test.pdb",make_TER=False,include_CONECT=True)


