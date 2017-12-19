import molar.pdb as pdb
import vtk

cer = pdb.Pdb()
glu = pdb.Pdb()

cer.ReadFile("/Users/alinar/Dropbox/The Project/MicroBio/molecule_pool/lip1.pdb")
glu.ReadFile("/Users/alinar/Dropbox/The Project/MicroBio/molecule_pool/glucose.pdb")

c5 =glu.GetAtomsByName("C5")[0]
o6 =glu.GetAtomsByName("O6")[0]
c6 =cer.GetAtomsByName("C6")[0]
c7 =cer.GetAtomsByName("C7")[0]
o1 =cer.GetAtomsByName("O1")[0]
h49=cer.GetAtomsByName("H49")[0]

v1=o6.pos - c5.pos
v2=o1.pos - h49.pos 
v_pos = c7.pos - c6.pos

trans = pdb.RotateToParallel(v2,v1)
glu.ApplyTransform(trans)
trans.Identity()
trans.Translate(o1.pos - o6.pos)
#trans.Translate(1.0 * v_pos) #fine_tune
glu.RemoveAtom(o6)
glu.ApplyTransform(trans)
cer.RemoveAtom(h49)
result=pdb.MergePdb([cer,glu])

#result.Show()
result.WriteOnFile("/Users/alinar/Dropbox/The Project/MicroBio/molecule_pool/cerebroside.pdb")
