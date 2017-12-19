import molar.pdb as p
import numpy as np
import random
m=p.Pdb("cernp24.pdb")

c1=m.GetAtomByName("C2F")
c2=m.GetAtomByName("C6F")

c3=m.GetAtomByName("C4S")
c4=m.GetAtomByName("C8S")

a1=m.GetAtomByName("C2S")
a2=m.GetAtomByName("C3S")

a3=m.GetAtomByName("NF")
a4=m.GetAtomByName("C1F")

a5=m.GetAtomByName("C1F")
a6=m.GetAtomByName("C2F")

###############

def Change(amp = 10):
    alpha    =  amp * 2 * (random.random() - 0.5 )
    beta      =  amp * 2 * (random.random() - 0.5 )
    gamma =  amp * 2 * (random.random() - 0.5 )

    m.SelectBond(a2,a1)
    m.RotateBond(alpha)

    m.SelectBond(a3,a4)
    m.RotateBond(beta)

    m.SelectBond(a5,a6)
    m.RotateBond(gamma)
    
    return [alpha,beta,gamma]

#################

def UndoChange(rotation_list):
    alpha      = rotation_list[0]
    beta        = rotation_list[1]
    gamma   = rotation_list[2]
    m.SelectBond(a2,a1)
    m.RotateBond(-1*alpha)

    m.SelectBond(a3,a4)
    m.RotateBond(-1*beta)

    m.SelectBond(a5,a6)
    m.RotateBond(-1*gamma)
 
################

def Cost():
  #print np.dot((c2.pos-c1.pos),(c4.pos-c3.pos))
  return np.dot((c2.pos-c1.pos),(c4.pos-c3.pos))
  ## return np.degrees(np.arcsin( np.linalg.norm( np.cross((c2.pos-c1.pos),(c4.pos-c3.pos)) ) / (np.linalg.norm(c2.pos-c1.pos) * np.linalg.norm(c4.pos-c3.pos) ) ) )

################
s_alpha     = 0
s_beta      = 0
s_gamma = 0
cost = Cost()
print "angle: " ,  np.degrees( np.arccos(cost /  ( np.linalg.norm(c2.pos-c1.pos) * np.linalg.norm(c4.pos-c3.pos ) ) ) ), "cost: ",cost
max =  np.linalg.norm(c2.pos-c1.pos) * np.linalg.norm(c4.pos-c3.pos )
while cost < 0.99*max: 
    rotation_list_ = Change(amp = 90)
    cost_aux = Cost() 
    if cost_aux > cost :
        cost = cost_aux
        s_alpha    += rotation_list_[0]
        s_beta      += rotation_list_[1]
        s_gamma += rotation_list_[2]
    else:
        UndoChange(rotation_list_) 
    print "angle: " , np.degrees( np.arccos(cost /  ( np.linalg.norm(c2.pos-c1.pos) * np.linalg.norm(c4.pos-c3.pos ) ) ) ) , "cost: ",cost, " -->> ",max


########################
m.Show("sphere")
m.WriteOnFile(file_name_str="test.pdb",make_TER=False,include_CONECT=True)
print "alpha: ", s_alpha, "beta:  ", s_beta , "gamma: ", s_gamma
