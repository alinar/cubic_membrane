#!/usr/bin/python
# By Christian
import sys
import operator
import numpy as np
import MDAnalysis
from math import sqrt, pow
import argparse

#Data-object with trajectory
universe = MDAnalysis.Universe('only_h1.gro', 'only_h1.trr')


parser=argparse.ArgumentParser(description='Insert molecules into ceramide structure.')
parser.add_argument('--calculate', action="store_true", help="Calculate velocities")

args = parser.parse_args()

def PLANE_FROM_POINTS(pointvector):
    #This function is taken from http://www.ilikebigbits.com/blog/2015/3/2/plane-from-points

    #Normal vector to the plane will have components a,b,c
    a=0.0
    b=0.0
    c=0.0

    #Calculate centroid of point cloud
    number_of_points = len(pointvector)
    centroid = np.array([0.0,0.0,0.0])
    for point in pointvector:
        centroid += point
    centroid /= number_of_points


    #Calculate full 3x3 covariance matrix, excluding symmetries
    xx=0.0
    yy=0.0
    zz=0.0
    xy=0.0
    yz=0.0
    xz=0.0

    for point in pointvector:
        r=point-centroid
        xx += r[0] * r[0]
        xy += r[0] * r[1]
        xz += r[0] * r[2]
        yy += r[1] * r[1]
        yz += r[1] * r[2]
        zz += r[2] * r[2]

    det_x = yy*zz - yz*yz
    det_y = xx*zz - xz*xz
    det_z = xx*yy - xy*xy

    det_max = max(det_x,det_y,det_z)


    #We assume that one of the components is non-zero, find out which one is the best
    if det_max == det_x:
        a = 1.0
        b = (xz*yz - xy*zz) / det_x
        c = (xy*yz - xz*yy) / det_x
    elif det_max == det_y:
        a = (yz*xz - xy*zz) / det_y
        b = 1.0
        c = (xy*xz - yz*xx) / det_y
    else:
        a = (yz*xy - xz*yy) / det_z
        b = (xz*xy - yz*xx) / det_z
        c = 1.0

    #Return normal of plane
    return np.array([a,b,c])
 
def ATOM_DISTANCE(atom1, atom2, simulation_box):
    #Return distance between two atoms, atom1 is reference. Also correct xyz-coordinates for atom2 if across periodic boundaries

    #Check distance in each dimension
    for i in range(3):
        distance = atom1[i]-atom2[i]
        #Periodic box problems?
        if distance > simulation_box[i]/2:
            atom2[i] += simulation_box[i]
        elif distance < -simulation_box[i]/2:
            atom2[i] -= simulation_box[i]

    return sqrt(pow(atom1[0]-atom2[0],2) + pow(atom1[1]-atom2[1],2) + pow(atom1[2]-atom2[2],2))

def VECTOR_NORM(vector):
    return sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2))

def VECTOR_PROJECTION(vector1, vector2):
    #Caclulate the projection of vector1 onto vector2

    #The projection of a vector v1 along the direction of vector v2 is given by
    #Proj_v2(v1) = v2 * scalar product(v1,v2) / (norm(v2))^2
    
    scalar_product = vector1 * vector2
    vector2_norm = VECTOR_NORM(vector2)
    projection = vector2*scalar_product/pow(vector2_norm,2)

    return projection

def COMPUTE_NEIGHBORS(atom, nr_neighbors):
    #Calculates the nr_neighbors number of closest neighbours at this frame to atom with index atom_index
    neighbor_dict = {}
    cutoff = 15 #Angstrom
    neighbors = universe.select_atoms("around %s bynum %s" %(cutoff, atom.index+1))
    
    for neighbor in neighbors:
        distance = ATOM_DISTANCE(atom.position,neighbor.position,atomgroup.dimensions)
        neighbor_dict[int(neighbor.index)] = [distance,neighbor.position]

    #Sort neighbors according to distance to reference-atoms
    neighbor_dict_sorted = sorted(neighbor_dict.items(), key=operator.itemgetter(1))
    closest_neighbors={}

    for i in range(min(nr_neighbors,len(neighbor_dict_sorted))):        
        closest_neighbors[neighbor_dict_sorted[i][0]] = neighbor_dict_sorted[i][1][1]

    return closest_neighbors

def VELOCITY_ALONG_PLANE(atom, plane):
    #Calculate the velocity parallell to the plane
    #The velocity of a vector along the plane is obtained by subtracting the vector
    #component orthogonal to the plane. This is done by subtracting the part of the
    #velocity that is parallell to the normal direction of the plane.
    velocity_along_plane = atom.velocity - VECTOR_PROJECTION(atom.velocity, plane)
    
    return velocity_along_plane
    
######################################
#####  MAIN FUNCTION STARTS HERE  ####
######################################
NUMBER_OF_NEIGHBORS = 5 #Number of neighbors to use in calculation of local plane
pointvector = np.empty((NUMBER_OF_NEIGHBORS,3))
EMPTY_VELOCITY_ARRAY = ["NAN", "NAN", "NAN"]

#Create atom selection
atomgroup = universe.select_atoms('name H1 H2 H3')

data_structure = []

universe.trajectory[0]
if args.calculate:
    with open("velocity_data.txt", "w") as f:
        #Loop over frames    
        for frame in universe.trajectory:
            velocity_data = []
            for atom in atomgroup:
                #Compute NUMBER_OF_NEIGHBORS closest neighbors used to create local plane
                close_neighbors = COMPUTE_NEIGHBORS(atom, NUMBER_OF_NEIGHBORS)
                if len(close_neighbors) < 3:
                    for i in range(3):
                        f.write("%s " %(EMPTY_VELOCITY_ARRAY[i]))
                    velocity_data.append(EMPTY_VELOCITY_ARRAY)
                    continue
                else:
                    #Create pointvector for calculation of local plane
                    for i,neighbor in enumerate(close_neighbors):
                        pointvector[i] = close_neighbors[neighbor]

                    #Calculate normal of local plane
                    plane = PLANE_FROM_POINTS(pointvector)

                    #Calculate velocity parallell to plane
                    velocity = VELOCITY_ALONG_PLANE(atom, plane)
                    for i in range(3):
                        f.write("%s " %(velocity[i]))
                    velocity_data.append(velocity)
            f.write("\n")
            data_structure.append(velocity_data)

origin_reset_counter = 0
diffusion = []
for element in range(50):
    diffusion.append(0)
time_origin = 0

#Read data from file
data_structure = []
frame_data = []
with open("velocity_data.txt", "r") as f:
    atomindex = 0
    counter = 0
    velocity = []
    velocity_lines = f.readlines()    
    for j, frame in enumerate(universe.trajectory):
        data = velocity_lines[j].split()
        for value in data:
            if value != "NAN":
                velocity.append(float(value))
            else:
                velocity.append(value)
            counter += 1
            if counter == 3:
                frame_data.append(velocity)
                velocity = []
                counter = 0
        data_structure.append(frame_data)
        frame_data = []

#print data_structure[0]
#Calculate Green-Kubo diffusion
for frame in universe.trajectory:
    correlation = 0
    ignored_values = 0
    #Calculate correlation values
    for atom in atomgroup:
        if data_structure[frame.frame][atom.index][0] != "NAN":
            if data_structure[time_origin][atom.index][0] != "NAN":
                correlation += np.dot(data_structure[frame.frame][atom.index], data_structure[time_origin][atom.index])
            else:
                ignored_values += 1
        else:
            ignored_values += 1
    #correlation = np.sum(np.array(data_structure[frame.frame])*np.array(data_structure[time_origin]))

    #Fill output-array
    diffusion[origin_reset_counter] += correlation/(len(atomgroup)-ignored_values)

    origin_reset_counter += 1
    #Should we reset counter?
    if origin_reset_counter == 50:
        time_origin = frame.frame + 1
        origin_reset_counter = 0

#Normalize
diffusion[:] = [x / (len(universe.trajectory)/50) for x in diffusion]
diffusion[:] = [x / max(diffusion) for x in diffusion]

with open("diffusion.txt","w") as f:
    for i in range(50):
        f.write("%s  %s\n" %(i, diffusion[i]))

#for atom in atomgroup:
#    print data_structure[atom.index]
#    correlation = np.correlate(data_structure[atom.index], data_structure[atom.index], mode='same')
#    N = len(correlation)
#    autocorrelation = correlation[N//2:]

    #np.correlate uses the unstandardized definition of correlation.                                                                                                                                
    #As the lag gets bigger, the number of points in the overlap between the two signals gets smaller.                                                                                              
    #So the magnitude of the correlations decreases.                                                                                                                                                
    #Correct this by dividing through by the lengths                                                                                                                                                

#    lengths = range(N, N//2, -1)
#    autocorrelation /= lengths

    #Normalize                                                                                                                                                                                      
#    autocorrelation /= autocorrelation[0]
