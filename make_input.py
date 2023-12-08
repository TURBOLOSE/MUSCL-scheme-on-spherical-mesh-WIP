import numpy as np
import pandas as pd


dim = 4
face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
N=len(face_centers[0])


face_centers=np.array(face_centers)

rho=np.ones(N)
l=[]


omega=np.array([0,0,1])

for face_num, R in enumerate(face_centers):
    #if( R[2] >=0):
       #omega=np.array([0,0,1])
       #rho[face_num]*=2
    #else:
        #omega=np.array([0,0,-1])
    l.append(rho[face_num]*np.cross(R,np.cross(omega,R)))


l=np.array(l)




pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2]]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)



