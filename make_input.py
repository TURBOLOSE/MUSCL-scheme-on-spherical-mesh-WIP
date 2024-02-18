import numpy as np
import pandas as pd


def make_input_4(): #no energy as separate variable
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)


    l=[]

  

    theta_face_centers=-np.arccos(face_centers[:,2])+np.pi/2


    #rho=np.ones(N)
    omega=np.array([0,0,2])
    rho=np.exp(-1/3*(np.linalg.norm(omega)**2)*np.sin(-np.arccos(face_centers[:,2]))**3)


    for face_num, R in enumerate(face_centers):
        #if( R[2] >0):
        #    omega=np.array([0,0,2])
        #elif ( R[2] <0):
        #    omega=np.array([0,0,-2])
        #else:
        #    omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
    l=np.array(l)

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2]]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)


def make_input_5():  #adds energy
    gam=1.4
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)

    #rho=np.ones(N)



    p=np.ones(N)

    #rho[3]=3
    #p[3]=3

    l=[]

    omega=np.array([0,0,0])
    rho=np.ones(N)

    #rho=np.exp(-(np.linalg.norm(omega)**2)*np.sin(np.arccos(face_centers[:,2]))**3)


    for face_num, R in enumerate(face_centers):
        if( R[2] >0):
            omega=np.array([0,0,2])
        elif ( R[2] <0):
            omega=np.array([0,0,-2])
        else:
            omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
    l=np.array(l)
    
    v=np.cross(omega,np.array(face_centers))/np.array(
        [np.linalg.norm(face_centers, axis=1),np.linalg.norm(face_centers, axis=1),np.linalg.norm(face_centers, axis=1)]
        ).transpose()
    E=1/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)







make_input_4()






