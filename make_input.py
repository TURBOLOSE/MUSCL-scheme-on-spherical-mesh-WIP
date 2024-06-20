import numpy as np
import pandas as pd


def make_input_4(): #no energy as separate variable
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)


    l=[]



    theta_face_centers=-np.arccos(face_centers[:,2]/np.linalg.norm(face_centers, axis=1))+np.pi/2


    #rho=np.ones(N)
    omega=np.array([0,0,2])

    #rho=np.ones(N)
    r=np.sqrt(face_centers[:,1]**2+face_centers[:,2]**2+face_centers[:,0]**2)
    rho=np.exp(-1/2*(np.linalg.norm(omega)**2)*np.sin(-np.arccos(face_centers[:,2]/r)+np.pi/2)**2)

    for face_num, R in enumerate(face_centers):
        #if( R[2] >0):
        #   omega=np.array([0,0,0.5])
        #elif ( R[2] <0):
        #    omega=np.array([0,0,-0.5])
        #else:
        #    omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
    l=np.array(l)

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2]]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)

def make_input_4_cos_bell(): #no energy as separate variable
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)


    l=[]
    theta=-np.arccos(face_centers[:,2])+np.pi/2 #lat
    phi=np.arctan2(face_centers[:,1],face_centers[:,0]) #long

    rho=np.ones(N)
    omega=np.array([0,0,0])




    for face_num, R in enumerate(face_centers):
        if( abs(theta[face_num])<np.pi/6):
            if(abs(phi[face_num])<np.pi/6):
                rho[face_num]=1+10*np.cos(2.2*np.sqrt(theta[face_num]**2+phi[face_num]**2))
                if(rho[face_num]<0):
                    print('check input rho')
                omega=np.array([0,0,2])

        #elif ( R[2] <0):
        #    omega=np.array([0,0,-2])
        #else:
        #    omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        omega=np.array([0,0,0])
    l=np.array(l)

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2]]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)



def make_input_5():
    gam0=4./3
    gam=2-1/gam0
    #face_centers=pd.read_table('results/face_centers_ico_6.dat', header=None, delimiter=r"\s+")
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)

    rho=np.ones(N)
    omega=np.array([0,0,0])
    p=np.ones(N)
    #p=1+(np.linalg.norm(omega)**2*rho/2*np.sin(-np.arccos(face_centers[:,2]))**2)
    l=[]
    v=[]

    
    theta=-np.arccos(face_centers[:,2]/np.linalg.norm(face_centers, axis=1))+np.pi/2 #lat
    phi=np.arctan2(face_centers[:,1]/np.linalg.norm(face_centers, axis=1),face_centers[:,0]/np.linalg.norm(face_centers, axis=1)) #long
    a=1
    R0=a/3
    #r=a*np.arccos(np.sin(0)*np.sin(theta)+np.cos(theta)*np.cos(0)*np.cos(phi-np.pi*3/2))
    theta_c=0.3
    phi_c=0.5
    r=a*np.arccos(np.sin(theta_c)*np.sin(theta)+np.cos(theta)*np.cos(theta_c)*np.cos(phi-phi_c))


    for face_num, R in enumerate(face_centers):
        #if( R[2] >0):
        #    omega=np.array([0,0,0.5])
        #elif ( R[2] <0):
        #    omega=np.array([0,0,-0.5])
        #else:
        #    omega=np.array([0,0,0])
        #if(np.abs(theta[face_num]-np.pi/2)<0.1):
        #   p[face_num]=5
        if r[face_num]<R0:
            p[face_num]=50
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R)/np.linalg.norm(R))
        #print(np.cross(R,-np.cross(R, l[face_num]))/(np.linalg.norm(R)**2))

        
    l=np.array(l)
    v=np.array(v)

    #E=1/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2
    E=gam/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2
    if(E.any()<0):
        print('Energy<0!!')

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)

def make_input_5_sp_layer():
    gam0=4./3
    gam=2-1/gam0

    #V_conv=0.01
    V_conv=1

    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])

    rho=np.ones(N) #10^7 g/cm^2
    c_s=2*10**(-3)/V_conv #c
    p=c_s**2*np.ones(N)/gam
    print(p[0])
    omega_ns=0.01/V_conv#c
    omega=np.array([0,0,0.01])/V_conv #c
    

    
    #face_centers=pd.read_table('results/face_centers_ico_6.dat', header=None, delimiter=r"\s+")
    


    face_centers=np.array(face_centers)

    
    

    
    #p=1+(np.linalg.norm(omega)**2*rho/2*np.sin(-np.arccos(face_centers[:,2]))**2)
    l=[]
    v=[]

    theta=np.arccos(face_centers[:,2]/np.linalg.norm(face_centers, axis=1)) 

    for face_num, R in enumerate(face_centers):
        #if( R[2] >0):
        #    omega=np.array([0,0,0.5])
        #elif ( R[2] <0):
        #    omega=np.array([0,0,-0.5])
        #else:
        #    omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R)/np.linalg.norm(R))


        #print(np.linalg.norm(v0)**2 -np.linalg.norm(v[face_num])**2)
        #print(np.cross(R,-np.cross(R, l[face_num]))/(np.linalg.norm(R)**2))

    
    l=np.array(l)
    v=np.array(v)
    #E=gam/(gam-1)*p+rho*((np.linalg.norm(v, axis=1)**2)/2-np.ones(N)*omega_ns**2*(np.sin(theta)**2)/2)
    E=1/(gam-1)*p+rho*((np.linalg.norm(v, axis=1)**2)/2-np.ones(N)*omega_ns**2*(np.sin(theta)**2)/2)
    #E=gam/(gam-1)*p+rho*((np.linalg.norm(v, axis=1)**2)/2)
    if(E.any()<0):
        print('Energy<0!!')

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False, float_format="%.15f")



def make_input_5_new_p():
    gam0=4./3
    gam=2-1/gam0
    
    #face_centers=pd.read_table('results/face_centers_ico_6.dat', header=None, delimiter=r"\s+")
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])


    face_centers=np.array(face_centers)

    rho=np.ones(N)
    omega_ns=5
    omega=np.array([0,0,5])
    

    p=np.ones(N)
    #p=1+(np.linalg.norm(omega)**2*rho/2*np.sin(-np.arccos(face_centers[:,2]))**2)
    l=[]
    v=[]

    theta=np.arccos(face_centers[:,2]/np.linalg.norm(face_centers, axis=1)) 

    for face_num, R in enumerate(face_centers):
        #if( R[2] >0):
        #    omega=np.array([0,0,0.5])
        #elif ( R[2] <0):
        #    omega=np.array([0,0,-0.5])
        #else:
        #    omega=np.array([0,0,0])
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R)/np.linalg.norm(R))


        #print(np.linalg.norm(v0)**2 -np.linalg.norm(v[face_num])**2)
        #print(np.cross(R,-np.cross(R, l[face_num]))/(np.linalg.norm(R)**2))

    
    l=np.array(l)
    v=np.array(v)

    #E=gam/(gam-1)*p+rho*((np.linalg.norm(v, axis=1)**2)/2-np.ones(N)*omega_ns**2*(np.sin(theta)**2)/2)
    E=gam/(gam-1)*p+rho*((np.linalg.norm(v, axis=1)**2)/2)
    if(E.any()<0):
        print('Energy<0!!')


    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False, float_format="%.15f")


def make_input_5_cos_bell():  #adds energy
    gam=1.4
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])

    face_centers=np.array(face_centers)


    l=[]
    
    theta=-np.arccos(face_centers[:,2])+np.pi/2 #lat
    phi=np.arctan2(face_centers[:,1],face_centers[:,0]) #long


    omega=np.array([0,0,0.5])
    rho=np.ones(N)
    



    v=[]
    a=1
    R0=a/3
    r=a*np.arccos(np.sin(0)*np.sin(theta)+np.cos(theta)*np.cos(0)*np.cos(phi-np.pi*3/2))

    for face_num, R in enumerate(face_centers):
        
        
        if r[face_num]<R0:
            rho[face_num]=1+0.5*np.cos(np.pi*r[face_num]/R0)


        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R/np.linalg.norm(R)))
    

    l=np.array(l)
    v=np.array(v)

    p=1+(np.linalg.norm(omega)**2*rho/2*np.sin(-np.arccos(face_centers[:,2]))**2)
    E=1/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)

def make_input_5_const_entr():  


    gam0=4/3
    gam=2-1/gam0
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])

    face_centers=np.array(face_centers)

    l=[]
    theta=-np.arccos((face_centers[:,2])/np.linalg.norm(face_centers, axis=1)) 

    omega=np.array([0,0,2])
    rho_0=1
    p_0=1
    a_0=np.sqrt(gam*p_0/rho_0)

    M_0=np.linalg.norm(omega)/a_0

    #M_0=(gam-1)/p_0 *np.linalg.norm(omega)**2

    rho=rho_0*(1+(gam-1)/2*M_0**2*np.sin(theta)**2)**(1/(gam-1))
    p=p_0*(1+(gam-1)/2*M_0**2*np.sin(theta)**2)**(gam/(gam-1))



    v=[]


    for face_num, R in enumerate(face_centers):
        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R/np.linalg.norm(R)))
    

    l=np.array(l)
    v=np.array(v)

    
    E=gam/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)


def make_input_6_cos_bell(): 
    gam=1.4
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
    N=len(face_centers[0])

    face_centers=np.array(face_centers)
    
    
    

    omega=np.array([0,0,0.5])


    rho=np.ones(N)
    tracer=np.zeros(N)


    l=[]
    theta_0=-np.arccos(face_centers[:,2]) 

    rho_0=1
    p_0=1
    a_0=np.sqrt(gam*p_0/rho_0)

    M_0=np.linalg.norm(omega)/a_0

    #M_0=(gam-1)/p_0 *np.linalg.norm(omega)**2

    rho=rho_0*(1+(gam-1)/2*M_0**2*np.sin(theta_0)**2)**(1/(gam-1))
    p=p_0*(1+(gam-1)/2*M_0**2*np.sin(theta_0)**2)**(gam/(gam-1))


    v=[]
    a=0.05
    R0=a/3
    phi=np.arctan2(face_centers[:,1],face_centers[:,0]) #long
    theta=-np.arccos(face_centers[:,2])+np.pi/2 #lat
    r=a*np.arccos(np.sin(0)*np.sin(theta)+np.cos(theta)*np.cos(0)*np.cos(phi-np.pi*3/2))

    for face_num, R in enumerate(face_centers):
        if r[face_num]<R0:
            tracer[face_num]=500*(1+1*np.cos(np.pi*r[face_num]/R0))


        l.append(rho[face_num]*np.cross(R,np.cross(omega,R))/(np.linalg.norm(R)**2))
        v.append(np.cross(omega,R/np.linalg.norm(R)))
    

    l=np.array(l)
    v=np.array(v)

    

    E=1/(gam-1)*p+rho*np.linalg.norm(v, axis=1)*np.linalg.norm(v, axis=1)/2

    pd.DataFrame(data=np.array([rho, l[:,0],l[:,1],l[:,2],E, tracer]).transpose()).to_csv('input/input.dat',index=False, sep=' ', header=False)



#make_input_5_new_p()
#make_input_5()
#make_input_5_const_entr()
make_input_5_sp_layer()



#make_input_4()
#make_input_6_cos_bell()




