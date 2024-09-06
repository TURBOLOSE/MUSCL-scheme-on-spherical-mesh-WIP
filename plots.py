import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
from tqdm import tqdm


def projection_plots(value, print_residuals:bool=False): #value = rho,p,omega
    skipstep=1
    

    if(value=='rho'):
        data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")
        label_pr='Density'
    elif(value=='p'):
        data_rho=pd.read_table('results/p.dat', header=None, delimiter=r"\s+")
        label_pr='Pressure'
    elif(value=='omega'):
        data_rho=pd.read_table('results/omega.dat', header=None, delimiter=r"\s+")
        label_pr='Omega_z'
    elif(value=='vort'):
        data_rho=pd.read_table('results/curl.dat', header=None, delimiter=r"\s+")
        label_pr='Vorticity'
    else:
        print("wrong type of plot value")
        return

    maxstep=len(data_rho.loc[:,0])

    if(print_residuals):
        for i in range(1,maxstep):
            data_rho.loc[i,:]-=data_rho.loc[0,:]
        data_rho.loc[0,:]-=data_rho.loc[0,:]
        label_pr+=" residuals"
        
        
    #data_p=pd.read_table('results/p.dat', header=None, delimiter=r"\s+")
    #data_omega=pd.read_table('results/omega.dat', header=None, delimiter=r"\s+")
   



    data_faces=pd.read_table('results/faces.dat', header=None, delimiter=r"\s+", names=['col' + str(x) for x in range(6) ])
    face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")

    data=pd.read_table('results/vertices.dat', header=None, delimiter=r"\s+")
    vertices=np.array(data.loc[:,:])
    faces=np.array(data_faces.loc[:,:])


    #==============================================================================================

    # theta=-np.arccos(np.array(face_centers)[:,2]/np.linalg.norm(np.array(face_centers), axis=1)) 
    # gam=2-3./4
    # omega=np.array([0,0,5])
    # rho_0=1
    # p_0=1
    # a_0=np.sqrt(gam*p_0/rho_0)
    # M_0=np.linalg.norm(omega)/a_0
    # rho_aa=rho_0*(1+(gam-1)/2*M_0**2*np.sin(theta)**2)**(1/(gam-1))



    # for i in range(maxstep):
    #     data_rho.loc[i,1:len(faces)]=(data_rho.loc[i,1:len(faces)]-rho_aa)/rho_aa
    #     #data_rho.loc[i,1:len(faces)]=data_p.loc[i,1:len(faces)]/data_rho.loc[i,1:len(faces)]**gam

    #==============================================================================================

    faces_new=[]

    for face_num, face in enumerate(faces): #trick for variable length of each face (needed for hex meshes)
        faces_new.append(face[~np.isnan(face)].astype(int))

    faces=faces_new


    theta=-np.arccos(vertices[:,2])+np.pi/2
    phi=np.arctan2(vertices[:,1],vertices[:,0])


    x_plot=phi/(np.sqrt(2)) #projection
    y_plot=np.sqrt(2)*np.sin(theta)


    x_plot_full=[]
    y_plot_full=[]


    for face in faces:
        temp_x=[]
        temp_y=[]
        for face_el in face:
            temp_x.append(x_plot[face_el])
            temp_y.append(y_plot[face_el])
        x_plot_full.append(temp_x)
        y_plot_full.append(temp_y)




    for face_num,face in enumerate(faces): #fix x
        sign_arr=np.sign(x_plot_full[face_num])
        if( (not (0 in sign_arr)) and (1 in sign_arr) and (-1 in sign_arr) and (np.min(np.abs(x_plot_full[face_num])) > 1)):
            for i,element in enumerate(x_plot_full[face_num]):
                if(element < 0):
                    x_plot_full[face_num][i]+=2*np.pi/np.sqrt(2)



    #=====================================================
    # face_centers=np.array(face_centers)
    # omega=np.array([0,0,2])
    # #for i in range(maxstep):
    # #    data_rho.loc[i,1:len(faces)]-=1+(np.linalg.norm(omega)**2*1./2*np.sin(-np.arccos(face_centers[:,2]))**2)
                    
    # data_rho.loc[:,1:len(faces)]=data_p.loc[:,1:len(faces)]/data_rho.loc[:,1:len(faces)]**(1.4)
    #=====================================================
    colorm = plt.get_cmap('viridis')
    min_rho=np.min( data_rho.loc[:maxstep,1:len(x_plot)])
    #max_rho=np.max( data_rho.loc[:maxstep,1:len(x_plot)])
    max_rho=4.5

    norm = mpl.colors.Normalize(vmin=min_rho, vmax=max_rho)
    mpl.rcParams.update({'font.size': 22})

    
    for i in tqdm(range(maxstep)): #dens
        if((i % skipstep)==0 ):
            fig, ax = plt.subplots(figsize=(16, 10), layout='constrained', nrows=2,height_ratios=[15,1])
            #fig.tight_layout()
            plt.subplots_adjust(hspace=10)
            rho=(np.array(data_rho.loc[i,1:len(faces)])-min_rho)/(max_rho-min_rho)
            fig.suptitle('t='+"{:10.4f}".format(data_rho.loc[i,0]))
            ax[0].set_xlabel(r'$\lambda / \sqrt{2}$', fontsize=25)
            ax[0].set_ylabel(r'$\sqrt{2}  \sin(\varphi )$', fontsize=25)
            for face_num,face in enumerate(faces):
                ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]),edgecolor =colorm(rho[face_num]))
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),cax=ax[1], orientation='horizontal', label=label_pr)
            fig.savefig('plots/fig'+"{0:0>4}".format(i)+'.png', bbox_inches='tight')
            plt.clf()
            plt.close()
    




projection_plots("vort", print_residuals=False)




def light_curve(data_p, face_centers):


    observer_vector=np.array([9.4*10**3*3*10**16/15000,0,0]) #dist=9400pc (in R_ns)
    fc=np.array(face_centers)
    #theta_fc=-np.arccos(fc[:,2]/np.linalg.norm(fc))+np.pi/2
    phi_fc=np.arctan2(fc[:,1]/np.linalg.norm(fc),fc[:,0]/np.linalg.norm(fc))
    flux=[]
    t=np.array(data_p.loc[:,0])


    for n_step,t_step in enumerate(t):
        flux.append(0)
        for face_num,face_center in enumerate(fc):
            if(phi_fc[face_num] <np.pi/2 and phi_fc[face_num] >=-np.pi/2  ):
                d_vec=np.dot(observer_vector,face_center)
                cos_alpha=np.linalg.norm(d_vec)/(np.linalg.norm(observer_vector)*np.linalg.norm(face_center))
                flux[n_step]+=data_p.loc[n_step,1+face_num]*cos_alpha
    
    fig=px.line(x=t*3.33*10**(-5), y=flux,  labels={"x": "t, sec", "y":"Flux"})
    fig.update_layout(font=dict(size=40))
    fig.show()





# face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
# data_p=pd.read_table('results/p.dat', header=None, delimiter=r"\s+")
# light_curve(data_p, face_centers)




    

