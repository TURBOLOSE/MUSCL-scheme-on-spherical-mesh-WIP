import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import plotly.express as px
from tqdm import tqdm
from scipy.interpolate import griddata


def projection_plots(value, print_residuals:bool=False, print_log:bool=False): #value = rho,p,omega
    skipstep=1
    
    gam=1.25

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
        #label_pr='Bernoulli integral -1 /R'
    elif(value=='c_s'):
        data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")
        data_p=pd.read_table('results/p.dat', header=None, delimiter=r"\s+")
        label_pr='Speed of sound'
        data_rho.loc[:,1:]=data_p.loc[:,1:]/data_rho.loc[:,1:]
        data_rho.loc[:,1:]=np.sqrt(1.25*data_rho.loc[:,1:])
    elif(value=='vel_abs'):
        label_pr='Speed'
        data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")
        data_Lx=pd.read_table('results/Lx.dat', header=None, delimiter=r"\s+")
        data_Ly=pd.read_table('results/Ly.dat', header=None, delimiter=r"\s+")
        data_Lz=pd.read_table('results/Lz.dat', header=None, delimiter=r"\s+")
        face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")
        maxstep=len(data_rho.loc[:,0])
        n_faces=len(data_rho.loc[0,:])-1
        for i in range(maxstep):
            for j in range(1,n_faces):
                L=np.array([data_Lx.loc[i,j],data_Ly.loc[i,j],data_Lz.loc[i,j]])
                fc=np.array(face_centers.loc[j-1]/np.linalg.norm(face_centers.loc[j-1]))
                rho0=data_rho.loc[i,j]
                data_rho.loc[i,j]=np.linalg.norm(np.cross(fc,L))/rho0
    # elif(value=='mach'):
    #     data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")
    #     data_p=pd.read_table('results/p.dat', header=None, delimiter=r"\s+")
    #     maxstep=len(data_rho.loc[:,0])
    #     for i in range(1,maxstep):
    #         data_rho.loc[i,:]=np.sqrt(gam*data_p.loc[0,:]/data_rho.loc[0,:])

    #     label_pr='Mach number'
    else:
        print("wrong type of plot value")
        return

    maxstep=len(data_rho.loc[:,0])

    if(print_residuals):
        for i in range(1,maxstep):
            data_rho.loc[i,:]-=data_rho.loc[0,:]
        data_rho.loc[0,:]-=data_rho.loc[0,:]
        label_pr+=" residuals"

    if(print_log):
        data_rho.loc[:,1:]=np.log10(data_rho.loc[:,1:])
        
        
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


        patches = []

    for face_num,face in enumerate(faces):
        polygon = Polygon(np.vstack([x_plot_full[face_num], y_plot_full[face_num]]).T,closed=True)
        patches.append(polygon)

   
    #=====================================================
    # face_centers=np.array(face_centers)
    # omega=np.array([0,0,2])
    # #for i in range(maxstep):
    # #    data_rho.loc[i,1:len(faces)]-=1+(np.linalg.norm(omega)**2*1./2*np.sin(-np.arccos(face_centers[:,2]))**2)
                    
    # data_rho.loc[:,1:len(faces)]=data_p.loc[:,1:len(faces)]/data_rho.loc[:,1:len(faces)]**(1.4)
    #=====================================================
    colorm = plt.get_cmap('viridis')
    min_rho=np.min( data_rho.loc[:maxstep,1:len(x_plot)])
    max_rho=np.max( data_rho.loc[:maxstep,1:len(x_plot)])
    #min_rho=0
    #max_rho=2e-4
    max_rho=1.8

    norm = mpl.colors.Normalize(vmin=min_rho, vmax=max_rho)
    mpl.rcParams.update({'font.size': 22})

   


    for i in tqdm(range(maxstep)): #dens
        if((i % skipstep)==0 ):

            collection = PatchCollection(patches)
            colors=colorm(norm(data_rho.loc[i,1:len(faces)]))

            fig, ax = plt.subplots(figsize=(16, 10), layout='constrained', nrows=2,height_ratios=[15,1])
            #fig.tight_layout()
            plt.subplots_adjust(hspace=10)
            #rho=(np.array(data_rho.loc[i,1:len(faces)])-min_rho)/(max_rho-min_rho)
            fig.suptitle('t='+"{:.4f}".format(data_rho.loc[i,0]*3.3e-5))
            ax[0].set_xlabel(r'$\lambda / \sqrt{2}$', fontsize=25)
            ax[0].set_ylabel(r'$\sqrt{2}  \sin(\varphi )$', fontsize=25)
            #for face_num,face in enumerate(faces):
                #ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]),edgecolor =colorm(rho[face_num]))
                #ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]),edgecolor =colorm(rho[face_num]))
            #ax[0].fill(x_plot_full, y_plot_full,facecolor=colorm(norm(data_rho.loc[i,1:len(faces)])), edgecolor=colorm(norm(data_rho.loc[i,1:len(faces)])))
            ax[0].add_collection(collection)
            collection.set_color(colors)
            #ax[0].autoscale_view()
            ax[0].set_xlim([-2.5, 3.4])
            ax[0].set_ylim([-1.5, 1.5])

            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),cax=ax[1], orientation='horizontal', label=label_pr)
            fig.savefig('plots/fig'+"{0:0>4}".format(i)+'.png', bbox_inches='tight')
            plt.clf()
            plt.close()
    



projection_plots("rho", print_residuals=False, print_log=False)



def light_curve(data_p, face_centers):


    observer_vector=np.array([9.4*10**3*3*10**16/15000,0,0]) #dist=9400pc (in R_ns)
    fc=np.array(face_centers)
    #theta_fc=-np.arccos(fc[:,2]/np.linalg.norm(fc, axis=1))+np.pi/2
    phi_fc=np.arctan2(fc[:,1]/np.linalg.norm(fc, axis=1),fc[:,0]/np.linalg.norm(fc, axis=1))
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



def vel_plot():
    gam=1.25

    #turn_angle=np.pi/2

    turn_angle=0
    path_to_res='results/'
    #path_to_res='plots/big_sim/'
    #path_to_res='plots/0.4c crashes/time series/'
    

    data_full=pd.read_table(path_to_res+'final_state.dat', header=None, delimiter=r"\s+")
    data_faces=pd.read_table(path_to_res+'faces.dat', header=None, delimiter=r"\s+", names=['col' + str(x) for x in range(6) ])
    face_centers=pd.read_table(path_to_res+'face_centers.dat', header=None, delimiter=r"\s+")
    data=pd.read_table(path_to_res+'vertices.dat', header=None, delimiter=r"\s+")
    vertices=np.array(data.loc[:,:])
    faces=np.array(data_faces.loc[:,:])

    vel=[]
    p=[]


    turn_matrix=np.matrix([[np.cos(turn_angle),0, np.sin(turn_angle)],[0,1,0],[-np.sin(turn_angle),0, np.cos(turn_angle)]])


    for num, el in enumerate(data_full.loc[:,4]):
        l=np.array([data_full.loc[num,1],data_full.loc[num,2],data_full.loc[num,3]])
        R=np.array(face_centers.loc[num])
        vel.append(np.cross(R,l)/(-np.linalg.norm(R)*data_full.loc[num,0]))
        p.append(data_full.loc[num,0])
        #p.append((data_full.loc[num,4]-data_full.loc[num,0]/2 * np.linalg.norm(vel[num])**2)*(gam-1))


    for ver_num, vertice in enumerate(vertices):
        vertices[ver_num]=np.matmul(turn_matrix,vertices[ver_num])

    #for face_num, face in enumerate(faces):
    #    faces[face_num]=np.matmul(turn_matrix,faces[face_num])


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


    vel=np.array(vel)

    for vel_num, ve in enumerate(vel):
            vel[vel_num]=np.matmul(turn_matrix,vel[vel_num])

    face_centers=np.array(face_centers)

    for face_cent_num, face_cent in enumerate(face_centers):
            face_centers[face_cent_num]=np.matmul(turn_matrix,face_centers[face_cent_num])



    theta_fc=np.arccos(face_centers[:,2]/np.linalg.norm(face_centers, axis=1))
    phi_fc=np.arctan2(face_centers[:,1]/np.linalg.norm(face_centers, axis=1),
                    face_centers[:,0]/np.linalg.norm(face_centers, axis=1))

    x_fc=phi_fc/np.sqrt(2)
    y_fc=np.sin(-theta_fc+np.pi/2)*np.sqrt(2)

    for face_num,face in enumerate(faces): #fix x
        sign_arr=np.sign(x_plot_full[face_num])
        if( (not (0 in sign_arr)) and (1 in sign_arr) and (-1 in sign_arr) and (np.min(np.abs(x_plot_full[face_num])) > 1)):
            for i,element in enumerate(x_plot_full[face_num]):
                if(element < 0):
                    x_plot_full[face_num][i]+=2*np.pi/np.sqrt(2)
            #x_fc[face_num]+=2*np.pi/np.sqrt(2)



    rd=np.sin(theta_fc)*np.cos(phi_fc)*vel[:,0]+np.sin(theta_fc)*np.sin(phi_fc)*vel[:,1]+np.cos(theta_fc)*vel[:,2]
    theta_d=np.cos(theta_fc)*np.cos(phi_fc)*vel[:,0]+np.cos(theta_fc)*np.sin(phi_fc)*vel[:,1]-np.sin(theta_fc)*vel[:,2]
    phi_d=(-np.sin(phi_fc)*vel[:,0]+np.cos(phi_fc)*vel[:,1])/np.sin(theta_fc)

    for num,ph in enumerate(phi_d):
        if(ph>1e5 or np.isinf(ph)):
            phi_d[num]=0



    xd=phi_d/np.sqrt(2)
    yd=theta_d*np.sqrt(2)*np.cos(theta_fc)


    #ax[0].quiver(x_fc[::3], y_fc[::3],xd[::3],yd[::3])
    X_gr, Y_gr=np.meshgrid(np.linspace(-2.2,2.2, 100),np.linspace(-1.4, 1.4, 100))

    mask=np.logical_or(np.isnan(xd, where=False),np.isnan(xd, where=False))


    xd_gr=griddata(np.stack([x_fc[mask].T, y_fc[mask].T]).T, xd[mask],(X_gr,Y_gr), method='nearest')
    yd_gr=griddata(np.stack([x_fc[mask].T, y_fc[mask].T]).T, yd[mask],(X_gr,Y_gr), method='nearest')    








    ##======================================================
    # order=np.argsort(x_fc)
    # x_fc=np.sort(x_fc)
    # y_fc=y_fc[order]

    # xd=xd[order]
    # yd=yd[order]

    # plt.scatter(x_fc, y_fc)
    # plt.show()


    # x_fc_mg, y_fc_mg=np.meshgrid(x_fc, y_fc)
    # xd_fc_mg, yd_fc_mg=np.meshgrid(xd, yd)



    # plt.streamplot(x_fc_mg, y_fc_mg,xd_fc_mg, yd_fc_mg)
    # plt.show()
    # #======================================================



    colorm = plt.get_cmap('viridis')
    min_p=np.min(p)
    max_p=np.max(p)

    norm = mpl.colors.Normalize(vmin=min_p, vmax=max_p)
    mpl.rcParams.update({'font.size': 22})



    colorm2 = plt.get_cmap('inferno')
    v=np.sqrt(xd_gr**2+yd_gr**2)

    #print(xd_gr,xd_gr)
    #print(np.min(v),np.max(v))
    norm2 = mpl.colors.Normalize(vmin=np.min(v), vmax=np.max(v))



    fig, ax = plt.subplots(figsize=(18, 10), layout='constrained', nrows=3,height_ratios=[16,1,1])

    plt.subplots_adjust(hspace=10)
    rho=(np.array(p-min_p)/(max_p-min_p))
    fig.suptitle('t='+"{:10.4f}".format(1.4235))
    ax[0].set_xlabel(r'$\varphi / \sqrt{2}$', fontsize=25)
    ax[0].set_ylabel(r'$\sqrt{2}  \sin(\theta )$', fontsize=25)
    for face_num,face in enumerate(faces):
        ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]),edgecolor =colorm(rho[face_num]))

    #ax[0].quiver(x_fc[::3], y_fc[::3],xd[::3],yd[::3])
        #if(face_num % 20 == 0):
    #for face_num,face in enumerate(faces):
        #ax[0].arrow(x_fc[face_num],y_fc[face_num],1e-1*xd[face_num],1e-1*yd[face_num],width=0.007, color='grey', alpha=0.9)

    #ax[0].streamplot(X_gr,Y_gr,xd_gr,yd_gr,color=np.sqrt(xd_gr*2*+yd_gr**2), arrowsize=3)
    ax[0].streamplot(X_gr,Y_gr,xd_gr,yd_gr,color=v,norm=norm2, cmap=colorm2, arrowsize=3)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),cax=ax[1], orientation='horizontal', label="Density")
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm2, cmap=colorm2),cax=ax[2], orientation='horizontal', label="Speed")
    fig.savefig('plots/vel_plot_1.png', bbox_inches='tight',dpi=400)
    plt.clf()
    plt.close()



vel_plot()

