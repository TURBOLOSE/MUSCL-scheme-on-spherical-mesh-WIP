import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import imageio
import os



def make_gif(path):
    _, _, files = next(os.walk(path))
    images = []
    for filename in files:
        images.append(imageio.imread(path+"/"+filename))
    imageio.mimsave('plots/res.gif', images, duration=5)





maxstep=600
skipstep=100


data=pd.read_table('results/vertices.dat', header=None, delimiter=r"\s+")
vertices=np.array(data.loc[:,:])



data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")

data_faces=pd.read_table('results/faces.dat', header=None, delimiter=r"\s+", names=['col' + str(x) for x in range(6) ])

faces=np.array(data_faces.loc[:,:])

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





colorm = plt.get_cmap('viridis')


min_rho=np.min( data_rho.loc[:,1:len(x_plot)])
max_rho=np.max( data_rho.loc[:,1:len(x_plot)])

[min_rho,max_rho]
norm = mpl.colors.Normalize(vmin=min_rho, vmax=max_rho)




# rho=(np.array(data_rho.loc[i*skipstep,1:len(faces)])-min_rho)/(max_rho-min_rho)
# fig, ax = plt.subplots(figsize=(16, 9), layout='constrained', nrows=2,height_ratios=[15,1])


# for face_num,face in enumerate(faces):
#     ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]))
# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),
#              cax=ax[1], orientation='horizontal', label='Density')

face_centers=pd.read_table('results/face_centers.dat', header=None, delimiter=r"\s+")

theta_fc=-np.arccos(face_centers.loc[:,2])+np.pi/2
np.max(np.abs(np.array(data_rho.loc[0,1:len(faces)])-np.cos(theta_fc)))


for i in range(maxstep):
    if(i % skipstep==0 ):
        rho=(np.array(data_rho.loc[i,1:len(faces)])-min_rho)/(max_rho-min_rho)
        fig, ax = plt.subplots(figsize=(16, 9), layout='constrained', nrows=2,height_ratios=[15,1])
        fig.suptitle('t='+str(data_rho.loc[i,0]))
        for face_num,face in enumerate(faces):
            ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]))
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),cax=ax[1], orientation='horizontal', label='Density')
        fig.savefig('plots/fig'+"{0:0>3}".format(i)+'.png', bbox_inches='tight')

#fig.show()



make_gif("plots/quad_4_true_longer")

