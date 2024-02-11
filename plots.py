import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
import imageio
import os



def make_gif(path):
    _, _, files = next(os.walk(path))
    images = []
    for filename in files:
        images.append(imageio.imread(path+"/"+filename))
    imageio.mimsave('plots/res.gif', images, duration=5)


skipstep=500




data=pd.read_table('results/vertices.dat', header=None, delimiter=r"\s+")
vertices=np.array(data.loc[:,:])



data_rho=pd.read_table('results/rho.dat', header=None, delimiter=r"\s+")



data_faces=pd.read_table('results/faces.dat', header=None, delimiter=r"\s+", names=['col' + str(x) for x in range(6) ])

maxstep=len(data_rho.loc[:,0])-1


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


min_rho=np.min( data_rho.loc[:maxstep,1:len(x_plot)])
max_rho=np.max( data_rho.loc[:maxstep,1:len(x_plot)])

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
    if(((i+1) % skipstep)==0 ):
        rho=(np.array(data_rho.loc[i,1:len(faces)])-min_rho)/(max_rho-min_rho)
        fig, ax = plt.subplots(figsize=(16, 9), layout='constrained', nrows=2,height_ratios=[15,1])
        fig.suptitle('t='+str(data_rho.loc[i,0]))
        for face_num,face in enumerate(faces):
            ax[0].fill(x_plot_full[face_num], y_plot_full[face_num],facecolor=colorm(rho[face_num]))
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorm),cax=ax[1], orientation='horizontal', label='Density')
        fig.savefig('plots/fig'+"{0:0>3}".format(i)+'.png', bbox_inches='tight')



#fig.show()




rho=np.array(data_rho.loc[maxstep,1:len(faces)])

axis_dist=[]




# for face_center in np.array(face_centers):
#     axis_dist.append(np.sqrt(face_center[0]**2+face_center[1]**2))


# fig=px.scatter(x=axis_dist, y=rho,  labels={"x": "R", "y":"rho"})
# fig.show()


omega=np.array([0,0,2])
face_centers=np.array(face_centers)
rho_analytic=np.exp(-1/3*(np.linalg.norm(omega)**2)*np.sin(-np.arccos(face_centers[:,2]))**3)



for i in range(5):

    rho=np.array(data_rho.loc[maxstep-i*100,1:len(faces)])
    fig=px.scatter(x=theta_fc, y=rho_analytic,  labels={"x": "theta", "y":"rho"})
    fig.update_traces(marker=dict(color='red'))
    fig.add_traces(list(px.scatter(x=theta_fc, y=rho,  labels={"x": "theta", "y":"rho"}).select_traces()))
    fig.update_layout(title_text="t="+str((maxstep-i*100)*0.002),showlegend=False)
    fig.update_yaxes(range = [0.5,4.3])
    fig.show()


rho=np.array(data_rho.loc[0,1:len(faces)])
fig=px.scatter(x=theta_fc, y=rho_analytic,  labels={"x": "theta", "y":"rho"})
fig.update_traces(marker=dict(color='red'))
fig.add_traces(list(px.scatter(x=theta_fc, y=rho,  labels={"x": "theta", "y":"rho"}).select_traces()))
fig.update_layout(title_text="t="+str((0)*0.002),showlegend=False)
fig.update_yaxes(range = [0.5,4.3])
fig.show()


fig=px.scatter(x=theta_fc, y=rho,  labels={"x": "theta", "y":"rho"})
fig.show()




# mesh_step=[1,2,3,4,5,1,2,3,4,2,3,4,5]
# mesh_type=['quad','quad','quad','quad','quad','ico','ico','ico','ico','hex','hex','hex','hex']
# mesh_errs=[0.00012,5e-6,0.0018,0.0011,0.00055,0.0003,0.0003,0.0004,0.0003, 0.00014,2.33e-5, 2.8e-5,1.67e-5]


# df=pd.DataFrame(data={'mesh_step':mesh_step, 'mesh_type':mesh_type, 'mesh_errs': mesh_errs})


# fig=px.scatter(df, x='mesh_step', y='mesh_errs', 
#                color='mesh_type', title='Mesh rho errors depending on type on 1 rotation', log_y=True)
# fig.update_layout(font=dict(size=30))
# fig.update_traces(mode="lines+markers")
# fig.update_traces(marker_size=15)
# fig.show()


# # make_gif("plots/quad_4_true_longer")



data=pd.read_table('errors_p.txt', header=None, delimiter=r"\s+")

T_rot=2*np.pi/8

fig=px.scatter(x=data[0]/T_rot, y=data[1],  labels={"x": "N_of_rotations", "y":"relative rho error"}, log_y=True)
fig.update_layout(font=dict(size=30))
fig.show()


