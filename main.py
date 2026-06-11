import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import LinearLocator
from PIL import Image
import time
import os

start=time.time()                   #start the time count for the numerical method computation

# Create output folder if it doesn't exist
os.makedirs('output', exist_ok=True)

# In this section I am setting the domain of solution and the discretised grid
A=0                                 #starting point (lower bound of the domain)
B=1                                 #upper bound of space domain
C=1                                 #upper bound of time domain
N=100                               #space nodes in each direction
h=(B-A)/(N-1)                       #grid size
k=0.001                             #time step
k_check=h/np.sqrt(2)                #time step to check for convergence

if k>k_check:                       #check the convergence 
    print('the mesh size is too small, causing divergence in the solution. The progam will be aborted')    #display an error message
    exit()                          #quit the program if convergnece can't be achieved

x=np.arange(A,B+h,h)                #x direction grid
y=np.arange(A,B+h,h)                #y direction grid
X, Y = np.meshgrid(x, y)            #2D space mesh grid
t=np.arange(A,C+k,k)                #time interval grid
Nt=len(t)                           #number of time nodes


# In this section I am defining arrays I would need (if neeeded)
u=np.zeros((len(t),len(x),len(y)))         #setting my array u, the wave's amplitude
fx=np.zeros((len(x),len(y)))               #setting an array to assign the sarting condition for the grid


# In this section I am setting the boundary conditions/initial values
def f(x,y):
    return (0.5*np.exp(-700*((x-0.5)**2+(y-0.15)**2)))      #function to describe the initial condition

for i in range(len(x)):
    for j in range(len(y)):
        u[0,i,j] = f(x[i],y[j])                                   #each grid point gets assigned a value based on the starting function
        
gx=0                                                              



# In this section I am implementing the numerical method

for p in range(1,np.shape(u)[0]):                                     #iterating through time. 
                                                                      #the initial condition is set already (hence start from p=1)
    for q in range(1,np.shape(u)[1]-1):                               #iterating through one space dimension 
                                                                      #the BC is set (hence start from q=1 and finish at np.shape(u)[1]-1)
        for n in range(1,np.shape(u)[2]-1):                           #iterating through the other space dimension 
                                                                      #the BC is set (hence start from q=1 and finish at np.shape(u)[2]-1)
            if p==1:                                                                     #implementing the numerical method for the first time step
                u[p,q,n]=0.25*(u[0,q+1,n]+u[0,q-1,n]+u[0,q,n+1]+u[0,q,n-1])+k*2*gx
            if p!=1:                                                                     #implementing the numerical method for the follwoing time period
                u[p,q,n]=0.5*(u[p-1,q-1,n]+u[p-1,q+1,n]+u[p-1,q,n-1]+u[p-1,q,n+1])-u[p-2,q,n]

                        
end=time.time()                                                        #ending the time count for the numerical method to be computed
print('Time elapsed to compute the numerical method is ',round(end-start,2),'s')  #displaying the elapsed time 


# In this section I am showing the results
print('The mesh size is h =', round(h,5))                           #print info on mesh size
print('The time step is k =', round(k,5))                           #print info on time step
print('The max time step is k_max =', round(k_check,5))             #print info on time step convergent

contour_size_inches = 6                                             #size to match contour gif dimensions
#sliced plane development plot
pts=[0,60,100,180,250]                                             #some instants are chosen to be plotted
for item in pts:
    fig = pl.figure(figsize=(contour_size_inches, contour_size_inches) if item == 60 else pl.rcParams['figure.figsize'])  #match gif size for t=60 only
    pl.plot(y,u[item,int(0.5/h),:])                                #plot the surface sliced in the mid plane
    pl.grid(which='both')                                          #plot the grid
    if item != 60:
        pl.title('Acoustic wave profile along midplane \
parallel to x')                                                    #plot title (skip for t=60)
    pl.xlabel('y')                                                 #label the x axis
    pl.ylabel('Amplitude')                                         #label the y axis (renamed from 'u')
    pl.xlim(0, 1)                                                  #set x axis range from 0 to 1
    ax = pl.gca()
    ax.yaxis.set_ticklabels([])                                    #hide vertical axis tick values
    pl.savefig(f"output/midplane_t{item}.png")                     #save the plot to output folder
    pl.close()                                                     #close the plot



#create a gif to see the evolution in 2D
frames = 200                                                   #set the number of frames
nt=0                                                           #start the gif at time 0
for n in range(frames):                                        #Generate each frame
    fig = pl.figure(figsize=(contour_size_inches, contour_size_inches))   #match size of contour gif
    ax = fig.add_subplot(111, projection='3d')                 #create subplot
    ax.plot_surface(X, Y, u[nt], cmap='Blues',                 
                       linewidth=0, antialiased=False)         #plot the surface
    ax.set_zlim(-1,1)                                          #set the axis limit of the function axis
    ax.zaxis.set_major_locator(LinearLocator(10))              #formatting of the z axis
    ax.zaxis.set_major_formatter('{x:.02f}')                   #formatting of the z axis
    ax.set_xlabel('x')                                         #label the x axis
    ax.set_ylabel('y')                                         #label the y axis
    ax.set_zlabel('Amplitude')                                 #label the z axis (renamed from 'u')
    ax.set_xticklabels([])                                     #hide x axis tick values
    ax.set_yticklabels([])                                     #hide y axis tick values
    ax.set_zticklabels([])                                     #hide z (amplitude) axis tick values
    pl.savefig(f"output/{n}.png")                              #save each frame to output folder
    pl.close()                                                 #close the plot
    
    fig, ax2d = pl.subplots(figsize=(contour_size_inches, contour_size_inches))  #square figure matching 3D gif size
    ax2d.contourf(X, Y, u[nt,:,:], cmap='Blues')               #plot the contour
    ax2d.set_aspect('equal')                                    #make both axes the same dimension
    ax2d.set_xlabel('x')                                        #label the x axis
    ax2d.set_ylabel('y')                                        #label the y axis
    # no title for contour plot
    pl.savefig(f"output/c{n}.png")                             #save each frame to output folder
    pl.close()                                                  #close the plot
    
    nt+=1                                                      #increment the time for the frame

images = [Image.open(f"output/{n}.png") for n in range(frames)]       #Use pillow to iterate through all frames
images[0].save('output/wave.gif', save_all=True, append_images=images[1:], duration=60, loop=0)    #saving all frames as an animation in a gif file

    
images_c = [Image.open(f"output/c{n}.png") for n in range(frames)]       #Use pillow to iterate through all frames
images_c[0].save('output/contour.gif', save_all=True, append_images=images_c[1:], duration=60, loop=0)    #saving all frames as an animation in a gif file



# In this section I am celebrating
print('Done')