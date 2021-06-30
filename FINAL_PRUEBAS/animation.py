import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
from collections import deque
import imageio

#Inputs
#nparticles = 10

fig, ax = plt.subplots()
line, = plt.plot([], [], 'ro')
#xdata = np.zeros((nparticles, N))
#xdata = [[] for i in range (nparticles)]
#ydata = [[] for i in range (nparticles)]

xdatafile = open ("xanimation_data.txt", "r")
ydatafile = open ("yanimation_data.txt", "r")
xdatafile.seek(0)
ydatafile.seek(0)
xlineaux = xdatafile.readline()
xlineaux = xlineaux.split('\t')
xlineaux[-1] = xlineaux[-1].strip()
nparticles = len(xlineaux)-1
print (nparticles)
xdatafile.seek(0)
ydatafile.seek(0)
xlineaux = xdatafile.readlines()
ndata = len(xlineaux)
print (ndata)
xdatafile.seek(0)
ydatafile.seek(0)

#xdata = [[0.0 for j in range (ndata)] for i in range (nparticles)]
#ydata = [[0.0 for j in range (ndata)] for i in range (nparticles)]

xdata = np.zeros((nparticles, ndata))
ydata = np.zeros((nparticles, ndata))


for ii in range (ndata):
        xline = xdatafile.readline()
        yline = ydatafile.readline()
        xline = xline.split('\t')
        yline = yline.split('\t')
        xline[-1] = xline[-1].strip()
        yline[-1] = yline[-1].strip()
        for jj in range (nparticles):
                xdata[jj][ii] = float(xline[jj])
                ydata[jj][ii] = float(yline[jj])
xdatafile.close()
ydatafile.close()


history_x, history_y = deque(maxlen=ndata), deque(maxlen=ndata)
def init(): #Inicializa, poniendo los limites de la grafica
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    return line,

thisx = np.zeros((nparticles, nparticles, ndata))
thisy = np.zeros((nparticles, nparticles, ndata))

def update(frame):
    for ii in range (nparticles):
        thisx[ii] = xdata[ii][frame]
        thisy[ii] = ydata[ii][frame]
        #thisy = [ydata[0][frame], ydata[1][frame], ydata[2][frame], ydata[3][frame]]

    if frame == 0:
        history_x.clear()
        history_y.clear()    
        
    history_x.appendleft(thisx[0])
    history_y.appendleft(thisy[0])
    line.set_data(thisx, thisy)    
    return line,
    
anim = FuncAnimation(fig, update, frames=ndata,
                    init_func=init, interval=300,blit=True) #interval = 20
#plt.show()

Writer = writers['ffmpeg']
#Writer = writers['imagemagick']
writer = Writer (fps=250, metadata={'artist':'Grupo3'},bitrate=1800) #bitrate = 1800
anim.save('Prueba.mp4', writer) #Prueba.mp4

