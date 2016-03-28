# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 12:59:38 2016

@author: abigail_thayer
"""

#import numpy as np
import landlab as ll
import matplotlib.pyplot as plt
from landlab.plot.imshow import imshow_node_grid
from pylab import show, figure

"""initialize"""

#constants
k = 2. #w/mK thermal conductivity
rho = 2000. #kg/m3 density
Cp = 2000. #J/kgK heat capacity
kappa = k/(rho*Cp) 

#temperature
geotherm = 25. #degrees C/km
Ts = -5. #degrees C surface temp
Tbasal = 1500. #degrees C basal flux
qm = 0.04 #w/m2 heat flux

#space array
num_cols = 201 #x array
num_rows = 201 #z array
dx = 40. #why is this also dz?

#time array
nyears = 100; #run time in years
tmax = nyears*24*3600*365. #run time in seconds
#tmax = np.pi*(10**7)*200000. #total run time in seconds from Rachel?
dt = 0.2*(dx**2/kappa) #time step necessary to prevent instability

#create grid
mg = ll.RasterModelGrid(num_rows,num_cols,dx)

#fill data fields with initial zeros
z = mg.add_zeros('node', 'depth')
x = mg.add_zeros('node', 'distance')
T = mg.add_zeros('node', 'temperature')
q = mg.add_zeros('node', 'heat flux')
dTdz = mg.add_zeros('link', 'temperature_gradient')

#set boundary conditions (bottom, left, top, right)
mg.set_closed_boundaries_at_grid_edges(True, True, True, False) 

#initialize arrays
z[:]=mg.node_y
x[:]=mg.node_x
T[:] = Ts + z*qm/k

#this tells the computer that we have a square of magma at the bottom
for i in range (int(mg.number_of_nodes)):
    if 2000 < int(x[i]):
        if int(x[i]) < 6000:
            if int(z[i]) > 7000:
                T[i] = Tbasal
           #  T[z>7000] = Tbasal #can the upper lines be more efficient?
            
"""run"""
for j in range (int(tmax/dt)):
    dTdz[mg.active_links]=mg.calculate_gradients_at_active_links(T)
    Q=-k*dTdz
    dQdz= mg.calculate_flux_divergence_at_nodes(Q[mg.active_links])
    T[mg.core_nodes] += (-dQdz[mg.core_nodes])/(rho*Cp) * dt
   
figure('Temperature profile')
imshow_node_grid(mg, 'temperature')
figure('plot')
imshow_node_grid(mg, z)
show()
             
#plt.figure(1)
#im = plt.imshow(T, cmap=plt.cm.jet, extent=[0,400,0,400], origin='lower')  # display a colored image
#plt.colorbar(im)
#plt.title('Temperature Profile of Cooling Magma Body')

#plt.show()














