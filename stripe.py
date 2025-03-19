# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
import scipy.interpolate as si


#relevant functions from paper "Batista2020"
def pendulum (y, s, landa):
    theta, omega = y
    dydt = [omega, -landa*np.sin(theta)]
    return dydt

def integrand(theta):
    return 1./np.sqrt(np.cos(theta)-np.cos(alpha[j]))

# experimental data 
L_up = [9.3, 8.6, 6.7, 5, 6.2, 7.9, 8.7]
H_up = [1.3, 2.5, 3.0, 4.8, 2.8, 2.1, 1.3]
L_down = [10.5, 9.3, 8.2, 6.4, 7.5, 9.0, 10.0]
H_down = [0.6, 1.8, 3.6, 4.3, 3.4, 2.5, 0.9]
length_up = 9.8 # mm
length_down = 10.8 # mm

#start code
l0 = [length_up, length_down]

alpha = np.linspace(1., 99.5, 10)
alpha = np.asarray(alpha)
alpha = alpha*np.pi/180

figure, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)

for length in l0:
    
    curv = []
    x_fin = []
    y_s_max = []

    for j in range (len(alpha)):
            
        landa = (8./length**2)*(integrate.quad(integrand, 0, alpha[j])[0])**2
        
        omega0 = np.sqrt(2.0*landa*(1.0 - np.cos(alpha[j])))
        
        y0 = [0, omega0]
        
        s = np.linspace(0, length, 100)
        sol = odeint(pendulum, y0, s, args=(landa,))
        
        x_s = []
        y_s = []
        
        interpolation_x_prime = si.interp1d(s, np.cos(sol[:, 0]), kind = "slinear")
        
        for j in s:
            x_s.append(integrate.quad(interpolation_x_prime, 0, j)[0])
        
        
        interpolation_y_prime = si.interp1d(s, np.sin(sol[:, 0]), kind = "slinear")
        
        for j in s:
            y_s.append(integrate.quad(interpolation_y_prime, 0, j)[0])
        
        x_s = np.asarray(x_s)
        y_s = np.asarray(y_s)
        
        if (length == length_up):
            ax1.plot(x_s, y_s)
            with open("Piegature_up.txt", "ab") as f:
                np.savetxt(f, (x_s, y_s))
            ax1.set_xlabel('L [mm]', fontsize=20)
            ax1.set_ylabel('H [mm]', fontsize=20) 
            ax1.xaxis.set_tick_params(labelsize=15)
            ax1.yaxis.set_tick_params(labelsize=15)
        else:
            ax5.plot(x_s, y_s)
            with open("Piegature_down.txt", "ab") as f:
                np.savetxt(f, (x_s, y_s))
            ax5.set_xlabel('L [mm]', fontsize=20)
            ax5.set_ylabel('H [mm]', fontsize=20) 
            ax5.xaxis.set_tick_params(labelsize=15)
            ax5.yaxis.set_tick_params(labelsize=15)
            
       
        x_fin.append(x_s[-1])
        
        y_s_max.append(y_s[len(s)/2])
        
        df1 =np.gradient(y_s, x_s)
        
        df2 = np.gradient(df1, x_s)
        
        curv.append(np.abs(df2[len(s)/2])/(1.0 + df1[len(s)/2]**2)**1.5) #curvature
        
        
    x_fin = np.asarray(x_fin)
    y_s_max = np.asarray(y_s_max)
    curv = np.asarray(curv)
    
    radius = 1./curv
    
    if (length == length_up):
        
        ax2.plot(x_fin[1:], radius[1:], '-o', color = 'blue')
        np.savetxt("L-R_up.txt", (x_fin[1:], radius[1:]))
        interp_L_R = si.interp1d(x_fin[1:], radius[1:], kind = "slinear")
        ax2.set_xlabel('L [mm]', fontsize=20)
        ax2.set_ylabel('R [mm]', fontsize=20) 
        ax2.xaxis.set_tick_params(labelsize=15)
        ax2.yaxis.set_tick_params(labelsize=15)
        
        ax3.plot(y_s_max[1:], radius[1:], '-o', color = 'red')
        np.savetxt("Hmax-R_up.txt", (y_s_max[1:], radius[1:]))
        ax3.set_xlabel('Hmax [mm]', fontsize=20)
        ax3.set_ylabel('R [mm]', fontsize=20)  
        ax3.xaxis.set_tick_params(labelsize=15)
        ax3.yaxis.set_tick_params(labelsize=15)
        
        ax4.plot(x_fin, y_s_max, '-o', color = 'black', linewidth = 4, markersize = 10, label = 'theory')
        np.savetxt("L-Hmax_up.txt", (x_fin, y_s_max))
        ax4.set_xlabel('L [mm]', fontsize=20)
        ax4.set_ylabel('Hmax [mm]', fontsize=20)  
        ax4.xaxis.set_tick_params(labelsize=15)
        ax4.yaxis.set_tick_params(labelsize=15)
        
        ax4.plot(L_up, H_up, 'o', color = 'red', label = 'exp_up', markersize = 10)
        ax4.legend()
        
        R_up = interp_L_R(L_up)
        
        np.savetxt("L-H-R_up.txt", (L_up, H_up, R_up))
    
    else:
        
        ax6.plot(x_fin[1:], radius[1:], '-o', color = 'blue')
        np.savetxt("L-R_down.txt", (x_fin[1:], radius[1:]))
        interp_L_R = si.interp1d(x_fin[1:], radius[1:], kind = "slinear")
        ax6.set_xlabel('L [mm]', fontsize=20)
        ax6.set_ylabel('R [mm]', fontsize=20) 
        ax6.xaxis.set_tick_params(labelsize=15)
        ax6.yaxis.set_tick_params(labelsize=15)
        
        ax7.plot(y_s_max[1:], radius[1:], '-o', color = 'red')
        np.savetxt("Hmax-R_down.txt", (y_s_max[1:], radius[1:]))
        ax7.set_xlabel('Hmax [mm]', fontsize=20)
        ax7.set_ylabel('R [mm]', fontsize=20)  
        ax7.xaxis.set_tick_params(labelsize=15)
        ax7.yaxis.set_tick_params(labelsize=15)
        
        ax8.plot(x_fin, y_s_max, '-o', color = 'black', linewidth = 4, markersize = 10, label = 'theory')
        np.savetxt("L-Hmax_down.txt", (x_fin, y_s_max))
        ax8.set_xlabel('L [mm]', fontsize=20)
        ax8.set_ylabel('Hmax [mm]', fontsize=20)  
        ax8.xaxis.set_tick_params(labelsize=15)
        ax8.yaxis.set_tick_params(labelsize=15)
        
        ax8.plot(L_down, H_down, 'o', color = 'blue', label = 'exp_down', markersize = 10)
        ax8.legend()
        
        R_down = interp_L_R(L_down)
        
        np.savetxt("L-H-R_down.txt", (L_down, H_down, R_down))
        
