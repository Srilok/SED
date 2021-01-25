# -*- coding: utf-8 -*-
"""
Created on Sat May 14 18:50:17 2016

@author: Srilok Srinivasan(srilok@iastate.edu)
"""
"""
Created on Wed Apr  6 13:02:42 2016
Calculates the spectral energy density (SED) from the velocity data in dump
files from lammps. 

Restrictions: The atom IDs in the dump file must be sorted. 

Dependencies: "dump.py" from the Pizza tools for lammps must be placed in the 
               working folder 

Input: 1. dump file containing velocity, 2. mapping in http://www.plt.imshow/fo
Output: SED(k,omega)


To do: Generalize code for other dispersion along other symmetry lines 
       (currently does only along gamma-M)               
"""
#%%
from sys import argv
from scipy.fftpack import fft,fftfreq
from dump import dump 
from tauk import tauk
from eigen import Eigen
import numpy as np 
import matplotlib.pyplot as plt

#import matplotlib.cm as cm
from matplotlib.colors import LogNorm,Colormap,Normalize
#import os
#%%
#Funtion definitions

def get_eigen_index(omega,ffreq):
    
    freq_range = int(len(ffreq)/2)
    flag=[]
    for i in range(len(omega)): 
        for j in range(freq_range):
            
            if j < freq_range-1:  
                
                for k in range(len(omega[i])):               
                    if abs(ffreq[j]) <= omega[i][k] < abs(ffreq[j+1]):         
                        flag.append([i,j])
                        
            elif j == freq_range-1:
                
                for k in range(len(omega[i])):
                    
                    if abs(ffreq[j]) <= omega[i][k]:
                        flag.append([i,j])
                        
    return np.array(flag)
    
#%%    
def plot_SED(SED,ffreq,qx,file_name):
    
    freq_range = int(len(ffreq)/2)
    freq = ffreq[0:freq_range-1]   
    #PlottingSED
    fig = plt.figure()
    fig.suptitle('SED', fontsize=14, 
                 fontweight = 'bold')
    plt.imshow(SED.T[0:freq_range,:],extent=(np.amin(qx)*2,
               np.amax(qx)*2,np.amin(freq),np.amax(freq)),
               origin='lower',aspect='auto',norm = LogNorm())          
    plt.xlabel("$\kappa/\kappa_M$")
    plt.ylabel('$\omega$/2$\pi$ (THz)')    
    fig.savefig(file_name,dpi=600,facecolor='w',edgecolor='w',
                orientation = 'landscape',
                fotmat='pdf')
    fig.clf()
    plt.close()
    
#%%    
def C_v(omega,T): 
    
    omega = 2*np.pi*omega*1e+12
    h_bar = 6.582e-16 # Planck constant (eV.s) 
    kb = 8.617e-5 # Boltzman constant (eV.K^-1)
    
    x = (h_bar*omega)/(kb*T)
    
    return ((h_bar*omega)**2*(x/(x-1)**2))/(kb*(T**2))

#%%
def thermal_cond(qx,thickness,num_branch,phonon_vel,tau):
    #Calculating thermal conductivity
    ev_2_J = 1.602e-19  
    k_scaling = ev_2_J*1e+22/(4*np.pi*thickness)
    k_s = np.zeros(num_branch)
    k_accumulation = np.zeros([num_branch*len(qx),3],dtype=np.float32)
    dq = qx[1]-qx[0]
    k=0
    for i in range(num_branch):
        
        for j in range(len(qx)):
            
            k_accumulation[k,0] = omega[j,i]
            k_accumulation[k,1] = abs(2*np.pi*phonon_vel[0,j,i])*abs(tau[j,i])
            k_per_mode = (phonon_vel[0,j,i]**2) * abs(
                        tau[j,i])*C_v(omega[j,i],T) * qx[j][0] * dq[0] 
            k_accumulation[k,2] = k_per_mode
            k_s[i] += k_per_mode
            k = k+1
    
    k_s = k_s * k_scaling 
    k_accumulation[:,2] = k_accumulation[:,2]*k_scaling
    
    return [k_s, k_accumulation]
     

#End of function definitions
#%%

script, input_file1, input_file2, input_file3, output_file = argv

if len(argv)!= 5: 
    print "Invalid input" 
    raise StandardError,"Usage: python [script_name] [dump_file] [map_file] \
                         [Eigen output from GULP] [SED_output]"    

#a  = float(raw_input("Enter the value of lattice constant:"))
dt_MD = float(raw_input(
        "Enter the timestep used in the MD simulation (in ps):")) 

mass = raw_input(
        "Enter the mass of the basis atom in the order used in simulations:").split()
mass = [float(i) for i in mass]

T = float(raw_input("Enter the Temperature (in K):"))
thickness = float(raw_input("Enter the thickness of the structure (in Angstr):"))

#Read the dump file
print "Reading Snapshots......."
d=dump(input_file1)

#List of timesteps
time_steps=d.time()

#Calculating timestep to be used for fourier transform 
dt = (time_steps[1]-time_steps[0])*dt_MD*1e-12 #secs

#Total time of integration 
t_o = (time_steps[-1]-time_steps[0])*dt_MD*1e-12 #secs

#Read the mapping info
map_dtype = np.dtype([('nx', np.int32), ('ny', np.int32),('dummy',np.int32),
                      ('basis',np.int32),('atom_ID',np.int32)])
map_data  = np.loadtxt(input_file2, dtype = map_dtype, skiprows=2)

#Total number of atoms and basis 
num_atoms = np.max(map_data['atom_ID'])
num_basis = len(np.unique(map_data['basis']))

#Check if the number of basis matches the number of mass input
if num_basis != len(mass):
    raise StandardError,"Number of basis atoms do not match the masses!" 

#Determine the number of units cells along each direction 
Nx=len(np.unique(map_data['nx']))
Ny=len(np.unique(map_data['ny']))
NT=Nx*Ny
print "Number of atoms:",num_atoms

#Sclaing constant for SED 
amu_2_kg = 1.66054e-27
scaling_const = amu_2_kg/(4*np.pi*t_o*NT)

#%%
#Obtaining eigenvalue and eigenvector information
print "Obtaining Eigenvectors and eigenvalues..."

e=Eigen(input_file3)

[qx,omega,ev] = e.get_eigen_gamma_M()

num_branch = e.get_num_branches()

#%%
#Allocating memory for the system size
print "Allocating Memory....."
vx= np.zeros([num_atoms,len(time_steps)],dtype=np.float32)[:]
vy= np.zeros([num_atoms,len(time_steps)],dtype=np.float32)[:]
vz= np.zeros([num_atoms,len(time_steps)],dtype=np.float32)[:]

qx_basis_component = np.zeros([num_basis,len(time_steps)],
                               dtype = np.complex128)[:]
qy_basis_component = np.zeros([num_basis,len(time_steps)],
                               dtype = np.complex128)[:]
qz_basis_component = np.zeros([num_basis,len(time_steps)],
                               dtype = np.complex128)[:]
    
sigma_qx = np.zeros([len(qx),num_branch,len(time_steps)],
                     dtype = np.complex128)[:]
sigma_qy = np.zeros([len(qx),num_branch,len(time_steps)],
                     dtype = np.complex128)[:]
sigma_qz = np.zeros([len(qx),num_branch,len(time_steps)],
                     dtype = np.complex128)[:]
total_q = np.zeros([len(qx),num_branch,len(time_steps)],
                     dtype = np.complex128)[:]

SED_pol = np.zeros([len(qx),num_branch,len(time_steps)],dtype = np.float64)
SED_tot = np.zeros([len(qx),len(time_steps)],dtype = np.float64)

#%%
#Collecting velocity data 
print "Collecting velocity data...."
for i in range(len(time_steps)):   
    vx[:,i],vy[:,i],vz[:,i]=d.vecs(time_steps[i],"vx","vy","vz")[:]

#%%
#Obtaining frequencies for fourier transform 
ffreq  = fftfreq(len(vx[0,:]),dt)/1e+12
freq_range = int(len(ffreq)/2)
freq = ffreq[0:freq_range-1]

flag = []

#Fourier transform of the velocities
print "Performing fourier transform...." 

for i in range(len(qx)):

    #Zeroing for accumlation at each K-point
    qx_basis_component.fill(0)
    qy_basis_component.fill(0)
    qz_basis_component.fill(0)
    
    #Summing over unit cells
    for j in range(num_atoms):    
        
        basis_index = map_data['basis'][j]
            
        vkx = vx[j,:]*np.exp(2j*np.pi*qx[i][0]*map_data['nx'][j])
        qx_basis_component[basis_index] += vkx[:]
        
        vky = vy[j,:]*np.exp(2j*np.pi*qx[i][0]*map_data['nx'][j])
        qy_basis_component[basis_index] += vky[:] 
        
        vkz = vz[j,:]*np.exp(2j*np.pi*qx[i][0]*map_data['nx'][j])
        qz_basis_component[basis_index] += vkz[:]
    
    #Summing over basis atoms        
    for l in range(num_branch):        
        for k in range(num_basis): 
            
            sigma_qx[i,l,:] += np.sqrt(mass[k])*np.conj(ev[i,l,k,0])* \
                                qx_basis_component[k]
            sigma_qy[i,l,:] += np.sqrt(mass[k])*np.conj(ev[i,l,k,1])* \
                                qy_basis_component[k]
            sigma_qz[i,l,:] += np.sqrt(mass[k])*np.conj(ev[i,l,k,2])* \
                                qz_basis_component[k]
        
        #Summing over cartesian components
        total_q[i,l,:] = sigma_qx[i,l,:] + sigma_qy[i,l,:] + sigma_qz[i,l,:]
    
        #Performing the fourier transform        
        SED_pol[i,l,:] = np.absolute(fft(total_q[i,l,:]))**2


SED_pol = SED_pol * scaling_const


#Summing over polarization branches
for l in range(num_branch): 
    SED_tot[:,:] += SED_pol[:,l,:]

#%%
#Dumping memory
del vx 
del vy 
del vz
del sigma_qx
del sigma_qy
del sigma_qz
del total_q
del vkx 
del vky 
del vkz

#%%
#Plotting Total SED
print "Plotting SEDs...."
plot_SED(SED_tot,ffreq,qx,'SED.pdf')

#Plotting individual branches 
for l in range(num_branch):
    plot_SED(SED_pol[:,l,:],ffreq,qx,'SED_mode_'+str(l)+'.pdf')

#%%
#Saving SED data 
#Writing SEDs
print "Saving Spectral energy density data..."

h= open(output_file+'tot.data', 'w')
#output_array = [2*qx[:,0],SED_tot[:,0:freq_range]]
np.savetxt(h,SED_tot[:,0:freq_range])
h.close()

for l in range(num_branch):
    h = open(output_file+'_'+str(l)+'.data', 'w')
    #output_array = np.concatenate(SED_pol[:,l,0:freq_range],axis=1)
    np.savetxt(h,SED_pol[:,l,0:freq_range])
    h.close()


index_flag = get_eigen_index(omega,ffreq)

#%%
#Phonon velocity    
phonon_vel = e.phonon_vel_gamma_M()

h= open(output_file+'.phon_vel', 'w')
np.savetxt(h,phonon_vel[0],"# ZA TA LA ZO TO LO (A/ps)")
h.close()

h= open(output_file+'.omega', 'w')
np.savetxt(h,omega,comments="# Omega/2*pi from GULP (THz)")
h.close()

#Relaxation time calculation by fitting SED peaks to Lorentzian
print "Fitting Lorentzian...." 
[tau,fit_freq] = tauk(SED_pol,index_flag,ffreq)

#Reshaping tau and frequency for plotting 
tau = np.reshape(tau,np.shape(omega))
fit_freq = np.reshape(fit_freq,np.shape(omega))

#Saving tau data
h= open(output_file+'.tau', 'w')
output_array = np.concatenate((2*qx,tau),axis=1)
np.savetxt(h,output_array,comments="# Normalized k, Tau (ps)")
h.close()

#Saving fit freq data
h= open(output_file+'.fit_freq', 'w')
output_array = np.concatenate((2*qx,fit_freq),axis=1)
np.savetxt(h,output_array,comments="# Normalized k, Fit frequency (THz)" )
h.close()
#%%
#Calculating thermal conductivity
[k_s,k_accumulation] = thermal_cond(qx,thickness, num_branch, phonon_vel,tau)
k_tot = np.sum(k_s)
print "Thermal Condunctivity = " +str(k_tot) + " W/m.K"

#Saving k data
h= open(output_file+'.k', 'w')
h.write("Total thermal conductivity = " + str(k_tot) + " W/m.K")
h.write("\n")

for l in range(num_branch): 
    h.write("k from mode " + str(l) + "= " + str(k_s[l]) + " W/m.K")
    h.write("\n")

h.close()
#%%
k_accumulation_freq = k_accumulation[np.argsort(k_accumulation[:,0])]
k_accumulation_freq[:,2] = np.cumsum(k_accumulation[:,2],axis=0)

k_accumulation_mfp = k_accumulation[np.argsort(k_accumulation[:,1])]
k_accumulation_mfp[:,2] = np.cumsum(k_accumulation[:,2],axis=0)

output_array = np.zeros([num_branch*len(qx),4],dtype=np.float32)
output_array[:,0] = k_accumulation_freq[:,0]
output_array[:,1] = k_accumulation_freq[:,2]
output_array[:,2] = k_accumulation_mfp[:,1]/10
output_array[:,3] = k_accumulation_mfp[:,2]

h= open(output_file+'.k_accumulation', 'w')
np.savetxt(h,output_array,
           comments="# Freq (THz) K_accumulation (W/mk) MFP (nm) K_accumulation (W/mK)")
h.close()


         
        
        


