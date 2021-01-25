# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 18:34:48 2017

@author: Srilok Srinivasan 
"""
from sys import argv 
from scipy.fftpack import fft,fftfreq
import numpy as np
import matplotlib.pyplot as plt 

def padwithzeros(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector
    
def smoothfft(vector,time_series, t_damping): 
    vector = fft(vector*np.exp(-time_series/t_damping))
    return vector
    
script, num_files, input_files, output_file = argv 

if len(argv) != 4: 
   print "Invalid input" 
   raise StandardError, "Usage: python [script_name] [number of files] [input_files] [output_file]" 
       
vacf_dtype = np.dtype([('t', np.float32), ('VACFx', np.float32),('VACFy',np.float32),
                      ('VACFz',np.float32),('VACFt',np.float32)])
input_file_str=str(input_files)

#Estimating the size of the input data 
temp = np.loadtxt(input_file_str+"-1.dat",vacf_dtype).shape

#Allocating memory
vacf_data = np.zeros((int(num_files),temp[0]),dtype=vacf_dtype)
vacf = np.zeros(temp[0],dtype = vacf_dtype)

for i in range(int(num_files)):
    vacf_data[i] = np.loadtxt(input_file_str+"-"+str(i+1)+".dat",vacf_dtype)

for i in range(len(vacf_data)): 
    vacf['t']        += vacf_data[i]['t']/len(vacf_data)
    vacf['VACFx']    += vacf_data[i]['VACFx']/len(vacf_data)
    vacf['VACFy']    += vacf_data[i]['VACFy']/len(vacf_data)
    vacf['VACFz']    += vacf_data[i]['VACFz']/len(vacf_data)
    vacf['VACFt']    += vacf_data[i]['VACFt']/len(vacf_data)

dt = (vacf['t'][1] - vacf['t'][0])*1e-12 #sec

#Plot Vacf-autocorrelatoion
plt.plot(vacf['t'],vacf['VACFx'])
plt.plot(vacf['t'],vacf['VACFy'])
plt.plot(vacf['t'],vacf['VACFz'])
#plt.plot(vacf['t'],vacf_autocorr_t)


#Zero padding 
vacf=np.pad(vacf[:],(0,len(vacf)),padwithzeros)

freq_range = fftfreq(len(vacf['VACFx']),dt) # Hz
freq_range = freq_range/1e+12 
pw_len = len(freq_range)

pw_x = np.zeros(pw_len,dtype = np.complex64)
pw_y = np.zeros(pw_len,dtype = np.complex64)
pw_z = np.zeros(pw_len,dtype = np.complex64)
pw_tot = np.zeros(pw_len,dtype = np.complex64 )

#Performing Fourier transform 
  #Smoothing function 
flag = str(raw_input("Do you want to use a smoothing function (y/n):")) 
if flag == 'y':
    t_damping = float(raw_input("Damping constant (in ps):"))
    pw_x = smoothfft(vacf['VACFx'],vacf['t'],t_damping)
    pw_y = smoothfft(vacf['VACFy'],vacf['t'],t_damping)
    pw_z = smoothfft(vacf['VACFz'],vacf['t'],t_damping)
    pw_tot = smoothfft(vacf['VACFt'],vacf['t'],t_damping)    
else:
    pw_x = fft(vacf['VACFx'])
    pw_y = fft(vacf['VACFy'])
    pw_z = fft(vacf['VACFz'])
    pw_tot = fft(vacf['VACFt'])
    
#Talking only the Real part
pw_x_real = np.real(pw_x)
pw_y_real = np.real(pw_y)
pw_z_real = np.real(pw_z)
pw_tot_real = np.real(pw_tot)

index= np.argmax(freq_range<0) 

#Normalizing phonon density of states
pw_x_real[0:index] = pw_x_real[0:index]/np.sum(pw_tot_real[0:index])
pw_y_real[0:index] = pw_y_real[0:index]/np.sum(pw_tot_real[0:index])
pw_z_real[0:index] = pw_z_real[0:index]/np.sum(pw_tot_real[0:index])
pw_tot_real[0:index] = pw_tot_real[0:index]/np.sum(pw_tot_real[0:index])

#Calculating in-plane DOS (x+y)
pw_in_plane = (pw_x_real[0:index] + pw_y_real[0:index])/2



g= open(output_file, 'w')
np.savetxt(output_file, [freq_range[0:index], pw_x_real[0:index], pw_y_real[0:index],pw_z_real[0:index],pw_tot_real[0:index]])

fmax = float(raw_input("Enter the maximum frequency in the DOS plot (in THz):"))  

fig = plt.figure() 
fig.suptitle('DOS Total', fontsize=14,fontweight = 'bold')
plt.plot(freq_range[0:index],pw_tot_real[0:index])
plt.xlim(0,fmax)
fig.savefig('DOS_total.pdf',dpi=600,facecolor='w',edgecolor='w',
            orientation = 'landscape', fotmat='pdf')
            
fig1 = plt.figure() 
fig1.suptitle('DOS in-plane', fontsize=14,fontweight = 'bold')
plt.plot(freq_range[0:index],pw_in_plane[0:index])
plt.xlim(0,fmax)
fig1.savefig('DOS_in_plane.pdf',dpi=600,facecolor='w',edgecolor='w',
            orientation = 'landscape', fotmat='pdf')
                        
fig2 = plt.figure() 
fig2.suptitle('DOS out-of-plane', fontsize=14,fontweight = 'bold')
plt.plot(freq_range[0:index],pw_z_real[0:index])
plt.xlim(0,fmax)
fig2.savefig('DOS_out_of_plane.pdf',dpi=600,facecolor='w',edgecolor='w',
            orientation = 'landscape', fotmat='pdf')
                        
##############################################################################
                        #Specific Heat Calculation 
choice = str(raw_input("Do you want to calculate the Heat capacity of the material(y/n):"))

if choice == 'y': 
    num_atoms = int(raw_input("Enter the number of atoms in the simulation:"))
    temp = float(raw_input("Enter the temperature (in K) to be used for C_v calculation:")) 
    #C_V = 3N*K_b*Integral(P(omega)*x^2*(exp(x)/(exp(x)-1)^2)); 
    # where x = h*omegha/(K_b*T)
    k_b = 8.671e-5 #Boltzman constant in eV.K^-1
    h = 4.135e-15 # eV.s 
    x= (h*freq_range[1:index]*1e+12)/(k_b*temp)
    d_omega = (freq_range[1] - freq_range[0])*1e+12
    integrand = pw_tot_real[1:index]*(x**2)*(np.exp(x)/(np.exp(x)-1)**2)
    C_v = 3*num_atoms*k_b*np.sum(integrand) #eV.K^-1
    print "Heat capacity at " + str(temp) + "K =" + str(C_v) +" eV/K"
    f = open('heat_capacity.log','w')
    f.write("Number of atoms: " + str(num_atoms) + "\n")
    f.write("Temperature: " +str(temp)+ "K \n")
    f.write("Heat Capacity: " + str(C_v) + "eV/K \n")
    f.close()

 
