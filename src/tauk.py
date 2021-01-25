# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:21:57 2017

@author: meisu
"""
#%%
import numpy as np 
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lorentzian(x,gamma): 
    # p = [hwhm,intensity,peak_freq]
    y = I*(gamma**2)/(gamma**2+ (x - wc)**2)
    return y
#%%
#def residual(hwhm,y,x,peak_intensity,peak_freq):
#    err = y - lorentzian(x,[hwhm,peak_intensity,peak_freq])
#    print x
#    return err
#%%
def plot_fig(i,freq_list,sample_data,fit):
    
    fig = plt.figure()
    plt.semilogy(freq_list,sample_data,'wo')
    plt.semilogy(freq_list,fit,'r-',lw=2)
    plt.xlabel(r'$\omega$ (THz$^{-1}$)', fontsize=18)
    plt.ylabel('Intensity (a.u.)', fontsize=18)              
    fig.savefig("FitFig/Fit"+str(i),dpi=100,facecolor='w',edgecolor='w',
                orientation = 'landscape',
                fotmat='tiff')
    fig.clf()
    plt.close()

##%% 
#def get_SED_peaks(SED,flag,ffreq): 
#    
#    for i in range(len(SED)): 
        
        
    

def tauk(SED_pol,flag,ffreq):
    
    
    
    #Checking if fitting figures directory exists 
    if not os.path.exists("FitFig/"):
        os.chmod("./",0777)
        os.makedirs("FitFig/",0777)
    
    freq_range = int(len(ffreq)/2)
    freq = ffreq[0:freq_range]
    freq = np.pad(freq,(0,len(freq)-1),'linear_ramp', end_values = (0,2*freq[-1]))

    best_fit = np.zeros([len(flag),3])
    
    fit_freq = np.zeros(len(flag))

    for i in range(len(flag)):
        
        q_index = flag[i][0]
        peak_index    = flag[i][1]

        j = i % 6
        
        if peak_index >= 50: 
            sample_data = SED_pol[q_index ,j,peak_index -50: peak_index+51]
            freq_list = freq[peak_index-50 : peak_index+51]
        else: 
            sample_data = SED_pol[q_index , j, 0:peak_index +101]
            freq_list = freq[0:peak_index + 101]           

        #Initial guess
        global I, wc
      
        I = np.max(sample_data)
        wc = freq_list[np.argmax(sample_data)]
        fit_freq[i] = wc
        gamma = 0.1
        
#        print I, wc
        #curve fitting 
        best_vals,covar = curve_fit(lorentzian,freq_list,sample_data,p0=[gamma],
                                    maxfev = 10000)
        best_fit[i,:] = np.array(best_vals)[:]
        
        fit = lorentzian(freq_list,best_vals[0])#,best_vals[1],best_vals[2])
        
        plot_fig(i,freq_list,sample_data,fit)
        
        
    tau = 1/((best_fit[:,0]*2)*2*np.pi)
    return tau,fit_freq

        