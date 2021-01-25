# -*- coding: utf-8 -*-
"""
Created on Fri May 19 17:14:51 2017
@author: Srilok Srinivasan (srilok@iastate.edu)
"""

import numpy as np 

class Eigen: 
    """
    Eigen class
    Gets the following attributes from the eigen file from GULP
    
    Number of sample K points = n_sample_k 
    Number of basis atoms = n_basis 
    Number of branches = n_branches 
    Real part of the eigenvector = eigen_vector_Re 
    Imaginery part of the eigenvector = eigen_vector_Im 
    Total eigenvector = eigen_vector
    
    The eigen vectors are 4 dimensional arrays with the following axis 
    axis 1 = Index of the sample K point 
    axis 2 = Index of the phonon branch 
    axis 3 = Index of the basis atom 
    axis 4 = x,y,z components of the eigen vector 
    
   """
    
    def __init__(self,filename):
        check = filename.split()
        if len(check)!=1: 
           raise StandardError,"Invalid arguments for eigen" 
           
        try:
            with file(filename) as f: 
                s = f.read()
        except: 
            raise StandardError,"Cannot read the file " + str(filename)
        #Get lines in eigne file    
        lines = s.split('\n')
        self.words =[]
        for line in lines:
            self.words.append(line.split())
            
        #Number of basis atoms     
        self.n_basis = int(self.words[0][0])
        #Number of sampling K-points 
        self.n_sample_k = int(self.words[self.n_basis+1][0])
        
        #Number of phonon branches
        self.n_branches = int(self.words[self.n_basis+2][0])
        
        #Memory allocation 
        self.k_points = np.zeros([self.n_sample_k,3])
        self.eigen_vector_Re = np.zeros([self.n_sample_k,self.n_branches,self.n_basis,3])
        self.eigen_vector_Im = np.zeros([self.n_sample_k,self.n_branches,self.n_basis,3])
        self.eigen_value     = np.zeros([self.n_sample_k,self.n_branches])
        
        self.n_lines_per_mode  = (2 + self.n_basis)
        self.n_lines_per_block =  self.n_lines_per_mode * self.n_branches + 1
        
        
        skip_header = 3 + self.n_basis
        
        for i in range(self.n_sample_k):
            self.k_points[i] = np.array([float(k) for k in 
                          self.words[self.n_lines_per_block*i+skip_header][3:6]])
            for j in range(self.n_branches):
                index = self.n_lines_per_block*i + j*self.n_lines_per_mode + skip_header + 2
                self.eigen_value[i,j] = np.array([float(k) for k in self.words[index]]) 
                for l in range(self.n_basis):
                    if(len(self.words[index+l+1]))==3:
                        self.eigen_vector_Re[i,j,l] = np.array(
                        [float(k) for k in self.words[index+l+1]])
                    if(len(self.words[index+l+1]))==6: 
                        self.eigen_vector_Re[i,j,l] = np.array(
                        [float(k) for k in self.words[index+l+1][0:3]])
                        self.eigen_vector_Im[i,j,l] = np.array(
                        [float(k) for k in self.words[index+l+1][3:]])
                        
        self.eigen_vector = self.eigen_vector_Re[:] + 1j*self.eigen_vector_Im[:]
        self.eigen_value = self.eigen_value * 0.03 # cm^-1 to THz
    
    def get_eigenvalue(self): 
        return self.eigen_value
    
    def get_eigenvector(self): 
        return self.eigen_vector
    
    def get_num_branches(self):
        return self.n_branches
    
    def get_eigen_gamma_M(self):
        gamma_M = np.array([1,0,0])
        sample_kpoints =[]
        sample_eigenvalue=[]
        sample_eigenvector=[]
        
        for i in range(len(self.k_points)): 
            if np.all((self.k_points[i]*gamma_M)-self.k_points[i] == [0,0,0]): 
                sample_kpoints.append(self.k_points[i])
                sample_eigenvalue.append(self.eigen_value[i])
                sample_eigenvector.append(self.eigen_vector[i])
                
        sample_kpoints = np.array(sample_kpoints)
        sample_eigenvalue = np.array(sample_eigenvalue)
        sample_eigenvector = np.array(sample_eigenvector)
        
        return [sample_kpoints,sample_eigenvalue, sample_eigenvector]
    
    def phonon_vel_gamma_M(self):
        [k_points,eigenvalue,eigenvector] = self.get_eigen_gamma_M()
        del_k_x = k_points[1][0] - k_points[0][0]        
        phonon_vel = np.gradient(eigenvalue)/del_k_x
        return phonon_vel
            
     
                