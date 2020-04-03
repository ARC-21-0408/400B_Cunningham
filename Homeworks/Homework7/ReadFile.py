#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import astropy.units as u


def Read(filename): 
    #initialize
    data=[]
    time=0
    numpart=0
    
    file=open(filename, 'r') #open file
    line1=file.readline()
    label, timeval = line1.split() #get time value
    time=float(timeval)*u.Myr #give it units
    
    line2=file.readline() 
    label, numpart= line2.split() #get number of particles
    
    file.close()
    
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3) #get rest of data and put it in an array
    #print(time, numpart,data)
    
    return time,numpart,data #return the three values

