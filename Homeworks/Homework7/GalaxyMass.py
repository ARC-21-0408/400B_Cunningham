#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import astropy.units as u
import import_ipynb
from ReadFile import Read #import Read function from ReadFile

#define function to determine the component mass of a galaxy
def ComponentMass(filename,parttype):
    #get values from read file
    #a:time
    #b:total number of particles
    #c:data(m,x,y,z,vx,vy,vz)
    a,b,c=Read(filename)
    
    #create index for particle type
    index=np.where(c['type']==parttype)
    
    #get all the mass values for particle type
    # 10^10Mo * 10^12Mo / 10^2 10^10Mo
    #the m values are already in units of 10^10 so we need to multiply by 10^2 to get it in units of 10^12
    m=c['m'][index]*10**-2*u.solMass
    
    #sum all the mass components to get total mass
    mtot=np.sum(m)
    
    return mtot #return the total mass

    


