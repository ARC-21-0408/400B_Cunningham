#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Homework 4
# Center of Mass Position and Velocity
# Seen Cunningham


# In[1]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[2]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        
        # Function to compute the center of mass position or velocity generically
        # input: array (a,b,c) of positions or velocities and the mass
        # returns: 3 floats  (the center of mass coordinates)
        xsum=[]
        ysum=[]
        zsum=[]
        msum=np.sum(m)

        for i in range(len(a)):
            xsum.append(np.sum(a[i]*m[i]))
            ysum.append(np.sum(b[i]*m[i]))
            zsum.append(np.sum(c[i]*m[i]))
        
        xnum=np.sum(xsum)
        ynum=np.sum(ysum)
        znum=np.sum(zsum)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = xnum/msum
        # ycomponent Center of mass
        Bcom = ynum/msum
        # zcomponent Center of mass
        Ccom = znum/msum
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM =np.sqrt(XCOM**2+YCOM**2+ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = (self.x - XCOM)
        yNew = (self.y-YCOM)
        zNew = (self.z-ZCOM)
        RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2+YCOM2**2+ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM    

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("RMAX", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = (self.x-XCOM2)
            yNew = (self.y-YCOM2)
            zNew = (self.z-ZCOM2)
            RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        return np.round(COMP*u.kpc,2)
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x*u.kpc - COMX
        yV = self.y*u.kpc - COMY
        zV = self.z*u.kpc - COMZ
        RV = np.sqrt(xV**2+yV**2+zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV <= RVMAX)

        # determine the velocity and mass of those particles within the max radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        VXCOM2=np.round(VXCOM,2)
        VYCOM2=np.round(VYCOM,2)
        VZCOM2=np.round(VZCOM,2)
        
        
        # write your own code below
        COMV = [VXCOM2, VYCOM2, VZCOM2]

        # return the COM vector                                                                                        
        return COMV*(u.km/u.s)


# In[3]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt",2)
M33COM = CenterOfMass("M33_000.txt",2)


# In[4]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
#0.1kpc tolerance

MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])

#M31 and M33 position and velocity of COM 
M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])

#Print MW, M31, M33 COM position and velocity components
print("MW")
print(MW_COMP)
print(MW_COMV)
print("M31")
print(M31_COMP)
print(M31_COMV)
print("M33")
print(M33_COMP)
print(M33_COMV)


# In[ ]:


# now write your own code to answer questions for section 6


# In[5]:


#Current separation and velocity of MW and M31 (Use Distance Formula) [sqrt(x1-x2)^2+(y1-y2)^2+(z1-z2)^2]
#initialize
MW_M31PDiff=[]
MW_M31Dist=0
MW_M31VDiff=[]
MW_M31Vel=0

#subtract and square the componets of positions and velocity vecetors of MW and M31
MW_M31PDiff = np.subtract(MW_COMP,M31_COMP)**2
MW_M31Dist=np.sqrt(np.sum(MW_M31PDiff))

#velocity
MW_M31VDiff = np.subtract(MW_COMV,M31_COMV)**2
MW_M31Vel=np.sqrt(np.sum(MW_M31VDiff))

#print distance and veleocity between MW and M31
#print(MW_M31PDiff)
print(np.round(MW_M31Dist,3))#ANSWER FROM LECTURE 2 SLIDE 30: Distance ~ 770 Kpc
#print(MW_M31VDiff)
print(np.round(MW_M31Vel,3)) #ANSWER FROM LECTURE 2 SLIDE 33: Velocity ~ 110 km/s


# In[6]:


#current separation and velocity between M33 and M31
#initialize
M33_M31PDiff=[]
M33_M31Dist=0
M33_M31VDiff=[]
M33_M31Vel=0

#subtract and square the components of positions and velocity vectors of M33 and M31
M33_M31PDiff = np.subtract(M33_COMP,M31_COMP)**2
M33_M31Dist=np.sqrt(np.sum(M33_M31PDiff))

#Velocity
M33_M31VDiff = np.subtract(M33_COMV,M31_COMV)**2
M33_M31Vel=np.sqrt(np.sum(M33_M31VDiff))

#Print the distance and veleocity between M33 and M31
print(np.round(M33_M31Dist,3))
print(np.round(M33_M31Vel,3)) #ANSWER FROM LECTURE 2 SLIDE 33: Velocity ~ 202 km/s


# # Question 6.4
# 
# It is important to iterate this process becuase the COM's will constantly be moving and eventually merge. This will allow us to know when and where the MW and M31 will merge. As MW and M31 merge and begin to share mass, their COM will shift and our estimations will change. In order for the simulation to be accurate we need to know the COM of both galaxies until they merge. It won't be accurate if we assume a COM and then never update it.
# 
