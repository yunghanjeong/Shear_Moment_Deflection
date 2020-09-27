# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 15:32:11 2020

@author: Yung Han Jeong
"""
import numpy as np
import math as mt

#--------------------------------FUNCTIONS------------------------------------
def __init__(self):
# at init, basic geometric constants are caluclated. These are GENERALLY a fixed number    
    OD = 5.0
    wallthick  = 0.25
    ID = OD - (2*wallthick)
    flexmod = 29000000 # Modulus of Elasticity
    inrt = (3.1415*(pow(OD,4)-pow(ID,4)))/64 #Moment of Inertia
    self.geo_const = flexmod * inrt # AKA E*I
    self.axl = 48 #Axle length
    
#statically indeterminate beam deflection calculation
def simp_im_def (self,x,a1,b1,a2,b2,w1,w2):
    imconst1 = w1*b1/(6*self.geo_const*self.axl)
    imconst2 = w2*b2/(6*self.geo_const*self.axl)
    ans1=0
    ans2=0
    if x <= a1:
        ans1= -imconst1*b1*x*(pow(self.axl,2)-pow(x,2)-pow(b1,2))
    if x > a1:
        ans1= -imconst1*b1*((self.axl/b1*pow((x-a1),3))+(x*(pow(self.axl,2)-pow(b1,2)))-pow(x,3))
    if x <= a2:
        ans2= -imconst2*b2*x*(pow(self.axl,2)-pow(x,2)-pow(b2,2))
    if x > a2:
        ans2= -imconst2*b2*((self.axl/b2*pow((x-a2),3))+(x*(pow(self.axl,2)-pow(b2,2)))-pow(x,3))
    return ans1 + ans2

#distance formula
def dist_form(x1,x2):
    return (mt.sqrt(pow(x2[0]-x1[0],2)-pow(x2[1]-x1[1],2)))

#strain
def strain(L1,L2):
    return abs(L1-L2)/L1