# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 15:32:11 2020

@author: Yung Han Jeong
"""
from sympy import *
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

#solves tandem trailor weight distribution. takes in geometric values of the trailer
#returns the solves forces A = Landing Gear B = Front Axle C= Rear Axle and Moment at the front axle
def solvetrailer(dist_weight, chasis_length, axle_spread):
    totalweight = dist_weight * chassis_legnth
    
    x = symbols('x') #initialize X this is needed for integration
    #w = weight per length
    #l = length of the beam
    #a = space between axles
    #c = constants generated during integration
    w, l, a, c, Ma = symbols('w l a c Ma') #initialize geometric constants
    A, B, C = symbols('A B C') #initialize reactionary forces and moments
    
    print('Solving Statically Indeterminate Beam using Double Integration Piece Wise Fx','\n')
    #-----------------------------FIRST SECTION------------------------------------
    #set Moment A equation from Sum of Moment at C
    ma = A*l + B*a - w*l**2/2
    mexp1 = ma + A*x-w*x**2/2
    #----------------------------------FIRST INTEGRATION---------------------------
    mexpint1 = integrate(mexp1,x)
    mexpint1 = mexpint1+mexpint1.subs(x,0) #+c, which is 0 in this case
    #---------------------------------SECOND INTEGRATION--------------------------
    mexpint2 = integrate(mexpint1,x)
    mexpint2 = mexpint2+mexpint2.subs(x, 0) #+c, which is 0 in this case
    #if c in the second integration is equal to itself
    
    #---------------SOLVING FOR A in respect to B and w ---------------------------
    mexpint2 = expand(mexpint2)
    eqA = simplify(mexpint2.subs([(A,0),(x,(l-a))])/(-1*mexpint2.subs([(x,(l-a)),(B,0),(w,0),(A,1)])))
    
    #-------------------------------THIRD SECTION----------------------------------
    mexp2 = w*x**2/2 - C*x #initialize equation
    #----------------------------------FIRST INTEGRATION---------------------------
    mexp2int1 = integrate(mexp2,x)
    mexp2int1 = mexp2int1 + mexp2int1.subs(x,l)
    #----------------------------------SECOND INTEGRATION--------------------------
    mexp2int2 = integrate(mexp2int1,x)
    mexp2int2 = mexp2int2 + mexp2int2.subs(x,l)
    mexp2int2=expand(mexp2int2)
    eqC = mexp2int2.subs([(C,0),(x,(l-a))])/(-1*mexp2int2.subs([(x,(l-a)),(w,0),(C,1)]))
    
    #print('0 = w*l + ', eqA, ' + B ', eqC,'\n')
    sumf = w*l + eqA + eqC + B
    sumf = expand(sumf)
    eqB = simplify(sumf.subs(B,0)/sumf.subs([(w,0),(B,1)]))
    solA = eqA.subs(B,eqB)
    #-----------------------------Substitute Values to solve---------------------------#
    
    reactionA = N(solA.subs([(w,dist_weight),(l,chasis_length),(a,axle_spread)]),8)
    reactionC = N(eqC.subs([(w,dist_weight),(l,chasis_length),(a,axle_spread)]),8)
    reactionB = totalweight - reactionA - reactionC
    reactionMa = ma.subs([(A,reactionA),(B,reactionB),(w,dist_weight),(l,chasis_length),(a,axle_spread)])
    
    print('REACTIONARY FORCES')
    print('A = ',reactionA,"lbs")
    print('B = ',reactionB, "lbs")
    print('C = ',reactionC,"lbs")
    print('Ma = ',reactionMa,"lbs*ft")
    return [reactionA, reactionB, reactionC, reactionMa]