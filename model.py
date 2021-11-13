# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 13:43:49 2021

@author: Jean
"""

import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import streamlit as st

def dy(t,y,const=[1,1,1,1,1,1]):
    k1,k2,k3,k4,k5,k6=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*rd*a-k5*a*rt+k4*c+k6*c*rd
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k5*c*rd
    da= -k5*a*rt+k4*c
    
    return [drt,dc,drd,da]
    
y0=[1,1,1,1]

t_min=0
t_max=10
t=np.arange(t_min,t_max,0.0001)

const=[1,1,1,1,1,1]

r=solve_ivp(dy,[t_min,t_max],y0,args=[const],t_eval=t,dense_output=True)

y=r.sol(t)

for i in range(4):
    plt.subplot(2,2,i+1)
    plt.plot(t,y[i])
#st.pyplot(fig)

#%%

import scipy.optimize

def f(y,const=[1,1,1,1,1,1]):
    k1,k2,k3,k4,k5,k6=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*rd*a-k5*a*rt+k4*c+k6*c*rd
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k5*c*rt
    da= -k5*a*rt
    
    return [drt,dc,drd,da]

x = scipy.optimize.newton_krylov(f, [1,1,1,1],f_tol=1e-14)
print(x)
