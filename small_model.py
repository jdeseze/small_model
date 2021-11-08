# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 12:47:20 2021

@author: Jean
"""
import numpy as np
import matplotlib.pyplot as plt

def fp(d,r,k):
    return -(d-r+1/k)+((d-r+1/k)**2+4*r/k)**0.5
    
def fn(d,r,k):
    return -(d-r+1/k)-((d-r+1/k)**2+4*r/k)**0.5

d=np.arange(0,10,0.1)

plt.plot(d,[fp(x,1,0.5) for x in d])
