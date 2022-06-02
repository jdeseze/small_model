# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 13:43:49 2021

@author: Jean
"""

import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import streamlit as st
import altair as alt
import pandas as pd

st.set_page_config(page_title="Model", page_icon=":model:",layout="wide")

t_min=0
t_max=10
t=np.arange(t_min,100*t_max,0.01)

with st.sidebar:
    kon=st.slider('kon',-5.0,5.0,0.0)
    koff=st.slider('koff',-5.0,5.0,0.0)
    kcat=st.slider('kcat',-5.0,5.0,0.0)
    km=st.slider('km',-5.0,10.0,0.0)
    k1=st.slider('k1',-5.0,5.0,0.0)
    k2=st.slider('k2',-5.0,5.0,0.0)
    minit=st.slider('minit',-5.0,5.0,0.0)
    rinit=st.slider('rinit',-5.0,5.0,0.0)
    
const=[10**kon,10**koff,10**kcat,10**km,10**k1,10**k2]



c=st.columns(2)

def dy(t,y,const):
    kon,koff,kcat,km,k1,k2=const
    r,rg,rt,m,c=y
    
    g=1
    
    dr=-kon*r*g+koff*rg
    drg=kon*g*r-koff*rg-kcat*rg
    drt=kcat*rg-km*rt
    dm=km*rt-k1*m*g+k2*c
    dc=k1*m*g-k2*c
    
    return  [dr,drg,drt,dm,dc]

t_min=0
t_max=100
    
t2=np.arange(t_min,t_max,0.01)


r=solve_ivp(dy,[t_min,10*t_max],[10**rinit,0,0,10**minit,0],args=[const],t_eval=t2,dense_output=True,method='BDF')
y=r.sol(t2)

plot_name=['Rho','Rho+G','RhoT','Myosin','G+Myosin']

c=st.columns(2)
for i in range(5):
    with c[int(i%2)]:
        source=pd.DataFrame({'t':list(t2),plot_name[i]:list(y[i])})
        chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
        st.altair_chart(chart)
  
# =============================================================================
# i=4 
# with c[0]:
#     source=pd.DataFrame({'t':list(t2),plot_name[i]:list(y2[i])})
#     chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
#     st.altair_chart(chart)
# =============================================================================
# =============================================================================
#         fig=plt.figure()
#         plt.subplot(2,2,i+1)
#         plt.plot(t,y[i])
#         st.pyplot(fig)
# =============================================================================
#st.pyplot(fig)