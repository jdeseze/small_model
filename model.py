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

st.set_page_config(page_title="Segmentation", page_icon=":model:",layout="wide")

def dy(t,y,const):
    k1,k2,k3,k4,k5,k6=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*rd*a-k5*a*rt+k4*c+k6*c*rd
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k6*c*rd
    da= -k5*a*rt+k4*c
    
    return [drt,dc,drd,da]
    


t_min=0
t_max=10
t=np.arange(t_min,t_max,0.01)

with st.sidebar:
    k1=st.slider('k1 (rd->rt)',0.0,1.0,0.5)
    k2=st.slider('k2 (rt->rd)',0.0,1.0,0.5)
    k3=st.slider('k3 (a+rd->a+rt)',0.0,1.0,0.5)
    k4=st.slider('k4 (c->rt+a)',0.0,1.0,0.5)
    k5=st.slider('k5 (a+rt->c)',0.0,1.0,0.5)
    k6=st.slider('k6 (c+rd->c+rt)',0.0,1.0,0.5)
    a0=st.slider('[A]0^',-3.0,3.0,0.5)
const=[k1,k2,k3,k4,k5,k6]
y0=[1,0,5,10**a0]

r=solve_ivp(dy,[t_min,t_max],y0,args=[const],t_eval=t,dense_output=True)

y=r.sol(t)

plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']

c=st.columns(2)
for i in range(4):
    with c[int(i/2)]:
        source=pd.DataFrame({'t':t,plot_name[i]:y[i]})
        chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
        st.altair_chart(chart)
        
r=solve_ivp(dy,[t_min,t_max],[y[0][-1],y[1][-1],y[2][-1],5*y[3][-1]],args=[const],t_eval=t,dense_output=True)

y=r.sol(t)

plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']

c=st.columns(2)
for i in range(4):
    with c[int(i/2)]:
        source=pd.DataFrame({'t':t,plot_name[i]:y[i]})
        chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
        st.altair_chart(chart)
# =============================================================================
#         fig=plt.figure()
#         plt.subplot(2,2,i+1)
#         plt.plot(t,y[i])
#         st.pyplot(fig)
# =============================================================================
#st.pyplot(fig)

#%%

import scipy.optimize

def f(y,const=const):
    k1,k2,k3,k4,k5,k6=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*rd*a-k5*a*rt+k4*c+k6*c*rd
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k6*c*rd
    da= -k5*a*rt+k4*c
    
    return [drt,dc,drd,da]



x=np.logspace(-2,2,100)
y = [scipy.optimize.newton_krylov(f, [y0[0],y0[1],y0[2],i],f_tol=1e-10) for i in x]
z=[y[i][0]/(y[i][2]+y[i][0]+y[i][1]) for i in range(len(y))]
fig=plt.figure()
ax = fig.add_subplot(2, 1, 1)
ax.plot(x,y)
ax.set_xscale('log')
with c[0]:
    source=pd.DataFrame({'ARHGEF11 concentration':x,'RhoGTPeq/RhoTOT':z})
    chart=alt.Chart(source).mark_line().encode(
        alt.X('ARHGEF11 concentration',scale=alt.Scale(type='log')),
        alt.Y('RhoGTPeq/RhoTOT',scale=alt.Scale(type='log'))).interactive()
    st.altair_chart(chart)
