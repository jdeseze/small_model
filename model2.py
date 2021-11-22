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
    k1,k2,k3,k4,k5,k6,kon,koff,a0=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*rd*a-k5*a*rt+k4*c+k6*c*rd
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k6*c*rd
    #da= -k5*a*rt+k4*c
    da=kon*10**a0-koff*a-k5*a*rt+k4*c

    return [drt,dc,drd,da]

t_min=0
t_max=10
t=np.arange(t_min,10*t_max,0.01)

with st.sidebar:
    k1=st.slider('k1 (rd->rt)',0.0,2.0,0.5)
    k2=st.slider('k2 (rt->rd)',0.0,2.0,0.5)
    k3=st.slider('k3 (a+rd->a+rt)',0.0,2.0,0.5)
    k4=st.slider('k4 (c->rt+a)',0.0,2.0,0.5)
    k5=st.slider('k5 (a+rt->c)',0.0,2.0,0.5)
    k6=st.slider('k6 (c+rd->c+rt)',0.0,2.0,0.5)
    kon=st.slider('kon ',0.0,2.0,0.5)
    koff=st.slider('koff',0.0,2.0,0.5)
    a0=st.slider('[A]0^',-3.0,3.0,0.5)
    fact=st.slider('factor of increase',1,100,1)
const=[k1,k2,k3,k4,k5,k6,kon,koff,a0]
y0=[1,0,0,10**a0]

r=solve_ivp(dy,[t_min,10*t_max],y0,args=[const],t_eval=t,dense_output=True)

y=r.sol(t)

plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']

c=st.columns(2)
for i in range(4):
    with c[int(i/2)]:
        source=pd.DataFrame({'t':t,plot_name[i]:y[i]})
        chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
        #st.altair_chart(chart)

t1=np.arange(t_min,1,0.01)
r1=solve_ivp(dy,[t_min,1],[y[0][-1],y[1][-1],y[2][-1],y[3][-1]],t_eval=t1,args=[const],dense_output=True)
y1=r1.sol(t1)

t2=np.arange(1,t_max,0.01)
r2=solve_ivp(dy,[1,t_max],[y1[0][-1],y1[1][-1],y1[2][-1],fact*y1[3][-1]],args=[const],t_eval=t2,dense_output=True)
y2=r2.sol(t2)

plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']

c=st.columns(2)
for i in range(4):
    with c[int(i/2)]:
        source=pd.DataFrame({'t':list(t1)+list(t2),plot_name[i]:list(y1[i])+list(y2[i])})
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


def f(const):
    k1,k2,k3,k4,k5,k6,kon,koff,a0=const
    return  (k1*(1-a0*(1-(kon/koff))+k3*(1-a0+(kon*a0/koff))*(kon*a0/koff)+k6*(a0-(kon*a0/koff))*(1-a0+(kon*a0/koff)))/(k1+k2+k3*(kon*a0/koff)+k6*(a0-(kon*a0/koff)))

x=np.logspace(-2,2,100)
z=[f(const[0:-1]+[i]) for i in x]
with c[0]:
    source=pd.DataFrame({'ARHGEF11 concentration':x,'RhoGTPeq/RhoTOT':z})
    chart=alt.Chart(source).mark_line().encode(
        alt.X('ARHGEF11 concentration',scale=alt.Scale(type='log')),
        alt.Y('RhoGTPeq/RhoTOT',scale=alt.Scale(type='log'))).interactive()
    st.altair_chart(chart)
