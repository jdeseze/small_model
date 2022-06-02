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
    k1=st.slider('k1 (rd->rt)',0.0,2.0,0.5)
    k2=50*k1#st.slider('k2 (rt->rd)',0.0,2.0,0.5)
    k3=k1#st.slider('k3 (a+rd->a+rt)',0.0,2.0,0.5)
    k4=st.slider('k4 (c->rt+a)',0.0,2.0,0.5)
    k5=1/k4#st.slider('k5 (a+rt->c)',0.0,2.0,0.5)
    k6=k3/5#st.slider('k6 (c+rd->c+rt)',0.0,2.0,0.5)
    kon=0.1#st.slider('kon ',0.0,2.0,0.5)
    koff=1/10#st.slider('koff',0.0,2.0,0.5)
    a0exp=st.slider('[A]0^',-3.0,3.0,0.5)
    fact=st.slider('factor of increase',0.0,2.0,1.0)
a0=10**a0exp
const=[k1,k2,k3,k4,k5,k6,kon,koff,10**a0]
y0=[1,0,0,10**a0]

c=st.columns(2)

# =============================================================================
# def dy(t,y,const):
#     k1,k2,k3,k4,k5,k6,kon,koff,a0=const
#     rt,c,rd,a=y
#     
#     drt= k1*rd-k2*rt+k3*(rd*a)-k5*a*rt+k4*c+k6*c*rd
#     dc= k5*a*rt-k4*c
#     drd= -k1*rd+k2*rt-k3*a*rd-k6*c*rd
#     #da= -k5*a*rt+k4*c
#     da=kon*(10**a0)-koff*a-k5*a*rt+k4*c
#     
#     
#     return [drt,dc,drd,da]
# 
# r=solve_ivp(dy,[t_min,100*t_max],y0,args=[const],t_eval=t,dense_output=True)
# 
# y=r.sol(t)
# 
# plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']
# 
# c=st.columns(2)
# for i in range(4):
#     with c[int(i/2)]:
#         source=pd.DataFrame({'t':t,plot_name[i]:y[i]})
#         chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
#         #st.altair_chart(chart)
# 
# t1=np.arange(t_min,1,0.01)
# r1=solve_ivp(dy,[t_min,1],[y[0][-1],y[1][-1],y[2][-1],y[3][-1]],t_eval=t1,args=[const],dense_output=True)
# y1=r1.sol(t1)
# 
# t2=np.arange(1,t_max,0.01)
# r2=solve_ivp(dy,[1,t_max],[y1[0][-1],y1[1][-1],y1[2][-1],fact*y1[3][-1]],args=[const],t_eval=t2,dense_output=True)
# y2=r2.sol(t2)
# 
# plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']
# 
# for i in range(4):
#     with c[int(i/2)]:
#         source=pd.DataFrame({'t':list(t1)+list(t2),plot_name[i]:list(y1[i])+list(y2[i])})
#         chart=alt.Chart(source).mark_line().encode(x='t',y=plot_name[i]).interactive()
#         st.altair_chart(chart)
# =============================================================================
# =============================================================================
#         fig=plt.figure()
#         plt.subplot(2,2,i+1)
#         plt.plot(t,y[i])
#         st.pyplot(fig)
# =============================================================================
#st.pyplot(fig)


def f(const):
    k1,k2,k3,k4,k5,k6,kon,koff,a0=const
    rtot=1
    
    a=(-k6*k5*kon*a0/(k4*koff))*(1+k5*kon*a0/(koff*k4))
    b=-(k1+k2)-k5*k1*kon*a0/(k4*koff)-k3*a0*kon/koff-k5*k3*kon**2*a0**2/(k4*koff**2)+k6*k5*kon*a0*rtot/(k4*koff)
    c=(k1+k3*kon*a0/koff)*rtot
    
    return  [(-b-(b**2-4*a*c)**0.5)/(2*a),(-b+(b**2-4*a*c)**0.5)/(2*a)]

x=np.logspace(-2,2,100)
y=[f(const[0:-1]+[i])[0] for i in x]
z=[f(const[0:-1]+[i])[1] for i in x]


source=pd.DataFrame({'ARHGEF11 concentration':x,'x1 Rho-GTP':y,'x2':z})
chart=alt.Chart(source).mark_line().encode(
    alt.X('ARHGEF11 concentration',scale=alt.Scale(type='log')),
    alt.Y('x1 Rho-GTP')).interactive()
chart2=alt.Chart(source).mark_line().encode(
    alt.X('ARHGEF11 concentration',scale=alt.Scale(type='log')),
    alt.Y('x2')).interactive()

c[0].altair_chart(chart)
c[1].altair_chart(chart2)


def dy(t,y,const):
    k1,k2,k3,k4,k5,k6,kon,koff,a0=const
    rt,c,rd,a=y
    
    drt= k1*rd-k2*rt+k3*(rd*a)-k5*a*rt+k4*c+k6*c*rd 
    dc= k5*a*rt-k4*c
    drd= -k1*rd+k2*rt-k3*a*rd-k6*c*rd
    #da= -k5*a*rt+k4*c
    k=kon
    t_int=60
    t_ilid=20
    
    for i in range(5):
        k+=((10+i*t_int)<t<=(15+i*t_int))*(t-(10+i*t_int))*fact/5+(t>(15+i*t_int))*((fact))*np.exp((15+i*t_int-t)/t_ilid)
    da=k*a0-koff*a-k5*a*rt+k4*c
    
    
    return [drt,dc,drd,da]

t_min=0
t_max=100
    
const=[k1,k2,k3,k4,k5,k6,kon,koff,a0]
#y0=[1,0,0,10**a0,0]


plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH']

t2=np.arange(t_min,10*t_max,0.1)

rteq=f(const)[0]
aeq=kon*a0/koff
ceq=k5*aeq*rteq/k4 #1-rteq-rdeq
rdeq=1-rteq-ceq

r2=solve_ivp(dy,[t_min,10*t_max],[rteq,ceq,rdeq,aeq],args=[const],t_eval=t2,dense_output=True)
y2=r2.sol(t2)

plot_name=['RhoGTP','DH-PH-RhoGTP','RhoGDP','DH-PH','k']

c=st.columns(2)
for i in range(4):
    with c[int(i/2)]:
        source=pd.DataFrame({'t':list(t2),plot_name[i]:list(y2[i])})
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