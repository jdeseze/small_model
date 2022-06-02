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

st.set_page_config(page_title="Model", page_icon=":model:", layout="wide")

t_min = 0
t_max = 10
t = np.arange(t_min, 100*t_max, 0.01)

with st.sidebar:
    k1 = st.slider('k1 (rd->rt)', 0.0, 2.0, 0.38)
    k2 = 50*k1  # st.slider('k2 (rt->rd)',0.0,2.0,0.5)
    k3 = 10*k1  # st.slider('k3 (a+rd->a+rt)',0.0,2.0,0.5)
    k4 = st.slider('k4 (c->m+a)', 0.0, 2.0, 0.5)
    k5 = st.slider('k5 (rt+m->m+rt)',0.0,2.0,0.21)
    k6 = 1*st.slider('k6 (mt->m)', 0.0, 2.0, 0.18)
    k7 = 0.01*st.slider('k7 (mt+a->c)', 0.0, 2.0, 0.24)
    kon = 0.1  # st.slider('kon ',0.0,2.0,0.5)
    koff = 1/1  # st.slider('koff',0.0,2.0,0.5)
    a0exp = st.slider('[A]0^', -3.0, 3.0, 1.14)
    fact = st.slider('factor of increase', 0.0, 10.0, 0.52)
a0 = 10**a0exp
mtot=10
rtot=1
const = [k1, k2, k3, k4, k5, k6,k7, kon, koff, 10**a0,mtot]
y0 = [1, 0, 0, 10**a0]

c = st.columns(2)

x = np.logspace(-2, 2, 100)



def dy(t, y, const):
    k1, k2, k3, k4, k5, k6,k7, kon, koff, a0, mtot,rtot = const
    rt, c, a, mt = y

    drt = k1*(rtot-rt)-k2*rt+k3*((rtot-rt)*a)
    dc = k7*a*mt-k4*c
    dmt = k5*rt*(mtot-mt)-k6*mt-k7*a*mt
    #da= -k5*a*rt+k4*c
    k = kon
    t_int = 10
    t_ilid = 20

    for i in range(8):
        k += ((10+i*t_int) < t <= (15+i*t_int))*(t-(10+i*t_int))*fact / \
            5+(t > (15+i*t_int))*((fact))*np.exp((15+i*t_int-t)/t_ilid)
    da = k*a0-koff*a

    return [drt, dc, da, dmt]


t_min = 0
t_max = 200

const = [k1, k2, k3, k4, k5, k6,k7, kon, koff, a0, mtot,rtot]
# y0=[1,0,0,10**a0,0]


t2 = np.arange(t_min, t_max, 0.1)


aeq = kon*a0/koff
rteq = (k1*rtot+k3*aeq*rtot)/(k1+k2+k3*aeq)
rdeq = rtot-rteq
mteq = k5*rteq*mtot/(k5*rteq+k6+k7*aeq)
ceq = k7*aeq*mteq/k4

r2 = solve_ivp(dy, [t_min, 10*t_max], [rteq,ceq,aeq,mteq], args=[const], t_eval=t2, dense_output=True,method='BDF')
y2 = r2.sol(t2)

plot_name = ['RhoGTP', 'DH-PH-Myosin', 'DH-PH', 'Myosin']

c = st.columns(2)
for i in range(4):
    with c[i%2]:
        source = pd.DataFrame({'t': list(t2), plot_name[i]: list(y2[i]/y2[i][0])})
        chart = alt.Chart(source).mark_line().encode(
            x='t', y=plot_name[i]).interactive()
        st.altair_chart(chart)
