
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import streamlit as st

#with "proliferation" as b-d in dS term
def model2(x,t,tinf,aS,pS,dS,Bt,dI,pi,gam,ton,toff,flab):
             
    Su,Sl,Iu,Il,V= x
    
    #infection time
    if t > tinf:
        Iu+=0.001

    if t > ton and t < toff:
    
        f=flab
        #f=flab*1/(1-np.exp(rf*(t-ton))) #time vary??
        
        #proliferation moves Su to Sl (with fraction f), both types can die, both types can get infected
        dSudt = aS + pS*(1-f)*Su - dS*Su - Bt*Su*V
        dSldt = pS*f*Su + pS*f*Sl - dS*Sl - Bt*Sl*V
        
    else:

        #proliferation moves Sl to Su, flab=0 here (or could have it time vary??)
        f=0
        dSudt = aS + pS*(1-f)*Su - dS*Su - Bt*Su*V
        dSldt = -pS*Sl - dS*Sl - Bt*Sl*V
            
    dIudt = Bt*Su*V - dI*Iu
    dIldt = Bt*Sl*V - dI*Il
    
    dVdt = pi*(Iu+Il) - gam*V #assume label doesn't influence
    
    dydt = dSudt,dSldt,dIudt,dIldt,dVdt
    
    return dydt


# Streamlit interface
st.title('SIV infection with D2O labeling')

# Add explanation text below the title
st.write('This solves a mechanistic model of HIV dynamics as well as deuterium labeling')
st.write('Model equations (make this tex):')
st.write('While labeling, f=flab:')
st.latex(r'''\dot{S_u}=\alpha_S + \rho_S (1-f) S_u - \delta_S S_u - \beta S_u V''')
st.latex(r'''\dot{S_l}= \rho_S f S_u + \rho_S f S_l - \delta_S S_l - \beta S_l V''')
st.write('While NOT labeling, f=0:')
st.latex(r'''\dot{S_u}=\alpha_S + \rho_S S_u - \delta_S S_u - \beta S_u V''')
st.latex(r'''\dot{S_u}= - \rho_S S_l - \delta_S S_l- \beta S_u V''')
st.write('Rest of dynamics:')
st.latex(r'''\dot{I_u}= \beta S_u V - \delta_I I_u''')
st.latex(r'''\dot{I_l}= \beta S_l V - \delta_I I_l''')
st.latex(r'''\dot{V}= \pi (I_u+I_l) - \gamma V''')


# Add LaTeX expression

#fixed parameters
aS=40
Bt=1e-4
dI=0.8
pi=1e3
gam=23

# User inputs
#aS = st.slider("alpha_S -- Source Rate [cell/µL per day]", min_value=0.0, max_value=100.0, value=70.0, step=1.0)

pS = st.slider("p_S -- proliferation Rate [per day]", min_value=0.01, max_value=1.0, value=0.1, step=0.05)
clrS = st.slider("theta_S -- Clearance Rate [per day]", min_value=0.0, max_value=1.0, value=0.02, step=0.01)
#Bt = st.slider("beta -- Infection Rate [cell/virus per day]", min_value=1e-6, max_value=1e-3, value=1e-4, format="%.6f")
#dI = st.slider("delta_I -- Infected Cell Death Rate [per day]", min_value=0.0, max_value=1.0, value=0.8, step=0.01)
#pi = st.slider("pi -- Virus Production Rate [per day]", min_value=0.0, max_value=1e4, value=1e3, step=10.0)
#gam = st.slider("gamma -- Virus Clearance Rate [per day]", min_value=0.0, max_value=50.0, value=23.0, step=1.0)
tinf = st.slider("t_inf -- HIV inoculation [day]", min_value=0, max_value=40, value=14, step=1)
ton = st.slider("t_on -- Label Start [day]", min_value=0, max_value=40, value=14, step=1)
toff = st.slider("t_off -- Label Stop [day]", min_value=0, max_value=40, value=21, step=1)
flab = st.slider("f_label -- Labeling Fraction", min_value=0.0, max_value=1.0, value=0.3, step=0.01)

epsilonAP = st.slider("eps_AP -- APT efficacy (fold reduction of p_S)", min_value=1., max_value=10., value=2., step=0.1)

dS=pS-clrS #clrS = pS-dS 

# Initial conditions
Su0 = aS / clrS
Sl0 = 0
Iu0 = 0 #now this gets +1e-3 in the model on tinf day
Il0 = 0
V0 = 0
x0 = [Su0, Sl0, Iu0, Il0, V0]

#simulate first natural

ts = np.linspace(0,40,1000)

tinf=10
ton=0
toff=7
flab=0.3

aS=20
clrS=0.02
pS=0.1
dS=pS-clrS #clrS = pS-dS 
Bt=1e-4
dI=0.8
pi=1e3
gam=23

plist = tinf,aS,pS,dS,Bt,dI,pi,gam,ton,toff,flab

Su0=aS/clrS; Sl0=0; Iu0=0; Il0=0; V0=0;

x0 = [Su0,Sl0,Iu0,Il0,V0]

sol=odeint(model2, x0, ts, plist, mxstep=10000) 

Su = sol[:,0]; Sl = sol[:,1]; Iu = sol[:,2]; Il = sol[:,3]; V = sol[:,4]

u_tot = Su+Iu; l_tot = Sl+Il

#VL kinetics
fig = plt.figure(figsize=(12,7))

plt.subplot(221)
plt.semilogy(ts,V*1e3,lw=2,color='purple',label='Viral load')
plt.axvspan(ton, toff, facecolor='gray', alpha=0.5,label='Labeling period')
plt.ylabel('HIV RNA per mL')
plt.xlabel('Study Days')
plt.axvline(tinf,color='red',alpha=0.5,label='Inoculation')
plt.ylim([1e1,1e8])
plt.legend()

plt.subplot(222)
plt.plot(ts,Su+Sl,lw=2,label='S',color='green')
plt.axvspan(ton, toff, facecolor='gray', alpha=0.5)
#plt.legend()
plt.ylabel('Cells per µL')
plt.xlabel('Study Days')
plt.tight_layout()

#label kinetics
plt.subplot(223)
plt.plot(ts,u_tot,lw=2,label='Unlabeled CD4 cells')
plt.plot(ts,l_tot,lw=2,label='Labeled CD4 cells')
plt.legend()
plt.ylabel('Cells per µL')
plt.xlabel('Study Days')
plt.axvspan(ton, toff, facecolor='gray', alpha=0.5)
#plt.ylim([1,100])

plt.subplot(224)
plt.plot(ts,l_tot/(u_tot+l_tot)*100,color='k',lw=2,label='Normal')
plt.ylabel('Labeled fraction (%)')
plt.xlabel('Study Days')
plt.axvspan(ton, toff, facecolor='gray', alpha=0.5)
#plt.ylim([1,100])
#plt.savefig('plt.pdf',dpi=600)

#now simulate with APT

plist = tinf,aS,pS/epsilonAP,dS,Bt,dI,pi,gam,ton,toff,flab

Su0=aS/clrS; Sl0=0; Iu0=0; Il0=0; V0=0;

x0 = [Su0,Sl0,Iu0,Il0,V0]

sol=odeint(model2, x0, ts, plist, mxstep=10000) 

Su = sol[:,0]; Sl = sol[:,1]; Iu = sol[:,2]; Il = sol[:,3]; V = sol[:,4]

u_tot = Su+Iu; l_tot = Sl+Il

plt.subplot(221)
plt.semilogy(ts,V*1e3,lw=2,ls='--',color='purple')

plt.subplot(222)
plt.plot(ts,Su+Sl,lw=2,label='S',ls='--',color='green')

plt.subplot(223)
plt.plot(ts,u_tot,lw=2,ls='--',label='Unlabeled CD4 cells')
plt.plot(ts,l_tot,lw=2,ls='--',label='Labeled CD4 cells')

plt.subplot(224)
plt.plot(ts,l_tot/(u_tot+l_tot)*100,ls='--',color='k',lw=2,label='with APT')
plt.legend()

st.pyplot(fig)
