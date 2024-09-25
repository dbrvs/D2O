
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import streamlit as st

#with "proliferation" as b-d in dS term
def model2(x,t,aS,pS,dS,Bt,dI,pi,gam,ton,toff,flab):
             
    Su,Sl,Iu,Il,V= x
    
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

#radio buttons for which model!
option = st.radio('Select a model to investigate:', ['Basic (aS-thS*S)', 'Proliferation (aS + (pS-dS)*S)'])

if option == 'Basic (aS-thS*S)':
    st.write("Running code for Option 1")
    # Add the code you want to run for Option 1
    
    # Streamlit interface
    st.title("Model 2, with pS-dS = thS")
    
    #fixed parameters
    aS=70
    Bt=1e-4
    dI=0.8
    pi=1e3
    gam=23
    
    # User inputs
    #aS = st.slider("alpha_S -- Source Rate [cell/µL per day]", min_value=0.0, max_value=100.0, value=70.0, step=1.0)
    #Bt = st.slider("beta -- Infection Rate [cell/virus per day]", min_value=1e-6, max_value=1e-3, value=1e-4, format="%.6f")
    #dI = st.slider("delta_I -- Infected Cell Death Rate [per day]", min_value=0.0, max_value=1.0, value=0.8, step=0.01)
    #pi = st.slider("pi -- Virus Production Rate [per day]", min_value=0.0, max_value=1e4, value=1e3, step=10.0)
    #gam = st.slider("gamma -- Virus Clearance Rate [per day]", min_value=0.0, max_value=50.0, value=23.0, step=1.0)
    
    pS = st.slider("p_S -- proliferation Rate [cell/µL per day]", min_value=0.01, max_value=1.0, value=0.1, step=0.05)
    clrS = st.slider("theta_S -- Clearance Rate [per day]", min_value=0.001, max_value=1.0, value=0.02, step=0.01)
    ton = st.slider("t_on -- Label Start [day]", min_value=0, max_value=40, value=14, step=1)
    toff = st.slider("t_off -- Label Stop [day]", min_value=0, max_value=40, value=21, step=1)
    flab = st.slider("f_label -- Labeling Fraction", min_value=0.0, max_value=1.0, value=0.3, step=0.01)
    
    dS=pS-clrS #clrS = pS-dS 
    
    # Initial conditions
    Su0 = aS / clrS
    Sl0 = 0
    Iu0 = 1e-3
    Il0 = 0
    V0 = 0
    x0 = [Su0, Sl0, Iu0, Il0, V0]
    
    # Time points
    ts = np.linspace(0, 40, 1000)
    
    # Solve the ODEs
    plist = (aS, pS, dS, Bt, dI, pi, gam, ton, toff, flab)
    sol = odeint(model2, x0, ts, plist, mxstep=10000)
    
    # Extract the results
    Su = sol[:, 0]
    Sl = sol[:, 1]
    Iu = sol[:, 2]
    Il = sol[:, 3]
    V = sol[:, 4]
    
    # Total cell counts
    u_tot = Su + Iu
    l_tot = Sl + Il
    
    # Plot the results
    st.subheader("Simulation Results")
    
    # Viral load kinetics
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    
    axs[0].semilogy(ts, V * 1e3, lw=2, color='purple')
    axs[0].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[0].set_ylabel('HIV RNA per mL')
    axs[0].set_xlabel('Days post acquisition')
    axs[0].set_ylim([1e1, 1e8])
    
    # Source cells
    axs[1].plot(ts, Su + Sl, lw=2, label='S', color='green')
    axs[1].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[1].set_ylabel('Cells per µL')
    axs[1].set_xlabel('Days post acquisition')
    axs[1].set_ylim([0, 1e3])
    
    # Labeled vs Unlabeled CD4 cells
    axs[2].plot(ts, u_tot, lw=2, label='Unlabeled CD4 cells')
    axs[2].plot(ts, l_tot, lw=2, label='Labeled CD4 cells')
    axs[2].legend()
    axs[2].set_ylabel('Cells per µL')
    axs[2].set_xlabel('Days post acquisition')
    axs[2].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[2].set_ylim([0, 1e3])
    
    # Labeled fraction
    axs[3].plot(ts, l_tot / (u_tot + l_tot) * 100, color='k', lw=2)
    axs[3].set_ylabel('Labeled fraction (%)')
    axs[3].set_xlabel('Days post acquisition')
    axs[3].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    
    plt.tight_layout()
    st.pyplot(fig)

elif option == 'Proliferation (aS + (pS-dS)*S)':
    st.write("Running code for Option 2")
    # Add the code you want to run for Option 2
    
    # Streamlit interface
    st.title("Model 2, with pS-dS = thS")
    
    #fixed parameters
    aS=70
    Bt=1e-4
    dI=0.8
    pi=1e3
    gam=23
    
    # User inputs
    #aS = st.slider("alpha_S -- Source Rate [cell/µL per day]", min_value=0.0, max_value=100.0, value=70.0, step=1.0)
    #Bt = st.slider("beta -- Infection Rate [cell/virus per day]", min_value=1e-6, max_value=1e-3, value=1e-4, format="%.6f")
    #dI = st.slider("delta_I -- Infected Cell Death Rate [per day]", min_value=0.0, max_value=1.0, value=0.8, step=0.01)
    #pi = st.slider("pi -- Virus Production Rate [per day]", min_value=0.0, max_value=1e4, value=1e3, step=10.0)
    #gam = st.slider("gamma -- Virus Clearance Rate [per day]", min_value=0.0, max_value=50.0, value=23.0, step=1.0)
    
    pS = st.slider("p_S -- proliferation Rate [cell/µL per day]", min_value=0.01, max_value=1.0, value=0.1, step=0.05)
    clrS = st.slider("theta_S -- Clearance Rate [per day]", min_value=0.001, max_value=1.0, value=0.02, step=0.01)
    ton = st.slider("t_on -- Label Start [day]", min_value=0, max_value=40, value=14, step=1)
    toff = st.slider("t_off -- Label Stop [day]", min_value=0, max_value=40, value=21, step=1)
    flab = st.slider("f_label -- Labeling Fraction", min_value=0.0, max_value=1.0, value=0.3, step=0.01)
    
    dS=pS-clrS #clrS = pS-dS 
    
    # Initial conditions
    Su0 = aS / clrS
    Sl0 = 0
    Iu0 = 1e-3
    Il0 = 0
    V0 = 0
    x0 = [Su0, Sl0, Iu0, Il0, V0]
    
    # Time points
    ts = np.linspace(0, 40, 1000)
    
    # Solve the ODEs
    plist = (aS, pS, dS, Bt, dI, pi, gam, ton, toff, flab)
    sol = odeint(model2, x0, ts, plist, mxstep=10000)
    
    # Extract the results
    Su = sol[:, 0]
    Sl = sol[:, 1]
    Iu = sol[:, 2]
    Il = sol[:, 3]
    V = sol[:, 4]
    
    # Total cell counts
    u_tot = Su + Iu
    l_tot = Sl + Il
    
    # Plot the results
    st.subheader("Simulation Results")
    
    # Viral load kinetics
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    
    axs[0].semilogy(ts, V * 1e3, lw=2, color='purple')
    axs[0].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[0].set_ylabel('HIV RNA per mL')
    axs[0].set_xlabel('Days post acquisition')
    axs[0].set_ylim([1e1, 1e8])
    
    # Source cells
    axs[1].plot(ts, Su + Sl, lw=2, label='S', color='green')
    axs[1].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[1].set_ylabel('Cells per µL')
    axs[1].set_xlabel('Days post acquisition')
    axs[1].set_ylim([0, 1e3])
    
    # Labeled vs Unlabeled CD4 cells
    axs[2].plot(ts, u_tot, lw=2, label='Unlabeled CD4 cells')
    axs[2].plot(ts, l_tot, lw=2, label='Labeled CD4 cells')
    axs[2].legend()
    axs[2].set_ylabel('Cells per µL')
    axs[2].set_xlabel('Days post acquisition')
    axs[2].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    axs[2].set_ylim([0, 1e3])
    
    # Labeled fraction
    axs[3].plot(ts, l_tot / (u_tot + l_tot) * 100, color='k', lw=2)
    axs[3].set_ylabel('Labeled fraction (%)')
    axs[3].set_xlabel('Days post acquisition')
    axs[3].axvspan(ton, toff, facecolor='gray', alpha=0.5)
    
    plt.tight_layout()
    st.pyplot(fig)
