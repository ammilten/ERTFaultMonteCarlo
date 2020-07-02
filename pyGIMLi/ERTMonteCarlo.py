import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from ERTsim import PHert
import time

dep = 200
xtra = 1000
Q = 30
outfile = None

dip = [25, 165] #uniform
H = [120, 40] #normal
xpos = [250,40] #normal
rho_fault = [20, 10, .45] #lognormal
rho_back = [40, 8, .9] #lognormal

N = 3
MCfolder = '/home/ammilten/Programs/ERTFaultMonteCarlo/data/MC/'
showDists = False

np.random.seed()

# ------------------ Setting Distributions ---------------
dipd = st.uniform(loc=dip[0], scale=dip[1]-dip[0])
Hd = st.norm(loc=H[0], scale=H[1])
xposd = st.norm(loc=xpos[0], scale=xpos[1])
rho_faultd = st.lognorm(rho_fault[2],loc=rho_fault[0], scale=rho_fault[1])
rho_backd = st.lognorm(rho_back[2],loc=rho_back[0], scale=rho_back[1])

# ------------------- Plotting Distributions -------------
if showDists:
    fig, axs = plt.subplots(5,1)

    x = np.linspace(dipd.ppf(0.001)-10, dipd.ppf(0.999)+10, 1000)
    axs[0].plot(x,dipd.pdf(x))
    axs[0].set_xlabel('Dip')

    x = np.linspace(Hd.ppf(0.001), Hd.ppf(0.999), 1000)
    axs[1].plot(x,Hd.pdf(x))
    axs[1].set_xlabel('H')
    
    x = np.linspace(xposd.ppf(0.001), xposd.ppf(0.999), 1000)
    axs[2].plot(x,xposd.pdf(x))
    axs[2].set_xlabel('xpos')
    
    x = np.linspace(rho_faultd.ppf(0.001), rho_faultd.ppf(0.999), 1000)
    axs[3].plot(x,rho_faultd.pdf(x))
    axs[3].set_xlabel('rho_fault')
    
    x = np.linspace(rho_backd.ppf(0.001), rho_backd.ppf(0.999), 1000)
    axs[4].plot(x,rho_backd.pdf(x))
    axs[4].set_xlabel('rho_back')
    
    fig.tight_layout()
    plt.show()

# ------------------ Run Monte Carlo ---------------
dips = dipd.rvs(N)
Hs = Hd.rvs(N)
xposs = xposd.rvs(N)
rho_faults = rho_faultd.rvs(N)
rho_backs = rho_backd.rvs(N)

with open(MCfolder+'params.txt',"w") as pf:
    pf.write(str(dep)+"\t#dep\n")
    pf.write(str(xtra)+"\t#xtra\n")
    pf.write(str(Q)+"\t#Q\n")
    pf.write(str(dip[0])+", "+str(dip[1])+"\t#dip (uniform): lower, upper\n")
    pf.write(str(H[0])+", "+str(H[1])+"\t#H (normal): mean, std\n")
    pf.write(str(xpos[0])+", "+str(xpos[1])+"\t#xpos (normal): mean, std\n")
    pf.write(str(rho_fault[0])+", "+str(rho_fault[1])+", "+str(rho_fault[2])+"\t#rho_fault (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
    pf.write(str(rho_back[0])+", "+str(rho_back[1])+", "+str(rho_back[2])+"\t#rho_back (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n\n")
    pf.write(str(N)+"\t#N\n")
    pf.write("dip, H, xpos, rho_fault, rho_back\n")
    for i in range(N):
        pf.write(str(dips[i])+", "+str(Hs[i])+", "+str(xposs[i])+", "+str(rho_faults[i])+", "+str(rho_backs[i])+"\n")


tstart = time.time()
for i in range(N): 
    print("--------- Working on " + str(i+1) + " of " + str(N)+" ---------")  
    real = MCfolder + "data_" + str(i) + ".dat"
    PH = PHert(dip=dips[i], H=Hs[i], xpos=xposs[i], rho_fault=rho_faults[i], rho_back=rho_backs[i], dep=dep, xtra=xtra, Q=Q, outfile=real)
    PH.run_full()
    #PH.show_geom()

tend = time.time()
print("\nTime to simulate " + str(N) + " realizations: " + str(np.round((tend-tstart)/60,2)) + " minutes")



