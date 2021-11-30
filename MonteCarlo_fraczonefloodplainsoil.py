import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from ERTsim_fraczonefloodplainsoil import PHert
from multiprocessing import Pool
import os, sys, time
from os import path


# ------------- Wrapper to run a realization ---------------
def realize(params):
    print("Simulating realization " + str(params[0]))
    MCfolder = params[13]
    real = MCfolder + "data_" + str(params[0]) + ".dat"

    PH = PHert(dip=params[1], H=params[2], xpos=params[3], rho_fault=params[4], rho_back=params[5], shdep=params[6], rho_sh=params[7], fpthx=params[8], rho_fp=params[9],dep=params[10], xtra=params[11], Q=params[12], outfile=real)
    PH.run_full()
    return

# ------------- Parse Parameters ---------------------
def import_simulation(MCfolder, overwrite=False, parallel=False, nproc=None, showDists=False, hasSeed=True): 
    if MCfolder.endswith('/'):
        pfile = MCfolder + 'params.dat'
    else:
        pfile = MCfolder + '/params.dat'
    if not path.exists(pfile):
        print("Error: " + pfile + " does not exist. Aborting")
        sys.exit(2)

    with open(pfile,"r") as pf: # Might need to re-do to account for whitespace
        comm = '#'
        sep = ','

        dep = float(pf.readline().split(comm,1)[0])
        xtra = float(pf.readline().split(comm,1)[0])
        Q = float(pf.readline().split(comm,1)[0])

        dipln = pf.readline().split(comm,1)[0].split(sep)
        dip = [float(dipln[0]), float(dipln[1])]

        Hln = pf.readline().split(comm,1)[0].split(sep)
        H = [float(Hln[0]), float(Hln[1])]

        xposln = pf.readline().split(comm,1)[0].split(sep)
        xpos = [float(xposln[0]), float(xposln[1])]

        rho_faultln = pf.readline().split(comm,1)[0].split(sep)
        rho_fault = [float(rho_faultln[0]), float(rho_faultln[1]), float(rho_faultln[2])]

        rho_backln = pf.readline().split(comm,1)[0].split(sep)
        rho_back = [float(rho_backln[0]), float(rho_backln[1]), float(rho_backln[2])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        shdep = [float(ln[0]), float(ln[1])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_sh = [float(ln[0]), float(ln[1]), float(ln[2])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        fpthx = [float(ln[0]), float(ln[1])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_fp = [float(ln[0]), float(ln[1]), float(ln[2])]

        pf.readline() # THIS LINE ASSUMES THERE IS A SPACE
        if hasSeed:
            seed = int(pf.readline().split(comm,1)[0])
        else:
            seed = 0
        N = int(pf.readline().split(comm,1)[0])

        pf.readline() # Skip the header
        dips = [None]*N
        Hs = [None]*N
        xposs = [None]*N
        rho_faults = [None]*N
        rho_backs = [None]*N
        shdeps = [None]*N
        rho_shs = [None]*N
        fpthxs = [None]*N
        rho_fps = [None]*N
        for i in range(N):
            params = pf.readline().split(sep)
            dips[i] = float(params[0])
            Hs[i] = float(params[1])
            xposs[i] = float(params[2])
            rho_faults[i] = float(params[3])
            rho_backs[i] = float(params[4])
            shdeps[i] = float(params[5])
            rho_shs[i] = float(params[6])
            fpthxs[i] = float(params[7])
            rho_fps[i] = float(params[8])

        PARAMS = list(zip(list(range(N)), dips, Hs, xposs, rho_faults, rho_backs, shdeps, rho_shs, fpthxs, rho_fps, [dep]*N, [xtra]*N, [Q]*N, [MCfolder]*N))

    MC = MonteCarlo(dep, xtra, Q, dip, H, xpos, rho_fault, rho_back, shdep, rho_sh, fpthx, rho_fp, N, MCfolder, overwrite=overwrite, parallel=parallel, nproc=nproc, showDists=showDists, seed=seed, saveParams=False) 
    MC.PARAMS = PARAMS

    return MC

class MonteCarlo:
    def __init__(self, dep, xtra, Q, dip, H, xpos, rho_fault, rho_back, shdep, rho_sh, fpthx, rho_fp, N, MCfolder, overwrite=False, parallel=False, nproc=None, showDists=False, seed=0, saveParams=True):
        self.dep = dep
        self.xtra = xtra
        self.Q = Q
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.rho_fault = rho_fault
        self.rho_back = rho_back
        self.shdep =  shdep
        self.rho_sh = rho_sh
        self.fpthx = fpthx
        self.rho_fp = rho_fp 
        self.N = N
        self.overwrite = overwrite
        self.parallel = parallel
        self.nproc = nproc
        self.showDists = showDists  
        if MCfolder.endswith('/'):
            self.MCfolder = MCfolder
        else:
            self.MCfolder = MCfolder + '/'
        self.parameterfile = self.MCfolder + 'params.dat'
        self.seed = seed

        np.random.seed(self.seed)
        self.construct_dists()
        if saveParams:
            self.save_params()
        return

    # ---------------------- Construct Distributions ---------------------
    def construct_dists(self):
        self.dipd = st.uniform(loc=self.dip[0], scale=self.dip[1]-self.dip[0])
        self.Hd = st.norm(loc=self.H[0], scale=self.H[1])
        self.xposd = st.norm(loc=self.xpos[0], scale=self.xpos[1])
        self.rho_faultd = st.lognorm(self.rho_fault[2],loc=self.rho_fault[0], scale=self.rho_fault[1])
        self.rho_backd = st.lognorm(self.rho_back[2],loc=self.rho_back[0], scale=self.rho_back[1])
        self.shdepd = st.uniform(loc=self.shdep[0], scale=self.shdep[1]-self.shdep[0])
        self.rho_shd = st.lognorm(self.rho_sh[2],loc=self.rho_sh[0], scale=self.rho_sh[1])
        self.fpthxd = st.uniform(loc=self.fpthx[0], scale=self.fpthx[1]-self.fpthx[0])
        self.rho_fpd = st.lognorm(self.rho_fp[2],loc=self.rho_fp[0], scale=self.rho_fp[1])
        return

    # ------------- Sample & Save Parameters ---------------------
    def save_params(self):     
    
        dips = self.dipd.rvs(self.N)
        Hs = self.Hd.rvs(self.N)
        xposs = self.xposd.rvs(self.N)
        rho_faults = self.rho_faultd.rvs(self.N)
        rho_backs = self.rho_backd.rvs(self.N)
        shdeps = self.shdepd.rvs(self.N)
        rho_shs = self.rho_shd.rvs(self.N)
        fpthxs = self.fpthxd.rvs(self.N)
        rho_fps = self.rho_fpd.rvs(self.N)
        

        try:
            os.mkdir(self.MCfolder)
        except FileExistsError:
            if self.overwrite:
                print("MC folder already exists. Overwriting.")
            else:
                print("MC folder already exists. Only simulating incomplete realizations. Pass the switch '-o' to overwrite all realizations.")
    
        with open(self.parameterfile,"w") as pf:
            pf.write(str(self.dep)+"\t#dep\n")
            pf.write(str(self.xtra)+"\t#xtra\n")
            pf.write(str(self.Q)+"\t#Q\n")
            pf.write(str(self.dip[0])+", "+str(self.dip[1])+"\t#dip (uniform): lower, upper\n")
            pf.write(str(self.H[0])+", "+str(self.H[1])+"\t#H (normal): mean, std\n")
            pf.write(str(self.xpos[0])+", "+str(self.xpos[1])+"\t#xpos (normal): mean, std\n")
            pf.write(str(self.rho_fault[0])+", "+str(self.rho_fault[1])+", "+str(self.rho_fault[2])+"\t#rho_fault (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.rho_back[0])+", "+str(self.rho_back[1])+", "+str(self.rho_back[2])+"\t#rho_back (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.shdep[0])+", "+str(self.shdep[1])+"\t#shdep (uniform): lower, upper\n")
            pf.write(str(self.rho_sh[0])+", "+str(self.rho_sh[1])+", "+str(self.rho_sh[2])+"\t#rho_sh (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.fpthx[0])+", "+str(self.fpthx[1])+"\t#fpthx (uniform): lower, upper\n")
            pf.write(str(self.rho_fp[0])+", "+str(self.rho_fp[1])+", "+str(self.rho_fp[2])+"\t#rho_sh (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n\n")
            pf.write(str(self.seed)+"\t#random number seed\n")
            pf.write(str(self.N)+"\t#N\n")
            pf.write("dip, H, xpos, rho_fault, rho_back\n")
            for i in range(self.N):
                pf.write(str(dips[i])+", "+str(Hs[i])+", "+str(xposs[i])+", "+str(rho_faults[i])+", "+str(rho_backs[i])+", "+str(shdeps[i])+", "+str(rho_shs[i])+", "+str(fpthxs[i])+", "+str(rho_fps[i])+"\n")

        self.PARAMS = list(zip(list(range(self.N)), dips, Hs, xposs, rho_faults, rho_backs, shdeps, rho_shs, fpthxs, rho_fps, [self.dep]*self.N, [self.xtra]*self.N, [self.Q]*self.N, [self.MCfolder]*self.N))
        return self.PARAMS


    # ------------- Plots the 9 distributions --------------------
    def plot_dists(self, figsize=(6,18)):
        fig, axs = plt.subplots(9,1, figsize=figsize)

        x = np.linspace(self.dipd.ppf(0.001)-10, self.dipd.ppf(0.999)+10, 1000)
        axs[0].plot(x,self.dipd.pdf(x))
        axs[0].set_xlabel('Dip')
    
        x = np.linspace(self.Hd.ppf(0.001), self.Hd.ppf(0.999), 1000)
        axs[1].plot(x,self.Hd.pdf(x))
        axs[1].set_xlabel('H')
    
        x = np.linspace(self.xposd.ppf(0.001), self.xposd.ppf(0.999), 1000)
        axs[2].plot(x,self.xposd.pdf(x))
        axs[2].set_xlabel('xpos')
    
        x = np.linspace(self.rho_faultd.ppf(0.001), self.rho_faultd.ppf(0.999), 1000)
        axs[3].plot(x,self.rho_faultd.pdf(x))
        axs[3].set_xlabel('rho_fault')
    
        x = np.linspace(self.rho_backd.ppf(0.001), self.rho_backd.ppf(0.999), 1000)
        axs[4].plot(x,self.rho_backd.pdf(x))
        axs[4].set_xlabel('rho_back')
        
        x = np.linspace(self.shdepd.ppf(0.001)-1, self.shdepd.ppf(0.999)+1, 1000)
        axs[5].plot(x,self.shdepd.pdf(x))
        axs[5].set_xlabel('shdep')
        
        x = np.linspace(self.rho_shd.ppf(0.001), self.rho_shd.ppf(0.999), 1000)
        axs[6].plot(x,self.rho_shd.pdf(x))
        axs[6].set_xlabel('rho_sh')
        
        x = np.linspace(self.fpthxd.ppf(0.001)-1, self.fpthxd.ppf(0.999)+1, 1000)
        axs[7].plot(x,self.fpthxd.pdf(x))
        axs[7].set_xlabel('fpthx')
        
        x = np.linspace(self.rho_fpd.ppf(0.001), self.rho_fpd.ppf(0.999), 1000)
        axs[8].plot(x,self.rho_fpd.pdf(x))
        axs[8].set_xlabel('rho_fp')
    
        fig.tight_layout()
        plt.show()
        return

    def prepare_realization(self,i):
        real = self.MCfolder + "data_" + str(self.PARAMS[i][0]) + ".dat"
        PH = PHert(dip=self.PARAMS[i][1], H=self.PARAMS[i][2], xpos=self.PARAMS[i][3], rho_fault=self.PARAMS[i][4], rho_back=self.PARAMS[i][5], shdep=self.PARAMS[i][6], rho_sh=self.PARAMS[i][7], fpthx=self.PARAMS[i][8], rho_fp=self.PARAMS[i][9], dep=self.PARAMS[i][10], xtra=self.PARAMS[i][11], Q=self.PARAMS[i][12], outfile=real)
        return PH

    # ------------- Wrapper to run a realization ---------------
    def realize(self, i):
        print("Simulating realization " + str(self.PARAMS[i][0]))
        PH = self.prepare_realization(i)
        PH.run_full()
        return PH

    # ------------- Find realizations that have already been run -------
    def find_incomplete_reals(self):
        # find existing data files
        incomplete = []
        for i in range(self.N):
            real = self.MCfolder + 'data_' + str(i) + '.dat'
            if not path.exists(real):
                incomplete.append(i)
        return incomplete


    # ------------------ Runs full Monte Carlo ----------------
    def run(self):
        incomplete = self.find_incomplete_reals()
        if self.overwrite:
            print("All "+str(self.N) + " Realizations will be overwritten.")
            PARAMS = self.PARAMS
        else:
            if len(incomplete) > 0:
                print(str(len(incomplete)) + " incomplete realizations will be simulated.")
                PARAMS = [self.PARAMS[i] for i in incomplete]
            else:
                print("0 incomplete realizations. Exiting.")
                sys.exit()

    
        # Plot distributions (optional)
        if self.showDists:
            self.plot_dists()
    
        # Simulate
        print("-------------- Simulating " + str(len(PARAMS)) + " Realizations --------------")
        tstart = time.time()
        if self.parallel:
            if self.nproc is None:
                pool = Pool()
                pool.map(realize, PARAMS)
            else:
                pool = Pool(processes=self.nproc)
                pool.map(realize, PARAMS)
        else:
            for i in range(len(PARAMS)): 
                #print("--------- Working on " + str(i+1) + " of " + str(N)+" ---------") 
                realize(PARAMS[i])
        tend = time.time()
        print("\nTime to simulate " + str(len(PARAMS)) + " realizations: " + str(np.round((tend-tstart)/60,2)) + " minutes")
        return

