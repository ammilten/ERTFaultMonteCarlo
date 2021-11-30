import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from ERTsim_fourlayers import PHert
from multiprocessing import Pool
import os, sys, time
from os import path


# ------------- Wrapper to run a realization ---------------
def realize(params):
    print("Simulating realization " + str(params[0]))
    print(params)
    MCfolder = params[11]
    real = MCfolder + "data_" + str(params[0]) + ".dat"

    PH = PHert(dip=params[1], H=params[2], xpos=params[3], rho_fault=params[4], rho_back=params[5], dep=params[6], xtra=params[7], Q=params[8], outfile=real)
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

        ln = pf.readline().split(comm,1)[0].split(sep)
        thx_soil = [float(thx_ln[0]), float(thx_ln[1])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        thx_subsoil = [float(thx_ln[0]), float(thx_ln[1])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        thx_weathered = [float(thx_ln[0]), float(thx_ln[1])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_soil = [float(ln[0]), float(ln[1]), float(ln[2])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_subsoil = [float(ln[0]), float(ln[1]), float(ln[2])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_weathered = [float(ln[0]), float(ln[1]), float(ln[2])]
        
        ln = pf.readline().split(comm,1)[0].split(sep)
        rho_bedrock = [float(ln[0]), float(ln[1]), float(ln[2])]

        pf.readline() # THIS LINE ASSUMES THERE IS A SPACE
        if hasSeed:
            seed = int(pf.readline().split(comm,1)[0])
        else:
            seed = 0
        N = int(pf.readline().split(comm,1)[0])

        pf.readline() # Skip the header
        thx_soil = [None]*N
        thx_subsoil = [None]*N
        thx_weathered = [None]*N
        rho_soil = [None]*N
        rho_subsoil = [None]*N
        rho_weathered = [None]*N
        rho_bedrock = [None]*N
        for i in range(N):
            params = pf.readline().split(sep)
            thx_soils[i] = float(params[0])
            thx_subsoils[i] = float(params[1])
            thx_weathereds[i] = float(params[2])
            rho_soils[i] = float(params[3])
            rho_subsoils[i] = float(params[4])
            rho_weathereds[i] = float(params[5])
            rho_bedrocks[i] = float(params[6])
            
        PARAMS = list(zip(list(range(N)), thx_soils, thx_subsoils, thx_weathereds, rho_soils, rho_subsoils, rho_weathereds, rho_bedrocks, [dep]*N, [xtra]*N, [Q]*N, [MCfolder]*N))

    MC = MonteCarlo(dep, xtra, Q, thx_soil, thx_subsoil, thx_weathered, rho_soil, rho_subsoil, rho_weathered, rho_bedrock, N, MCfolder, overwrite=overwrite, parallel=parallel, nproc=nproc, showDists=showDists, seed=seed, saveParams=False) 
    MC.PARAMS = PARAMS

    return MC

class MonteCarlo:
    def __init__(self, dep, xtra, Q, thx_soil, thx_subsoil, thx_weathered, rho_soil, rho_subsoil, rho_weathered, rho_bedrock, N, MCfolder, overwrite=False, parallel=False, nproc=None, showDists=False, seed=0, saveParams=True):
        self.thx_soil = thx_soil
        self.thx_subsoil = thx_subsoil
        self.thx_weathered = thx_weathered
        self.rho_soil = rho_soil
        self.rho_subsoil = rho_subsoil
        self.rho_weathered = rho_weathered
        self.rho_bedrock = rho_bedrock
        self.dep = dep
        self.xtra = xtra
        self.Q = Q
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
        self.thx_soild = st.uniform(loc=self.thx_soil[0], scale=self.thx_soil[1]-self.thx_soil[0])
        self.thx_subsoild = st.uniform(loc=self.thx_subsoil[0], scale=self.thx_subsoil[1]-self.thx_subsoil[0])
        self.thx_weatheredd = st.uniform(loc=self.thx_weathered[0], scale=self.thx_weathered[1]-self.thx_weathered[0])


        self.rho_soild = st.lognorm(self.rho_soil[2],loc=self.rho_soil[0], scale=self.rho_soil[1])
        self.rho_subsoild = st.lognorm(self.rho_subsoil[2],loc=self.rho_subsoil[0], scale=self.rho_subsoil[1])
        self.rho_weatheredd = st.lognorm(self.rho_weathered[2],loc=self.rho_weathered[0], scale=self.rho_weathered[1])
        self.rho_bedrockd = st.lognorm(self.rho_bedrock[2],loc=self.rho_bedrock[0], scale=self.rho_bedrock[1])

        return

    # ------------- Sample & Save Parameters ---------------------
    def save_params(self):     
    
        thx_soils = self.thx_soild.rvs(self.N)
        thx_subsoils = self.thx_subsoild.rvs(self.N)
        thx_weathereds = self.thx_weatheredd.rvs(self.N)
        
        rho_soils = self.rho_soild.rvs(self.N)
        rho_subsoils = self.rho_subsoild.rvs(self.N)
        rho_weathereds = self.rho_weatheredd.rvs(self.N)
        rho_bedrocks = self.rho_bedrockd.rvs(self.N)
 

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
            pf.write(str(self.thx_soil[0]) + "," + str(self.thx_soil[1]) + "\t#thx_soil (uniform): lower, upper\n")
            pf.write(str(self.thx_subsoil[0]) + "," + str(self.thx_subsoil[1]) + "\t#thx_subsoil (uniform): lower, upper\n")
            pf.write(str(self.thx_weathered[0]) + "," + str(self.thx_weathered[1]) + "\t#thx_weathered (uniform): lower, upper\n")
            pf.write(str(self.rho_soil[0])+", "+str(self.rho_soil[1])+", "+str(self.rho_soil[2])+"\t#rho_soil (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.rho_subsoil[0])+", "+str(self.rho_subsoil[1])+", "+str(self.rho_subsoil[2])+"\t#rho_subsoil (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.rho_weathered[0])+", "+str(self.rho_weathered[1])+", "+str(self.rho_weathered[2])+"\t#rho_weathered (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.rho_bedrock[0])+", "+str(self.rho_bedrock[1])+", "+str(self.rho_bedrock[2])+"\t#rho_bedrock (lognormal): loc, scale, s (see scipy.stats.lognorm docs)\n")
            pf.write(str(self.seed)+"\t#random number seed\n")
            pf.write(str(self.N)+"\t#N\n")
            pf.write("thx_soil, thx_subsoil, thx_weathered, rho_soil, rho_subsoil, rho_weathered, rho_bedrock\n")
            for i in range(self.N):
                pf.write(str(thx_soils[i])+", "+str(thx_subsoils[i])+", "+str(thx_weathereds[i])+", "+str(rho_soils[i])+", "+str(rho_subsoils[i])+", "+str(rho_weathereds[i])+", "+str(rho_bedrocks[i])+"\n")

        self.PARAMS = list(zip(list(range(self.N)), thx_soils, thx_subsoils, thx_weathereds, rho_soils, rho_subsoils, rho_weathereds, rho_bedrocks, [self.dep]*self.N, [self.xtra]*self.N, [self.Q]*self.N, [self.MCfolder]*self.N))
        return self.PARAMS


    # ------------- Plots the 7 distributions --------------------
    def plot_dists(self, figsize=(6, 14)):
        fig, axs = plt.subplots(7,1, figsize=figsize)

        x = np.linspace(self.thx_soild.ppf(0.001)-0.1, self.thx_soild.ppf(0.999)+0.1, 1000)
        axs[0].plot(x,self.thx_soild.pdf(x))
        axs[0].set_xlabel('thx_soil')
        
        x = np.linspace(self.thx_subsoild.ppf(0.001)-0.1, self.thx_soild.ppf(0.999)+0.1, 1000)
        axs[1].plot(x,self.thx_subsoild.pdf(x))
        axs[1].set_xlabel('thx_subsoil')
        
        x = np.linspace(self.thx_weatheredd.ppf(0.001)-0.1, self.thx_weatheredd.ppf(0.999)+0.1, 1000)
        axs[2].plot(x,self.thx_weatheredd.pdf(x))
        axs[2].set_xlabel('thx_weathered')
        
        print(self.rho_soild)
        x = np.linspace(self.rho_soild.ppf(0.001)-5, self.rho_soild.ppf(0.999)+5, 1000)
        axs[3].plot(x,self.rho_soild.pdf(x))
        axs[3].set_xlabel('rho_soil')
        
        x = np.linspace(self.rho_subsoild.ppf(0.001)-5, self.rho_subsoild.ppf(0.999)+5, 1000)
        axs[4].plot(x,self.rho_subsoild.pdf(x))
        axs[4].set_xlabel('rho_subsoil')
        
        x = np.linspace(self.rho_weatheredd.ppf(0.001)-5, self.rho_weatheredd.ppf(0.999)+5, 1000)
        axs[5].plot(x,self.rho_weatheredd.pdf(x))
        axs[5].set_xlabel('rho_weathered')
        
        x = np.linspace(self.rho_bedrockd.ppf(0.001)-5, self.rho_bedrockd.ppf(0.999)+5, 1000)
        axs[6].plot(x,self.rho_bedrockd.pdf(x))
        axs[6].set_xlabel('rho_bedrock')

    
        fig.tight_layout()
        plt.show()
        return

    def prepare_realization(self,i):
        real = self.MCfolder + "data_" + str(self.PARAMS[i][0]) + ".dat"
        PH = PHert(
          thx_soil=self.PARAMS[i][1],
          thx_subsoil=self.PARAMS[i][2],
          thx_weathered=self.PARAMS[i][3],
          rho_soil=self.PARAMS[i][4],
          rho_subsoil=self.PARAMS[i][5],
          rho_weathered=self.PARAMS[i][6],
          rho_bedrock=self.PARAMS[i][7],
          dep=self.PARAMS[i][8],
          xtra=self.PARAMS[i][9],
          Q=self.PARAMS[i][10], 
          outfile=real)
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

