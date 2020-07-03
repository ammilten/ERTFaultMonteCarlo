import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from ERTsim import PHert
from multiprocessing import Pool
import time
import sys, getopt, os

# ---------------------- Construct Distributions ---------------------
def construct_dists(dip, H, xpos, rho_fault, rho_back):
    dipd = st.uniform(loc=dip[0], scale=dip[1]-dip[0])
    Hd = st.norm(loc=H[0], scale=H[1])
    xposd = st.norm(loc=xpos[0], scale=xpos[1])
    rho_faultd = st.lognorm(rho_fault[2],loc=rho_fault[0], scale=rho_fault[1])
    rho_backd = st.lognorm(rho_back[2],loc=rho_back[0], scale=rho_back[1])

    return dipd, Hd, xposd, rho_faultd, rho_backd

# ------------- Sample & Save Parameters ---------------------
def save_params(dep, xtra, Q, dip, H, xpos, rho_fault, rho_back, N, MCfolder, override):
    dipd, Hd, xposd, rho_faultd, rho_backd = construct_dists(dip, H, xpos, rho_fault, rho_back)
    
    dips = dipd.rvs(N)
    Hs = Hd.rvs(N)
    xposs = xposd.rvs(N)
    rho_faults = rho_faultd.rvs(N)
    rho_backs = rho_backd.rvs(N)

    try:
        os.mkdir(MCfolder)
    except FileExistsError:
        if override:
            print("MC folder already exists. Overriding.")
        else:
            print("MC folder already exists. Aborting. Pass the argument -o to override")
            sys.exit(2)

    with open(MCfolder+'params.dat',"w") as pf:
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

    PARAMS = list(zip(list(range(N)), dips, Hs, xposs, rho_faults, rho_backs, [dep]*N, [xtra]*N, [Q]*N))
    return PARAMS

# ------------- Wrapper to run the simulation with one argument ---------------
def run(params):
    print("Simulating realization " + str(params[0]))
    real = real = MCfolder + "data_" + str(params[0]) + ".dat"
    PH = PHert(dip=params[1], H=params[2], xpos=params[3], rho_fault=params[4], rho_back=params[5], dep=params[6], xtra=params[7], Q=params[8], outfile=real)
    PH.run_full()
    return

# ------------- Plots the 5 distributions --------------------
def plot_dists(dip, H, xpos, rho_fault, rho_back):
    dipd, Hd, xposd, rho_faultd, rho_backd = construct_dists(dip, H, xpos, rho_fault, rho_back)

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

# ---------------- Parse arguments ---------------------
def get_args(argv):
    Q = 14
    N = 1
    showDists = False
    parallel = False
    nproc = None
    MCfolder = ''
    foundMC = False
    override = False

    try:
        opts, args = getopt.getopt(argv,"hf:q:n:p:do",["mcfolder=","quality=","nreals=","showdists","parallel","override"])
    except getopt.GetoptError:
        print("ERTMonteCarlo.py -f <path_to_mcfolder> -q <min_angle> -n <num_reals> -p <num_processors> -d -o")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("\nUSAGE:\npython ERTFaultMonteCarlo.py -f <path_to_mcfolder> -n <num_reals> -q <quality>\n\nSWITCHES:\n  -f, --mcfolder <folder>: Path to Monte Carlo output folder (required)\n  -q, --quality <min_angle>: Minimum mesh angle. 14=low quality, 34=high quality (default=14)\n  -n, --nreals <num_reals>: Number of realizations (default=1)\n  -p, --parallel <num_processors>: If specified, run in parallel. Option to specify number of of processors. If number of processors is not specified, it is auto-selected.\n  -d, --showdists: If specified, plot the parameter distributions prior to simulations.\n  -o, --override: If specified, the Monte Carlo will overwrite the data in MCfolder if MCfolder already exists.")
            sys.exit(2)
        elif opt in ("-q", "--quality"):
            Q = float(arg)
        elif opt in ("-n", "--nreals"): 
            N = int(arg)
        elif opt in ("-f", "--mcfolder"):
            MCfolder = arg
            foundMC = True
        elif opt in ("-d", "--showdists"):
            showDists = True
        elif opt in ("-o", "--override"):
            override = True
        elif opt in ("-p", "--parallel"):
            parallel = True
            try:
                nproc = int(arg)
            except:
                nproc = None
             
    if not foundMC:
        print("MCfolder argument required: use '-f <path_to_mcfolder>' or '--mcfolder <path_to_mcfolder>'")
        sys.exit(2)

    print("MCfolder = "+MCfolder)
    print("Minimum angle = " + str(Q) + "    (14 = low quality mesh, 34 = high quality mesh)")
    print("Number of realizations = " + str(N))
    if parallel:
        if nproc is None:
            print("Parallelized = true (auto-select number of processors)")
        else:
            print("Parallelized = true (" + str(nproc) + " processors)")
    else:
        print("Parallelized = false")
    if showDists:
        print("Plotting distributions = true")
    else:
        print("Plotting distributions = false")

    print("")
    return MCfolder, Q, N, showDists, parallel, nproc, override



#------------------ MAIN SCRIPT -----------------------
if __name__ == "__main__":

    print("\n----------- Preparing Monte Carlo -----------")
    MCfolder, Q, N, showDists, parallel, nproc, override = get_args(sys.argv[1:])

#    MCfolder = '/home/ammilten/Programs/ERTFaultMonteCarlo/data/MC/'

    dip = [25, 165] #uniform
    H = [120, 40] #normal
    xpos = [250,40] #normal
    rho_fault = [20, 10, .45] #lognormal
    rho_back = [40, 8, .9] #lognormal

    dep = 200
    xtra = 1000

    np.random.seed()


    # Sample & save
    PARAMS = save_params(dep, xtra, Q, dip, H, xpos, rho_fault, rho_back, N, MCfolder, override)

    # Plot distributions (optional)
    if showDists:
        plot_dists(dip, H, xpos, rho_fault, rho_back)
    
    # Simulate
    print("-------------- Simulating --------------")
    tstart = time.time()
    if parallel:
        if nproc is None:
            pool = Pool()
            pool.map(run, PARAMS)
        else:
            pool = Pool(processes=nproc)
            pool.map(run, PARAMS)
    else:
        for i in range(N): 
            #print("--------- Working on " + str(i+1) + " of " + str(N)+" ---------") 
            run(PARAMS[i])
    tend = time.time()
    print("\nTime to simulate " + str(N) + " realizations: " + str(np.round((tend-tstart)/60,2)) + " minutes")


