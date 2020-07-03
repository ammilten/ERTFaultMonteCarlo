
# -----------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate 

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert

import pathlib
from os import fspath
DEFAULTPATH = fspath(pathlib.Path(__file__).parent.absolute())

# ----------Create Topo Polygon ----------------------
def createTopo(efile, xextra, dep):
    eloc = pd.read_csv(efile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    
    topo = eloc.values[:,(1,3)]
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[nelec-1,1]]])
    bots = np.array([[R_ext[0,0], e_bot-dep],[L_ext[0,0],e_bot-dep]])
    combo = np.concatenate((topo,R_ext,bots,L_ext))

    tpoly = mt.createPolygon(combo, isClosed=True, marker=2)

    return topo, tpoly

# --------- Create Fault Polygon ------------
def createFault(topo, xpos, dip, H, dep, xtra=500):
    topo2 = np.concatenate((np.array([[topo[0,0]-xtra, topo[0,1]]]), topo, np.array([[topo[topo.shape[0]-1,0]+xtra,topo[topo.shape[0]-1,1]]])))

    e_bot = np.floor(np.min(topo2[:,1]))
    zbot = np.floor(np.min(topo2[:,1])) - dep

    Z = interpolate.interp1d(topo2[:,0],topo2[:,1])
    zx = Z(xpos)
    Dz = zx-zbot

    Dx = Dz/np.tan(np.deg2rad(dip))
    xbot = xpos + Dx

    Dxbot = H/2/np.sin(np.deg2rad(dip))

    LR = (xbot+Dxbot,zbot)
    LL = (xbot-Dxbot,zbot)

    zfR = zbot + (LR[0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffR = interpolate.interp1d(zfR-topo2[:,1],topo2[:,0])
    UR = (float(diffR(0)), float(Z(diffR(0))) )

    zfL = zbot + (LL[0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffL = interpolate.interp1d(zfL-topo2[:,1],topo2[:,0])
    UL = (float(diffL(0)), float(Z(diffL(0))) )
    
    idxabove = (topo2[topo2[:,0] > UL[0],0]).argmin() + (topo2[topo2[:,0] < UL[0],0]).shape[0]
    idxbelow = (topo2[topo2[:,0] < UR[0],0]).argmax()

    middles = [(topo[j,0],topo[j,1]) for j in range(idxabove,idxbelow)]
    verts = [LL, UL] + middles + [UR,LR]
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=3, marker=1)

    return fpoly

def createFault_old(topo, xpos, dip, H, dep):
    e_bot = np.floor(np.min(topo[:,1]))
    zbot = np.floor(np.min(topo[:,1])) - dep

    Z = interpolate.interp1d(topo[:,0],topo[:,1])
    zx = Z(xpos)
    Dz = zx-zbot

    Dx = Dz/np.tan(np.deg2rad(dip))
    xbot = xpos + Dx

    Dxbot = H/2/np.sin(np.deg2rad(dip))

    LR = (xbot+Dxbot,zbot)
    LL = (xbot-Dxbot,zbot)

    zfR = zbot + (LR[0]-topo[:,0])*np.tan(np.deg2rad(dip))
    diffR = interpolate.interp1d(zfR-topo[:,1],topo[:,0])
    UR = (float(diffR(0)), float(Z(diffR(0))) )

    zfL = zbot + (LL[0]-topo[:,0])*np.tan(np.deg2rad(dip))
    diffL = interpolate.interp1d(zfL-topo[:,1],topo[:,0])
    UL = (float(diffL(0)), float(Z(diffL(0))) )
    
    idxabove = (topo[topo[:,0] > UL[0],0]).argmin() + (topo[topo[:,0] < UL[0],0]).shape[0]
    idxbelow = (topo[topo[:,0] < UR[0],0]).argmax()

    middles = [(topo[j,0],topo[j,1]) for j in range(idxabove,idxbelow)]
    verts = [LL, UL] + middles + [UR,LR]
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=3, marker=1)

    return fpoly
# ---------- Create mesh ---------------------
def createMesh(geom, topo, sfile, Q):

    scheme = pg.DataContainerERT()
    scheme.setSensorPositions(topo)

    srv = pd.read_csv(sfile,sep='\t',header=None)
    for i, elec in enumerate("abmn"):
        scheme[elec] = srv.values[:,i+1].astype(np.int)-1


    for p in scheme.sensors():
        geom.createNode(p)
        geom.createNode(p - [0, 0.1])

    mesh = mt.createMesh(geom, quality=Q)

    return scheme, mesh

# -------------- Simulate --------------------
def runSim(mesh, scheme, rho_fault, rho_back, outfile):
    rhomap = [[0, rho_back],
              [1, rho_fault],
              [2, rho_back]]

    data = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1, noiseAbs=1e-6, seed=1337, calcOnly=True)

    data.set('r', data('u') / data('i'))
    data.set('rhoa', data('k') * data('u') / data('i'))

#    pg.info('Simulated data', data)
#    pg.info('The data contains:', data.dataMap().keys())

#    pg.info('Simulated rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
#    pg.info('Selected data noise %(min/max)', min(data['err'])*100, max(data['err'])*100)

    if outfile is not None:
        data.save(outfile)

    return data, rhomap

class PHert:
    def __init__(self, dip=60, H=100, xpos=250, rho_fault=30, rho_back=50, efile=None, sfile=None, dep=200, xtra=500, Q=20, outfile=None):
        self.dep = dep
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.rho_fault = rho_fault
        self.rho_back = rho_back
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        if efile is None:
            self.efile = DEFAULTPATH+"/data/PH-2018-eloc.txt"
        else:
            self.efile = efile

        if sfile is None:
            self.sfile = DEFAULTPATH+"/data/PH-2018-srv.txt"
        else:
            self.sfile = sfile

    def create_geom(self):
        topo, tpoly = createTopo(self.efile,self.xtra,self.dep)
        fpoly = createFault(topo, self.xpos, self.dip, self.H, self.dep, xtra=self.xtra)
        return tpoly+fpoly 

    def show_geom(self):
        geom = self.create_geom()
        pg.show(geom)
        return 

    def run_full(self, showPlot=False):
        print('Creating topo...')
        topo, tpoly = createTopo(self.efile,self.xtra,self.dep)

        print('Creating fault...')
        fpoly = createFault(topo, self.xpos, self.dip, self.H, self.dep)

        print('Meshing...')
        scheme, mesh = createMesh(tpoly+fpoly, topo, self.sfile, self.Q)

        print('Simulating ERT...')
        data, rhomap = runSim(mesh, scheme, self.rho_fault, self.rho_back, self.outfile)

        if showPlot:
            srv = pd.read_csv(self.sfile,sep='\t',header=None)
            C = (srv.values[:,4]-srv.values[:,2])/2
            M = C + srv.values[:,2]
            R = np.array([data('r')[i] for i in range(data('r').size())])
            vmin = np.min((np.min(R),np.min(srv.values[:,5])))
            vmax = np.max((np.max(R),np.max(srv.values[:,5])))

            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

            pg.show(tpoly+fpoly, ax=ax1)
            ax1.set_title('Geometry')

            syn = ax2.scatter(M,-C,10,c=np.log(R))
            ax2.set_title('Synthetic')
            fig.colorbar(syn, ax=ax2, label='Resistance (Ohm)')
            syn.set_clim(np.log(vmin),np.log(vmax))

            pg.show(mesh, ax=ax3, data=rhomap, label=pg.unit('res'), showMesh=True)
            ax3.set_title('Mesh')

            tru = ax4.scatter(M,-C,10,c=np.log(srv.values[:,5]))
            ax4.set_title('True')
            fig.colorbar(tru, ax=ax4, label='Resistance (Ohm)')
            tru.set_clim(np.log(vmin),np.log(vmax))

            fig.tight_layout()
            plt.show()
        
        return data


# ------------------------------------------
if __name__ == "__main__":

    PH = PHert(xpos=300,dip=15,Q=14, xtra=1000)
    #PH.show_geom()
    PH.run_full(showPlot=True)






