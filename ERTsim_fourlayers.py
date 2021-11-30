
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
def createBackground(efile, xextra, botdep, shdep, marker=3):
    eloc = pd.read_csv(efile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    topo = eloc.values[:,(1,3)]

    topo[:,1] = topo[:,1] - shdep
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[topo.shape[0]-1,1]]])
        
    bots = np.array([[R_ext[0,0], e_bot-botdep+shdep],[L_ext[0,0],e_bot-botdep+shdep]])
    combo = np.concatenate((topo,R_ext,bots,L_ext))

    bpoly = mt.createPolygon(combo, isClosed=True, marker=marker)

    return topo, bpoly

# ----------Create Shallow Layer ----------------------
def createLayer(bot_topo, xextra, thx, marker):

    elev_bot = np.floor(np.min(bot_topo[:,1]))
    xL = np.floor(np.min(bot_topo[:,0]))
    xR = np.ceil(np.max(bot_topo[:,0]))

    L_ext = np.array([[xL-xextra, bot_topo[0,1]]])
    R_ext = np.array([[xR+xextra, bot_topo[bot_topo.shape[0]-1,1]]])
    topo_ext = np.concatenate((L_ext, bot_topo, R_ext))
    topo = np.copy(np.flipud(topo_ext))
    topo[:,1] = topo[:,1]+thx
    combo = np.concatenate((topo_ext, topo))

    markerpos = [topo[0,0], topo[0,1]-thx/2]
    lpoly = mt.createPolygon(combo, isClosed=True, marker=marker, markerPosition=markerpos)
    
    topo2 = np.copy(bot_topo)
    topo2[:,1] = topo2[:,1]+thx

    return topo2, lpoly
    
def createTopLayer(bot_topo, srffile, xextra, thx, marker):
    eloc = pd.read_csv(srffile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    srftopo = eloc.values[:,(1,3)]

    elev_bot = np.floor(np.min(bot_topo[:,1]))
    xL = np.floor(np.min(bot_topo[:,0]))
    xR = np.ceil(np.max(bot_topo[:,0]))

    L_ext = np.array([[xL-xextra, bot_topo[0,1]]])
    R_ext = np.array([[xR+xextra, bot_topo[bot_topo.shape[0]-1,1]]])
    topo_ext = np.flipud(np.concatenate((L_ext, bot_topo, R_ext)))  
    
    L_ext = np.array([[xL-xextra, srftopo[0,1]]])
    R_ext = np.array([[xR+xextra, srftopo[srftopo.shape[0]-1,1]]])
    srf_ext = np.concatenate((L_ext, srftopo, R_ext))

    combo = np.concatenate((topo_ext, srf_ext))

    markerpos = [srftopo[0,0], srftopo[0,1]-thx/2]
    lpoly = mt.createPolygon(combo, isClosed=True, marker=marker, markerPosition=markerpos)

    return srftopo, lpoly

# ---------- Create mesh ---------------------
def createMesh(geom, topo, sfile, Q, area=None):

    scheme = pg.DataContainerERT()
    scheme.setSensorPositions(topo)

    srv = pd.read_csv(sfile,sep='\t',header=None)
    for i, elec in enumerate("abmn"):
        scheme[elec] = srv.values[:,i+1].astype(np.int)-1


    for p in scheme.sensors():
        geom.createNode(p)
        geom.createNode(p - [0, 0.1])

    if area is None:
        mesh = mt.createMesh(geom, quality=Q)
    else:
        mesh = mt.createMesh(geom, quality=Q, area=area)

    return scheme, mesh

# ----------- Adds Geostatistical heterogeneity ---------
def addHeterogeneity(mesh, rhob, rhof, dip):
    marker = np.arange(mesh.cellCount())
    for cell in mesh.cells():
        marker[cell.id()] = cell.marker()
        
    res_b = rhob + 10 * pg.utils.generateGeostatisticalModel(mesh, I=[250,5], dip=0)
    res_f = rhof + 10 * pg.utils.generateGeostatisticalModel(mesh, I=[50,2], dip=180-dip)
    res_b[marker == 1] = res_f[marker == 1]
    
    return res_b

def addConst(mesh, rhob, rhof):
    res = np.ones(mesh.cellCount())
    return res
    

# -------------- Simulate --------------------
def runSim(mesh, scheme, rhomap, outfile):

    data = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1, noiseAbs=1e-6, seed=1337, calcOnly=True)

    data.set('r', data('u') / data('i'))
    data.set('rhoa', data('k') * data('u') / data('i'))

#    pg.info('Simulated data', data)
#    pg.info('The data contains:', data.dataMap().keys())

#    pg.info('Simulated rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
#    pg.info('Selected data noise %(min/max)', min(data['err'])*100, max(data['err'])*100)

    if outfile is not None:
        data.save(outfile)

    return data

class PHert:
    def __init__(self, thx_soil=0.5, thx_subsoil=0.5, thx_weathered=2.6, rho_soil=30, rho_subsoil=50, rho_weathered=70, rho_bedrock=90, efile=None, sfile=None, dep=200, xtra=500, Q=20, outfile=None, heterogeneity=False, area=None):
        self.dep = dep
        self.layerthx = [thx_soil, thx_subsoil, thx_weathered]
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        if efile is None:
            self.efile = DEFAULTPATH+"/data/PH-2018-eloc.txt"
        else:
            self.efile = efile
	
        if sfile is None:
            self.sfile = DEFAULTPATH+"/data/PH-2018-srv.txt"
        else:
            self.sfile = sfile

        self.heterogeneity = heterogeneity
        if not heterogeneity:
            self.rhomap = [[0, rho_soil],
                          [1, rho_subsoil],
                          [2, rho_weathered],
                          [3, rho_bedrock]]    
                
        

    def create_fault(self):
        topo, _ = createTopo(self.efile,self.xtra,self.dep)
        return createFault(topo, self.xpos, self.dip, self.H, self.dep, xtra=self.xtra)

    def create_geom(self):
        topo, bedrock = createBackground(self.efile, self.xtra, self.dep, sum(self.layerthx), marker=3)

        topo, weathered = createLayer(topo, self.xtra, self.layerthx[2], 2)

        topo, subsoil = createLayer(topo, self.xtra, self.layerthx[1], 1)   

        topo, topsoil = createTopLayer(topo, self.efile, self.xtra, self.layerthx[0], 0)
        self.topo = topo


        return bedrock+weathered+subsoil+topsoil  


    def show_geom(self):
        geom = self.create_geom()
        pg.show(geom)
        return 

    def show_mesh(self, showMesh=True, caxis=None, cmap="hsv"):
        try:
            if caxis is None:
                pg.show(self.mesh, data=self.rhomap, label=pg.unit('res'), showMesh=showMesh, cMap=cmap)
            else:
                pg.show(self.mesh, data=self.rhomap, label=pg.unit('res'), showMesh=showMesh, cMin=caxis[0], cMax=caxis[1],cMap=cmap)
        except:
            geom = self.create_geom()
            self.scheme, self.mesh = createMesh(geom, self.topo, self.sfile, self.Q, area=self.area)
            
            if self.heterogeneity:
                self.rhomap = addHeterogeneity(self.mesh, self.rho_back, self.rho_fault, self.dip)
            
            if caxis is None:
                pg.show(self.mesh, data=self.rhomap, label=pg.unit('res'), showMesh=showMesh,cMap=cmap)
            else:
                pg.show(self.mesh, data=self.rhomap, label=pg.unit('res'), showMesh=showMesh, cMin=caxis[0], cMax=caxis[1],cMap=cmap)

        return

    def run_full(self, showPlot=False, verbose=False, heterogeneity=False):
        if verbose:
            print('Creating geometry...')
        self.geom = self.create_geom()

        if verbose:
            print('Meshing...')
        self.scheme, self.mesh = createMesh(self.geom, self.topo, self.sfile, self.Q, area=self.area)

        if self.heterogeneity:
            if verbose:
                print('Adding heterogeneity...')
            self.rhomap = addHeterogeneity(self.mesh, self.rho_back, self.rho_fault, self.dip)

        if verbose:
            print('Simulating ERT...')
        self.data = runSim(self.mesh, self.scheme, self.rhomap, self.outfile)

        if showPlot:
            self.show_data()
        
        return self.data
    
    def show_data(self):
        srv = pd.read_csv(self.sfile,sep='\t',header=None)
        elec = pd.read_csv(self.efile,sep='\t',header=None)
        M = (elec.values[srv.values[:,2][:].astype(int)-1,1]+elec.values[srv.values[:,4][:].astype(int)-1,1])/2
        C = (elec.values[srv.values[:,1][:].astype(int)-1,3]+elec.values[srv.values[:,3][:].astype(int)-1,3])/2 - (elec.values[srv.values[:,4][:].astype(int)-1,1]-elec.values[srv.values[:,2][:].astype(int)-1,1])/2

#        C = (srv.values[:,4]-srv.values[:,2])/2
#        M = C + srv.values[:,2]
        R = np.array([self.data('r')[i] for i in range(self.data('r').size())])
        vmin = np.min((np.min(R),np.min(srv.values[:,5])))
        vmax = np.max((np.max(R),np.max(srv.values[:,5])))

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

        pg.show(self.geom, ax=ax1)
        ax1.set_title('Geometry')

        syn = ax2.scatter(M,C,10,c=np.log(R))
        ax2.set_title('Synthetic')
        fig.colorbar(syn, ax=ax2, label='Resistance (Ohm)')
        syn.set_clim(np.log(vmin),np.log(vmax))

        pg.show(self.mesh, ax=ax3, data=self.rhomap, label=pg.unit('res'), showMesh=True)
        ax3.set_title('Mesh')

        tru = ax4.scatter(M,C,10,c=np.log(srv.values[:,5]))
        ax4.set_title('True')
        fig.colorbar(tru, ax=ax4, label='Resistance (Ohm)')
        tru.set_clim(np.log(vmin),np.log(vmax))

        fig.tight_layout()
        plt.show()
        return

# ------------------------------------------
if __name__ == "__main__":

    PH = PHert(xpos=300,dip=15,Q=14, xtra=1000)
    #PH.show_geom()
    PH.run_full(showPlot=True)






