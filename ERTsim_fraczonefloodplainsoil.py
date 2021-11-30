
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
def createBackground(efile, xextra, botdep, shdep=0, marker=0):
    eloc = pd.read_csv(efile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    topo = eloc.values[:,(1,3)]
    
#    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
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
def createLayer(efile, xextra, shdep, marker=1):
    eloc = pd.read_csv(efile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    topo = eloc.values[:,(1,3)]
    
#    topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    
    e_bot = np.floor(np.min(topo[:,1]))
    xL = np.floor(np.min(topo[:,0]))
    xR = np.ceil(np.max(topo[:,0]))

    L_ext = np.array([[xL-xextra, topo[0,1]]])
    R_ext = np.array([[xR+xextra, topo[topo.shape[0]-1,1]]])
    topo_ext = np.concatenate((L_ext, topo, R_ext))
    topo_bot = np.copy(np.flipud(topo_ext))
    topo_bot[:,1] = topo_bot[:,1]-shdep
    combo = np.concatenate((topo_ext, topo_bot))

    lpoly = mt.createPolygon(combo, isClosed=True, marker=marker)

    return topo, lpoly
    
def createFPpoly(efile, fpthx, fplim, shdep=0):

    fpdep = shdep + fpthx

    eloc = pd.read_csv(efile,sep='\t',header=None)
    nelec = eloc.values.shape[0]
    topo = eloc.values[:,(1,3)]
    
    #topo = np.genfromtxt(efile ,delimiter=',',skip_header=True)
    topo[:,1] = topo[:,1] - shdep
    
    e_bot = np.floor(np.min(topo[:,1]))
    iL = np.flipud(np.argwhere(topo[:,0] < fplim[0]))[0][0]
    iR = np.argwhere(topo[:,0] > fplim[1])[0][0]

    topo_fp = topo[iL:iR,:]
    topo_bot = np.copy(np.flipud(topo_fp))
    topo_fp[:,1] = topo_fp[:,1]
    topo_bot[:,1] = topo_bot[:,1]-(fpdep-shdep)
    combo = np.concatenate((topo_fp, topo_bot))

    med = np.int(topo_fp.shape[0]/2) #Setting position for the marker
    markerpos = [topo_bot[med,0],topo_bot[med,1]+(fpdep-shdep)/2]
    
#    lpoly = mt.createPolygon(combo, isClosed=True, markerPosition=markerpos, marker=4)

    return combo, markerpos

def createFloodplain(verts, markerpos, marker=4):
    return mt.createPolygon(verts, isClosed=True, markerPosition=markerpos, marker=marker)

# --------- Create Fault Polygon ------------
def createFault(topo, fpverts, xpos, dip, H, dep, shdep=0, xtra=500, marker=3):
    fpverts_new = np.copy(fpverts)
    
    topo2 = np.concatenate((np.array([[topo[0,0]-xtra, topo[0,1]]]), topo, np.array([[topo[topo.shape[0]-1,0]+xtra,topo[topo.shape[0]-1,1]]])))
    topo2[:,1] = topo2[:,1]-shdep
    
    zbot = np.floor(np.min(topo2[:,1])) - dep+shdep

    Z = interpolate.interp1d(topo2[:,0],topo2[:,1])
    zx = Z(xpos)
    Dz = zx-zbot

    Dx = Dz/np.tan(np.deg2rad(dip))
    xbot = xpos + Dx

    Dxbot = H/2/np.sin(np.deg2rad(dip))

    LR = np.array((xbot+Dxbot,zbot))[np.newaxis,:]
    LL = np.array((xbot-Dxbot,zbot))[np.newaxis,:]

    zfR = zbot + (LR[0,0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffR = interpolate.interp1d(zfR-topo2[:,1],topo2[:,0])
    UR = np.array((float(diffR(0)), float(Z(diffR(0))) ))[np.newaxis,:]
    if UR[0,0] > np.min(fpverts[:,0]):
        midind = np.int(fpverts.shape[0]/2)
        botverts = np.flipud(fpverts[midind:,:])
        
        zfR = zbot + (LR[0,0]-botverts[:,0])*np.tan(np.deg2rad(dip))
        diffR = interpolate.interp1d(zfR-botverts[:,1],botverts[:,0])
        
        Zfp = interpolate.interp1d(botverts[:,0],botverts[:,1])
        
        UR = np.array((float(diffR(0)), float(Zfp(diffR(0))) ))
        botverts_new_L = np.flipud(botverts[botverts[:,0]<UR[0],:])
        botverts_new_R = np.flipud(botverts[botverts[:,0]>UR[0],:])
        botverts_new = np.concatenate((botverts_new_R, UR[np.newaxis,:], botverts_new_L))
        fpverts_new = np.concatenate((fpverts[:midind], botverts_new))
        
        inds_left_of_intersect = botverts[:,0] < UR[0]
        
        UR_fp = np.flipud(botverts[inds_left_of_intersect,:])
        UR_total = np.concatenate((UR[np.newaxis,:], UR_fp, fpverts[0,:][np.newaxis,:]))
        UR = np.flipud(UR_total)
        #sys.exit("Cannot handle fracture zone intersecting floodplain yet")

    zfL = zbot + (LL[0,0]-topo2[:,0])*np.tan(np.deg2rad(dip))
    diffL = interpolate.interp1d(zfL-topo2[:,1],topo2[:,0])
    UL = np.array((float(diffL(0)), float(Z(diffL(0))) ))[np.newaxis,:]    
    
    idxabove = (topo2[topo2[:,0] > UL[0,0],0]).argmin() + (topo2[topo2[:,0] < UL[0,0],0]).shape[0]
    idxbelow = (topo2[topo2[:,0] < UR[0,0],0]).argmax()

    middles = np.array([(topo[j,0],topo[j,1] - shdep) for j in range(idxabove-1,idxbelow)])

    #verts = [LL, UL] + middles + [UR,LR]
    verts = np.concatenate((LL,UL,middles,UR,LR))
    fpoly = mt.createPolygon(verts, isClosed=True, addNodes=0, marker=marker)

    return fpoly, fpverts_new


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
    def __init__(self, dip=60, H=100, xpos=250, rho_fault=30, rho_back=50, shdep=1, rho_sh=40, fpthx=3, rho_fp=80, fplim=[450,650], efile=None, sfile=None, dep=200, xtra=500, Q=20, outfile=None, heterogeneity=False, area=None):
        self.dep = dep
        self.dip = dip
        self.H = H
        self.xpos = xpos
        self.rho_fault = rho_fault
        self.rho_back = rho_back
        self.shdep = shdep
        self.rho_sh = rho_sh
        self.fpthx = fpthx
        self.rho_fp = rho_fp
        self.Q = Q
        self.outfile = outfile
        self.xtra = xtra
        self.area = area
        self.fplim = fplim
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
            self.rhomap = [[0, rho_back],
                          [1, rho_fault],
                          [2, rho_back],
                          [3, rho_sh],
                          [4, rho_fp]]    
                
        

    def create_fault(self):
        topo, _ = createTopo(self.efile,self.xtra,self.dep)
        return createFault(topo, self.xpos, self.dip, self.H, self.dep, xtra=self.xtra)

    def create_geom(self):
    
        _, bpoly = createBackground(self.efile, self.xtra, self.dep, self.shdep, marker=0)      
        
        self.topo, lpoly = createLayer(self.efile, self.xtra, self.shdep, marker=3)

        fpverts, markerpos = createFPpoly(self.efile, self.fpthx, self.fplim, shdep=self.shdep) 
        
        fpoly,fpverts_new = createFault(self.topo, fpverts, self.xpos, self.dip, self.H, self.dep, shdep=self.shdep, xtra=self.xtra, marker=1)
        fppoly = createFloodplain(fpverts_new, markerpos, marker=4)
        return bpoly+lpoly+fpoly+fppoly
        
#        _, bpoly = createBackground(self.efile,self.xtra, self.dep, self.shdep)
#        topo, lpoly = createLayer(self.efile, self.xtra, self.shdep)
#        fpoly = createFault(topo, self.xpos, self.dip, self.H, self.dep, self.shdep, xtra=self.xtra)
#        self.topo = topo
#        return bpoly+lpoly+fpoly

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






