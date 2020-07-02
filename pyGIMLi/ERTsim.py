# ---------------------- Parameters ----------------------
dep = 200

dip = 45 # from 0-180
H = 100
xpos = 250

rho_fault = 30
rho_back = 50
rho_air = None

efile = '/home/ammilten/Documents/LBNL/ERT/PH-2018-eloc.txt'
sfile = '/home/ammilten/Documents/LBNL/ERT/PH-2018-srv.txt'

# -----------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate 

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert


# ----------Create Topo Polygon ----------------------
print('Creating topo...')

eloc = pd.read_csv(efile,sep='\t',header=None)
srv = pd.read_csv(sfile,sep='\t',header=None)

e_top = np.ceil(np.max(eloc.values[:,3]))
e_bot = np.floor(np.min(eloc.values[:,3]))

topo = eloc.values[:,(1,3)]
L_ext = np.array([[-500, eloc.values[0,3]]])
R_ext = np.array([[eloc.values[320,1]+500, eloc.values[320,3]]])
bots = np.array([[R_ext[0,0], e_bot-dep],[L_ext[0,0],e_bot-dep]])
combo = np.concatenate((topo,R_ext,bots,L_ext))

tpoly = mt.createPolygon(combo, isClosed=True, marker=2)

# --------- Create Fault Polygon ------------
print('Creating fault...')
Z = interpolate.interp1d(topo[:,0],topo[:,1])
zx = Z(xpos)
Dz = zx-e_bot+dep

b = Dz/np.tan(np.deg2rad(dip))
xbot = xpos + b

a = H/2/np.sin(np.deg2rad(dip))

LR = (xbot+a,e_bot-dep)
LL = (xbot-a,e_bot-dep)

X = interpolate.interp1d(topo[:,0],topo[:,0]-(topo[:,1]-e_bot+dep)/np.tan(np.deg2rad(dip)))
UR = (X(xbot+a),Z(X(xbot+a)))
UL = (X(xbot-a),Z(X(xbot-a)))

fpoly = mt.createPolygon([UL,UR,LR,LL], isClosed=True, addNodes=3, marker=1)

# ---------- Create mesh ---------------------
print('Creating mesh...')
geom = tpoly + fpoly

scheme = pg.DataContainerERT()

pos = eloc.values[:,(1,3)]
scheme.setSensorPositions(pos)

for i, elec in enumerate("abmn"):
    scheme[elec] = srv.values[:,i+1].astype(np.int)-1


for p in scheme.sensors():
    geom.createNode(p)
    geom.createNode(p - [0, 0.1])

# Create a mesh for the finite element modelling with appropriate mesh quality.
mesh = mt.createMesh(geom, quality=34)

# Create a map to set resistivity values in the appropriate regions
# [[regionNumber, resistivity], [regionNumber, resistivity], [...]
rhomap = [[0, rho_back],
          [1, rho_fault],
          [2, rho_back]]

# -------------- Simulate --------------------
print('Simulating ERT...')
data = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1, noiseAbs=1e-6, seed=1337, calcOnly=True)

data.set('r', data('u') / data('i'))
data.set('rhoa', data('k') * data('u') / data('i'))

pg.info('Simulated data', data)
pg.info('The data contains:', data.dataMap().keys())

pg.info('Simulated rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
pg.info('Selected data noise %(min/max)', min(data['err'])*100, max(data['err'])*100)

# --------- Plot ------------------------
print('Plotting...')

C = (srv.values[:,4]-srv.values[:,2])/2
M = C + srv.values[:,2]
R = np.array([data('r')[i] for i in range(data('r').size())])

fig, ax = plt.subplots(nrows=2, ncols=1)
plt.subplot(211)
plt.scatter(M,-C,5,c=np.log(R))
plt.title('Synthetic')
plt.colorbar(label='Resistance (Ohm)')
plt.clim(-7,1)

plt.subplot(212)
plt.scatter(M,-C,5,c=np.log(srv.values[:,5]))
plt.title('True')
plt.colorbar(label='Resistance (Ohm)')
plt.clim(-7,1)

fig.tight_layout()

