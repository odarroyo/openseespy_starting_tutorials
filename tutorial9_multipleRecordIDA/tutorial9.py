# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:24:59 2022

@author: orlandoaram
"""

# hay que llamar al opensees
from openseespy.opensees import *
# hay que llamar al opsvis para poder plotear la visualización
import opsvis as opsv
# esta es una librería estándar para plotear otras cosasy poder crear figuras.
import matplotlib.pyplot as plt

import opseestools.analisis as an

import numpy as np

import multiprocessing

from joblib import Parallel, delayed

import os as os

import time

import warnings
warnings.filterwarnings("ignore")

wipe()
 
# Importar análisis de gravedad
# Esta rutina hace un IDA de un registro

#%% INFORMACION DE ENTRADA
# ======================================

# records es una lista de los registros de los nombres de los archivos txt de los sismos. Deben estar en la misma carpeta de este archivos
# SpectrumFactor incluye la lista de los factores para escalar los registros a un espectro determinado
# NSteps inclute el número de datos de cada terremoto.
# DTs incluye los incrementos de tiempo de cada registros
# nodes incluye el nodo del techo
# SFactor son los factores escalares (tomando como 1.0 el sismo escalado al espectro según SpectrumFactor)

records= ["GM01.txt", "GM02.txt", "GM03.txt", "GM04.txt", "GM05.txt", "GM06.txt", "GM07.txt", "GM08.txt", "GM09.txt", "GM10.txt", "GM11.txt", "GM12.txt", "GM13.txt", "GM14.txt", "GM15.txt",
         "GM16.txt", "GM17.txt", "GM18.txt", "GM19.txt", "GM20.txt", "GM21.txt", "GM22.txt", "GM23.txt", "GM24.txt", "GM25.txt", "GM26.txt", "GM27.txt", "GM28.txt", "GM29.txt", "GM30.txt",
          "GM31.txt", "GM32.txt", "GM33.txt", "GM34.txt", "GM35.txt", "GM36.txt", "GM37.txt", "GM38.txt", "GM39.txt", "GM40.txt", "GM41.txt", "GM42.txt", "GM43.txt", "GM44.txt"]
SpectrumFactor= [0.647426972,0.498726473,0.503999614,0.478161366,0.589546608,0.477704312,1.682467501,1.20620368,1.029770012,0.73124417,1.249838191,1.064764959,0.636741252,0.577270921,0.640404291,0.818394714,0.815874762,0.867515359,1.84137413,2.574316418,0.991379457,2.2349985,1.584584522,0.60577121,0.547166565,0.413121928,0.913989311,1.47252665,0.745693994,1.042030541,0.8533576,1.826641759,1.10171508,1.074092626,0.657668909,0.748815552,0.818734703,0.617128352,0.609401699,0.687171747,2.453660823,3.013166191,1.884958295,0.543088642]
Nsteps= [3000, 3000, 2000, 2000, 5590, 5590, 4535, 4535, 9995, 9995, 7810, 7810, 4100, 4100, 4100, 4100, 5440, 5440, 6000, 6000, 2200, 2200, 11190, 11190, 7995, 7995, 7990, 7990, 2680, 2300, 8000, 8000, 2230, 2230, 1800, 1800, 18000, 18000, 18000, 18000, 2800, 2800, 7270, 7270]
DTs= [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.005,0.005,0.01,0.01,0.01,0.01,0.005,0.005,0.05,0.05,0.02,0.02,0.0025,0.0025,0.005,0.005,0.005,0.005,0.02,0.02,0.005,0.005,0.01,0.01,0.02,0.02,0.005,0.005,0.005,0.005,0.01,0.01,0.005,0.005]
GMcode = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]
nodes = [2016]
SFactor = [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0]
#%% 

nrecords = len(DTs)
dtecho = []
dtechomax = []

#%% Creación de la función que realiza el análisis 
# ind es el número del registro, fact el factor escalar

# En donde dice DEFINICION DEL MODELO, se deberá definir el modelo HASTA las cargas de gravedad
logFile('log.out','-noEcho')
def rundyn(fact,ind):
    
    rec = str(ind+1) # linea que crea el indice para nombrar el registro
    factor = 9.81*fact
    # CARGAR MODELO DE GRAVEDAD
    # ================================================
    nombre = str(int(factor/9.81*100))
    
    # DEFINICION DEL MODELO
    # ============================================
    
    # execfile('parametricgenerationMVLEM_IDA.py') # Aquí debe estar el modelo que aplica las cargas de gravedad
    wipe()
    #
    # wipe() mientras no funcione bien el script es mejor tener el wipe antes
    # creación del modelo
    model('basic','-ndm',2,'-ndf',3)

    # %% INFORMACIÓN DE ENTRADA

    xloc = [0.0,9.15]
    yloc = [0.0,2.6,5.2,7.8,10.4,13.0,15.6,18.2,20.8,23.4,26.0,28.6,31.2,33.8,36.4,39.0,41.6]

    # Aquí se debe colocar las secciones de las vigas y columnas en alturas
    # colsecs = [sec30x30,sec30x30,sec30x30,sec30x30] # secciones de las columnas en altura
    # vigsecs = [sec30x40,sec30x40,sec30x40,sec30x40] # secciones de las vigas en altura

    diafragma = 1 # colocar 1 si se desea diafragma en cada piso

    pushlimit = 0.015 # limite del pushover
    pushtype = 2 # 1 para triangular, 2 para uniforme, 3 para proporcional al modo
    modepush = 1 # modo para las cargas del pushover

    w = 797.5*2/9.15 # carga uniforme sobre las vigas

    # muro 1 y 2
    rhobordeI = [0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012] # cuantia borde izquierdo en altura
    rhobordeD = [0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012] # cuantia borde derecho en altura
    rhoalma = [0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012,0.012] # cuantia alma en altura
    width = [0.3,1.175,1.175,1.175,1.175,1.175,1.175,0.3] # ancho de las macrofibras
    nfib = len(width)
    t1 = [7.5,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
    t2 = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,7.5]
    ta = [t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1,t1]
    tb = [t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2,t2]
    muro1 = [rhobordeI,rhobordeD,rhoalma,width,ta]
    muro2 = [rhobordeI,rhobordeD,rhoalma,width,tb]

    muros = [muro2,muro1]

    cMVLEM = 0.4 # parámetro

    # vigas del modelo
    bv = 0.3 # base
    hv = 0.5 # altura
    Av = bv*hv # area
    Iv = bv*hv**3/12 # inercia
    Ev = 24000000.0 # modulo de elasticidad 

    #%% NODOS Y RESTRICCIONES
    # =========================================


    ny = len(yloc)
    nx = len(xloc)

    for i in range(nx):
        for j in range(ny):
            nnode = 1000*(i+1)+j
            node(nnode,xloc[i],yloc[j])

    # plt.figure()
    # opsv.plot_model()

    print('Nodos generados')

    # apoyos
    empotrado = [1,1,1]
    grado2 = [1,1,0]

    # para colocarlos todos al nivel 0.0
    fixY(0.0,*empotrado)

    print('Restricciones asignadas')

    # %% ASIGNACIÓN DE DIAFRAGMAS
    # =========================================

    if diafragma == 1:
        for j in range(1,ny):
            for i in range(1,nx):
                masternode = 1000 + j
                slavenode = 1000*(i+1) + j
                equalDOF(masternode,slavenode,1)
        print('Diafragmas asignados')


    # %% ASIGNACION DE MASAS
    # =========================================

    # masas
    mass(1001,81.3,81.3,0.0)
    mass(1002,81.3,81.3,0.0)
    mass(1003,81.3,81.3,0.0)
    mass(1004,81.3,81.3,0.0)
    mass(1005,81.3,81.3,0.0)
    mass(1006,81.3,81.3,0.0)
    mass(1007,81.3,81.3,0.0)
    mass(1008,81.3,81.3,0.0)
    mass(1009,81.3,81.3,0.0)
    mass(1010,81.3,81.3,0.0)
    mass(1011,81.3,81.3,0.0)
    mass(1012,81.3,81.3,0.0)
    mass(1013,81.3,81.3,0.0)
    mass(1014,81.3,81.3,0.0)
    mass(1015,81.3,81.3,0.0)
    mass(1016,81.3,81.3,0.0)

    mass(2001,81.3,81.3,0.0)
    mass(2002,81.3,81.3,0.0)
    mass(2003,81.3,81.3,0.0)
    mass(2004,81.3,81.3,0.0)
    mass(2005,81.3,81.3,0.0)
    mass(2006,81.3,81.3,0.0)
    mass(2007,81.3,81.3,0.0)
    mass(2008,81.3,81.3,0.0)
    mass(2009,81.3,81.3,0.0)
    mass(2010,81.3,81.3,0.0)
    mass(2011,81.3,81.3,0.0)
    mass(2012,81.3,81.3,0.0)
    mass(2013,81.3,81.3,0.0)
    mass(2014,81.3,81.3,0.0)
    mass(2015,81.3,81.3,0.0)
    mass(2016,81.3,81.3,0.0)

    #%%  MATERIALES Y TRANSFORMACIONES
    # =========================================

    # Definición de material 
    #uniaxialMaterial('Elastic',1,E)
    fc = 35000.0
    ec = 0.002
    E = 1000*4300*(fc/1000)**0.5
    fcu = 0.1*fc
    ecu = 0.006
     
    k=1.3
    fcc=fc*k
    ecc= 2*fcc/E
    fucc=0.2*fcc
    eucc=0.02
     
    Fy=420000.0
    Es=210000000.0
    uniaxialMaterial('Concrete01', 2, -fc, -ec, -fcu, -ecu)
    uniaxialMaterial('Concrete01', 1, -fcc, -ecc, -fucc, -eucc)
    uniaxialMaterial('Steel01', 4, Fy, Es, 0.01)
    ey = Fy/Es

    # Acero 

    fu = 630000.0
    eult = 0.1
    uniaxialMaterial('Hysteretic',3,Fy,ey,fu,eult,0.05*Fy,0.11,-Fy,-ey,-fu,-eult,-0.05*Fy,-0.11,1.0,1.0,0.0,0.0)

    # Definición de elementos

    # material del shear
    tagshear = 1000
    G = E*0.4
    uniaxialMaterial('Elastic',tagshear,G*1.5)

    # %% Transformación

    lineal = 1
    geomTransf('Linear',lineal)

    pdelta = 2
    geomTransf('PDelta',pdelta)

    cor = 3
    geomTransf('Corotational',cor)


    # %% ELEMENTOS VERTICALES (MUROS)
    # ========================

    #element('MVLEM', eleTag, Dens, *eleNodes, m, c, '-thick', *thick, '-width', *widths, '-rho', *rho, '-matConcrete', *matConcreteTags, '-matSteel', *matSteelTags, '-matShear', matShearTag)

    concTags = [2,2,2,2,2,2,2,2]
    steelTags = [3,3,3,3,3,3,3,3]
    shearTags = [tagshear,tagshear,tagshear,tagshear,tagshear,tagshear,tagshear,tagshear]

    for i in range(nx):
        # se toma el muro actual y luego se sacan de el las propiedades
        mur = muros[i] # muro actual
        rbi = mur[0] # cuantia elemento borde izquierdo
        rbd = mur[1] # cuantia elemento borde derecho
        ra = mur[2] # cuantia del alma
        widths = mur[3] # ancho de las macrofibras
        thickness = mur[4] # espesor del muro
        rhou = [] 
        walltags = []      
        for j in range(ny-1):
            thick = thickness[j]
            nodeI = 1000*(i+1) + j
            nodeJ = 1000*(i+1) + (j+1)
            eltag = 1000*(i+1) + j
            walltags.append(eltag)
            for k in range(nfib):
                if k == 0:
                    rhou.append(rbi[j])
                elif k == nfib-1:
                    rhou.append(rbd[j])
                else:
                    rhou.append(ra[j])
            element('MVLEM',eltag, 0.0, nodeI,nodeJ, nfib, cMVLEM, '-thick', *thick, '-width', *widths, '-rho', *rhou, '-matConcrete', *concTags, '-matSteel', *steelTags, '-matShear',*shearTags)            

    # plt.figure()        
    # opsv.plot_model() 
    print('Muros generados')

    # %% ELEMENTOS HORIZONTALES (VIGAS)
    # ========================

    tagvigas=[]
    for j in range(1,ny):
        for i in range(nx-1):
            nodeI = 1000*(i+1) + j
            nodeJ = 1000*(i+2) + j
            eltag = 100000*(i+1) + j
            tagvigas.append(eltag)
            element('elasticBeamColumn',eltag,nodeI,nodeJ,Av,Ev,Iv,lineal)
            #element('forceBeamColumn',eltag,nodeI,nodeJ ,lineal,vigsecs[j-1]) # aqui basta con cambiar sec30x30 por una que se obtenga de una lista con las secciones de cada piso

    # plt.figure()
    # opsv.plot_model()
    print('vigas generadas')

    # %% CARGA SOBRE LOS MUROS
    # =========================================
    timeSeries('Linear', 1)
    pattern('Plain',1,1)
    # eleLoad('-ele',*tagvigas,'-type','-beamUniform',-w)
    # eleLoad('-ele',*tagvigas,'-type','beamPoint',-50.0,0.6)

    load(1001,0.0,-797.5,0.0)
    load(1002,0.0,-797.5,0.0)
    load(1003,0.0,-797.5,0.0)
    load(1004,0.0,-797.5,0.0)
    load(1005,0.0,-797.5,0.0)
    load(1006,0.0,-797.5,0.0)
    load(1007,0.0,-797.5,0.0)
    load(1008,0.0,-797.5,0.0)
    load(1009,0.0,-797.5,0.0)
    load(1010,0.0,-797.5,0.0)
    load(1011,0.0,-797.5,0.0)
    load(1012,0.0,-797.5,0.0)
    load(1013,0.0,-797.5,0.0)
    load(1014,0.0,-797.5,0.0)
    load(1015,0.0,-797.5,0.0)
    load(1016,0.0,-797.5,0.0)

    load(2001,0.0,-797.5,0.0)
    load(2002,0.0,-797.5,0.0)
    load(2003,0.0,-797.5,0.0)
    load(2004,0.0,-797.5,0.0)
    load(2005,0.0,-797.5,0.0)
    load(2006,0.0,-797.5,0.0)
    load(2007,0.0,-797.5,0.0)
    load(2008,0.0,-797.5,0.0)
    load(2009,0.0,-797.5,0.0)
    load(2010,0.0,-797.5,0.0)
    load(2011,0.0,-797.5,0.0)
    load(2012,0.0,-797.5,0.0)
    load(2013,0.0,-797.5,0.0)
    load(2014,0.0,-797.5,0.0)
    load(2015,0.0,-797.5,0.0)
    load(2016,0.0,-797.5,0.0)

  
    #%% ANALISIS
    # ========================

    an.gravedad()
    loadConst('-time',0.0)
    # DEFINICIÓN DE RECORDERS
    # ================================================
    
    # para que grabe cada uno distinto se debe usar el indice que se creo en rec. Por ejemplo, para el nodo de techo:
    # recorder('Node','-file',nombre+'/roof'+rec+'.out','-time','-node',*nodes,'-dof',1,'disp') 
    
    # EJECUCION DEL ANÁLISIS DINÁMICO
    # ================================================
    
    factoran = SpectrumFactor[ind]*factor
    
    t,techo = an.dinamicoIDA2(records[ind], DTs[ind], Nsteps[ind], 0.04, factoran, 0.025, 2016, 1,[0,2])
    wipe()
    return ind,fact,t,techo
    


# %% Ejemplo 2

num_cores = multiprocessing.cpu_count() # esta linea identifica el número de nucleos totales del PC.
# En equipos con SMT identifica los núcleos físicos y lógicos. La recomendación si se va a seguir usando el PC es dejar dos núcleos físicos libres


stime = time.time()

# resultados devuelve de momento cuatro cosas. La primera es el indice del terremoto, la segunda es el factor escalar, la tercera el tiempo del registro y la cuarta es el desplazamiento de techo
resultados = Parallel(n_jobs=num_cores)(delayed(rundyn)(ff,pp) for ff in SFactor for pp in GMcode) # loop paralelo

etime = time.time()

ttotal = etime - stime

print('tiempo de ejecucion: ',ttotal,'segundos')
# %% Procesamiento

maxtecho2 = np.zeros(len(resultados)) # para guardar el máximo desplazamiento de techo
maxt2 = np.zeros(len(resultados)) # para guardar el máximo tiempo

# este ciclo calcula el máximo desplazamiento de techo y tiempo del registro                    
for indice,result in enumerate(resultados):
    dtecho = result[3]
    tiempo = result[2]
    maxtecho2[indice] = np.max(dtecho)
    maxt2[indice] = np.max(tiempo)

Sa = np.array(SFactor)*0.5625

# for k in range(len(Sa)):

for k in range(44):
    s = maxtecho2[k::44]
    plt.plot(s,Sa,'o-')

plt.xlim([0,2])    
plt.ylabel('Sa(g)')
plt.xlabel('Roof displacement (m)')

plt.show()

etime2 = time.time()

ttotal2 = etime2 - stime

