# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:59:44 2022

@author: Orlando
"""
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ANALISIS DE GRAVEDAD
# =============================
def gravedad():
    
# Create the system of equation, a sparse solver with partial pivoting
    system('BandGeneral')

# Create the constraint handler, the transformation method
    constraints('Plain')

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
    numberer('RCM')

# Create the convergence test, the norm of the residual with a tolerance of
# 1e-12 and a max number of iterations of 10
    test('NormDispIncr', 1.0e-12, 10, 3)

# Create the solution algorithm, a Newton-Raphson algorithm
    algorithm('Newton')

# Create the integration scheme, the LoadControl scheme using steps of 0.1
    integrator('LoadControl', 0.1)

# Create the analysis object
    analysis('Static')

    ok = analyze(10)
    
    # if ok != 0:
    #     print('Análisis de gravedad fallido')
    #     sys.exit()
    # else:
    #     print('Análisis de gravedad completado')
        


# ANALISIS PUSHOVER
# =============================
            
def pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    
    # Eds = np.zeros((nels, Nsteps+1, 3)) # para grabar las rotaciones de los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V

#%% ================================ PUSHOVER CON REMOVAL =================================

def pushover2R(Dmax,Dincr,IDctrlNode,IDctrlDOF,columns,beams,ele,der,nodes_control,elements,norm=[-1,1],Tol=1e-8):
    # Dmax: Deriva máxima de techo --> Altura del edificio*0.05
    # Dincr: Deriva incremental 
    # IDctrlNode: Nodo de control
    # IDctrlDOF: 1
    # columns: Tag de las columnas
    # beams: Tag de las vigas
    # ele: Tag elementos que componen el muro (puntales)
    # der: Deriva máxima de cada muro. Despues de ese valor se rompe
    # nodes_control: Nodos pushover
    # elements: Tag columnas del primer piso para obtener fuerzas locales
    
    maxNumIter = 10                                                             # Máximo número de iteraciones
    
    # ------------------- Configuración básica del análisis -------------------
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    # ----------------------- Otras opciones de análisis ----------------------   
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    # -------------------------- Rutina del análisis --------------------------
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    dmuro1 = [nodeDisp(1001,IDctrlDOF)]
    dmuro2 = [nodeDisp(1002,IDctrlDOF)]
    Vbasal = [getTime()]
    
    nels = len(elements)
    Eds = np.zeros((nels, Nsteps+1, 6)) 
                                   
    
    nnodos = len(nodes_control)
    f_puntalB = np.zeros((len(ele),Nsteps,6))
    f_puntalA = np.zeros((len(ele),Nsteps,6))
    flag = np.zeros((len(ele),Nsteps))
    node_disp = np.zeros((Nsteps + 1, nnodos))                                  # Para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1))                                  # Para grabar la deriva de entrepiso
    
    kfails = 1e9*np.ones(len(ele))                                              # Para identificar el punto donde el muro "falla" (tiempo)
        
    nodeI=[]
    nodeJ=[]
    for n in range(len(ele)):
        nodes= eleNodes(ele[n][1])
        nodeI.append(nodes[0])
        nodeJ.append(nodes[1])
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
        for i in range(len(ele)):
            f1,a,flag1 = removalTH2(nodeI[i],nodeJ[i],ele[i],der[i])
            f_puntalA[i,k,:] = f1
            f_puntalB[i,k,:] = a
            if flag1 == 1:
                kfails[i] = np.min((k,kfails[i]))
                
        for el_i, ele_tag in enumerate(elements):
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    unicos = np.unique(kfails)                                                  # Para encontrar los valores únicos donde se producen fallas de la mamposteria
    unicos2 = unicos[unicos<1e8]                                                # Para quitar los valores de 1e9 que hicimos para el artificio
    
    return techo, V, Eds, f_puntalA, f_puntalB, node_disp, drift, unicos2

#%%Pushover con removal, considernado rotaciones y deformaciones unitarias 
#En este es importante cambiar las fibras en donde se quiere guardar registro
def pushover2R_Rot_def(Dmax,Dincr,IDctrlNode,IDctrlDOF,columns,beams,ele,der,nodes_control,elements,id_s,id_c,norm=[-1,1],Tol=1e-8):
    # Dmax: Deriva máxima de techo --> Altura del edificio*0.05
    # Dincr: Deriva incremental 
    # IDctrlNode: Nodo de control
    # IDctrlDOF: 1
    # columns: Tag de las columnas
    # beams: Tag de las vigas
    # ele: Tag elementos que componen el muro (puntales)
    # der: Deriva máxima de cada muro. Despues de ese valor se rompe
    # nodes_control: Nodos pushover
    # elements: Tag columnas del primer piso para obtener fuerzas locales
    
    maxNumIter = 10                                                             # Máximo número de iteraciones
    
    # ------------------- Configuración básica del análisis -------------------
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    # ----------------------- Otras opciones de análisis ----------------------   
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    # -------------------------- Rutina del análisis --------------------------
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    dmuro1 = [nodeDisp(1001,IDctrlDOF)]
    dmuro2 = [nodeDisp(1002,IDctrlDOF)]
    Vbasal = [getTime()]
    
    nels = len(elements)
    ncols = len(columns)
    nbeams = len(beams)
    Eds = np.zeros((nels, Nsteps+1, 6)) 
    Prot_cols = np.zeros((ncols, Nsteps+1, 3)) # para grabar las rotaciones de los elementos
    Prot_beams = np.zeros((nbeams, Nsteps+1, 3)) # para grabar las rotaciones de los elementos
    fiber_c = np.zeros((ncols, Nsteps+1, 2)) 
    fiber_s = np.zeros((ncols, Nsteps+1, 2)) 
    
    nnodos = len(nodes_control)
    f_puntalB = np.zeros((len(ele),Nsteps,6))
    f_puntalA = np.zeros((len(ele),Nsteps,6))
    flag = np.zeros((len(ele),Nsteps))
    node_disp = np.zeros((Nsteps + 1, nnodos))                                  # Para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1))                                  # Para grabar la deriva de entrepiso
    
    kfails = 1e9*np.ones(len(ele))                                              # Para identificar el punto donde el muro "falla" (tiempo)
        
    nodeI=[]
    nodeJ=[]
    for n in range(len(ele)):
        nodes= eleNodes(ele[n][1])
        nodeI.append(nodes[0])
        nodeJ.append(nodes[1])
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
        for i in range(len(ele)):
            f1,a,flag1 = removalTH2(nodeI[i],nodeJ[i],ele[i],der[i])
            f_puntalA[i,k,:] = f1
            f_puntalB[i,k,:] = a
            if flag1 == 1:
                kfails[i] = np.min((k,kfails[i]))
                
        for el_i, ele_tag in enumerate(elements):
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
        for el_i, col_tag in enumerate(columns):
           
           fiber_s[el_i , k+1, :] = eleResponse(col_tag,'section',1,'fiber',0.19,0.0,id_s,'stressStrain')
           fiber_c[el_i , k+1, :] = eleResponse(col_tag,'section',1,'fiber',-0.23,0.0,id_c,'stressStrain')
        
           Prot_cols[el_i , k+1, :] = [eleResponse(col_tag,'plasticDeformation')[0],
                                eleResponse(col_tag,'plasticDeformation')[1],
                                eleResponse(col_tag,'plasticDeformation')[2]]
           
        for el_i, b_tag in enumerate(beams):
        
           Prot_beams[el_i , k+1, :] = [eleResponse(b_tag,'plasticDeformation')[0],
                                eleResponse(b_tag,'plasticDeformation')[1],
                                eleResponse(b_tag,'plasticDeformation')[2]]
           
           
            
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    unicos = np.unique(kfails)                                                  # Para encontrar los valores únicos donde se producen fallas de la mamposteria
    unicos2 = unicos[unicos<1e8]                                                # Para quitar los valores de 1e9 que hicimos para el artificio
    
    return techo, V, Eds, f_puntalA, f_puntalB, node_disp, drift, unicos2, Prot_cols,Prot_beams, fiber_s, fiber_c


def pushover2Rot(Dmax,Dincr,IDctrlNode,IDctrlDOF,elements,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    nels = len(elements)
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    Prot = np.zeros((nels, Nsteps+1, 3)) # para grabar las rotaciones de los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for el_i, ele_tag in enumerate(elements):
            
            Prot[el_i , k+1, :] = [eleResponse(ele_tag,'plasticDeformation')[0],
                                  eleResponse(ele_tag,'plasticDeformation')[1],
                                  eleResponse(ele_tag,'plasticDeformation')[2]]
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V, Prot

def pushover2D(Dmax,Dincr,IDctrlNode,IDctrlDOF,nodes_control,norm=[-1,1],Tol=1e-8):
    '''Hace un pushover de la estructura extrayendo las derivas en los nodos indicados \n
    Los argumentos de entrada son: \n
    Dmax: desplazamiento objetivo \n
    Dincr: incremento de desplazamiento \n
    IDctrlNode: nodo de control del pushover\n
    IDctrlDOF: grado de libertad del pushover\n
    nodes_control: de cada piso para grabar las derivas\n
    norm: recibe un vector con la altura y el peso para normalizar\n
    Tol: tolerancia del análisis'''
    # creación del recorder de techo y definición de la tolerancia
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
         
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    nnodos = len(nodes_control)
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
      
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V, drift


def pushover2C(displ,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-4):
    '''Hace un Pushover Cíclico de una estructura. Recibe los siguientes argumentos: \n
        displ: desplazamientos partiendo de cero a los cuales deberá llegar el Pushover \n
        Dincr: incremento de desplazamiento \n
        IDctrlNode: nodo de control de la estructura para monitorear el desplazamiento \n
        IDctrlDOF: grado de libertad del nodo donde se controlará el desplazamiento \n
        norm: vector que recibe la altura y el peso de la estructura. Por defecto no normaliza \n
        Tol: tolerancia del análisis. Por defecto está en 1e-4 porque utiliza NormUnbalance
    '''
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('NormUnbalance', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    currdisp = 0 # LINEA NUEVA
    for dis in displ: # LINEA NUEVA
        Nsteps =  int(np.abs(dis - currdisp) / Dincr) # LINEA NUEVA PORQUE TOCO EDITAR LA OTRA
        if dis > 0: # LINEA NUEVA
            integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr) # LINEA NUEVA
        else: # LINEA NUEVA
            integrator('DisplacementControl', IDctrlNode, IDctrlDOF, -Dincr) # LINEA NUEVA
        for k in range(Nsteps):
            ok = analyze(1)
            # ok2 = ok;
            # En caso de no converger en un paso entra al condicional que sigue
            if ok != 0:
                print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
                for j in algoritmo:
                    if j < 4:
                        algorithm(algoritmo[j], '-initial')
        
                    else:
                        algorithm(algoritmo[j])
                    
                    # el test se hace 50 veces más
                    test('NormUnbalance', Tol, maxNumIter*50)
                    ok = analyze(1)
                    if ok == 0:
                        # si converge vuelve a las opciones iniciales de análisi
                        test('NormUnbalance', Tol, maxNumIter)
                        algorithm('Newton')
                        break
                        
            if ok != 0:
                print('Pushover analisis fallido')
                print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
                break
        
            
            dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
            Vbasal.append(getTime())
        currdisp = nodeDisp(IDctrlNode,IDctrlDOF) # LINEA NUEVA
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V


def pushover2T(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    fibras1 = [0]*8
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
       
        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER

def pushover3T(Dmax,Dincr,IDctrlNode,IDctrlDOF,elements,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    
    
    nels = len(elements)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    # Strains = np.zeros((Nsteps+1, 8, nels)) # # para grabar las deformaciones de los muros en las 8 fibras que tienen los elementos
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        
        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
        
        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER, Eds, Strains, cStress, sStress



def pushover3Tn(Dmax,Dincr,IDctrlNode,IDctrlDOF,elements,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    
    
    nels = len(elements)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    # Strains = np.zeros((Nsteps+1, 8, nels)) # # para grabar las deformaciones de los muros en las 8 fibras que tienen los elementos
    Strains = np.zeros((nels, Nsteps+1, 14))
    cStress = np.zeros((nels, Nsteps+1, 14)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 14)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        
        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7],
                                 eleResponse(ele_tag,'Fiber_Strain')[8],
                                 eleResponse(ele_tag,'Fiber_Strain')[9],
                                 eleResponse(ele_tag,'Fiber_Strain')[10],
                                 eleResponse(ele_tag,'Fiber_Strain')[11],
                                 eleResponse(ele_tag,'Fiber_Strain')[12],
                                 eleResponse(ele_tag,'Fiber_Strain')[13]]
                                 
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[8],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[9],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[10],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[11],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[12],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[13]]
                                 
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[8],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[9],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[10],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[11],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[12],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[13]]


        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER, Eds, Strains, cStress, sStress


# ANALISIS DINAMICO
# =============================   

# dinamico es el más sencillo de todos, corre un terremoto creando un recorder para el techo.
# dinamicoIDA crea el recorder en función del factor escalar.
# dinamicoAnim es dinamico pero guarda la información para animar el registro
# dinamicoIDA2 está modificado para ser utilizado cuando se desee correr en paralelo los cálculos. Devuelve solo el desplazamiento de techo
# dinamicoIDA3 PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS. SOLO PUEDEN SER LOS MUROS DE MOMENTO. También extrae desplazamientos de los nodos
# dinamicoIDA4 PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS. SOLO PUEDEN SER LOS MUROS DE MOMENTO. También extrae desplazamientos de los nodos, aceleraciones, derivas, velocidades, esfuerzos en concreto, acero y deformaciones unitarias de cada muro indicado en elements
# dinamicoIDA5 es lo mismo que IDA4, pero en lugar del nombre del registro, recibe una lista con las aceleraciones.

def dinamicoIDA2(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo


def dinamicoIDA4(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,Strains,cStress,sStress,node_disp,node_vel,node_acel,drift


def dinamicoIDA4T(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    # IGUAL A dinamicoIDA4T
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    Tf2 = [2*np.pi/w1]
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        # eigvalF = eigen(nmodes)
        
        # eigF = eigvalF[modes[0]]
        
        # wF = eigF**0.5
        # Tf= 2*np.pi/wF
        # Tf2.append(Tf)
        
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    eigvalF = eigen(nmodes)
    
    eigF = eigvalF[modes[0]]
    
    wF = eigF**0.5
    Tf= 2*np.pi/wF
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,node_disp,node_vel,node_acel,drift,Tf


def dinamicoIDA4P(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-4):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('NormUnbalance', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Prot = np.zeros((nels, Nsteps+1, 3)) # para grabar las rotaciones de los elementos
    
    
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormUnbalance', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormUnbalance', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
                      
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
        
            
            Prot[el_i , k+1, :] = [eleResponse(ele_tag,'plasticDeformation')[0],
                                  eleResponse(ele_tag,'plasticDeformation')[1],
                                  eleResponse(ele_tag,'plasticDeformation')[2]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
      
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    
    return tiempo,techo,Eds,node_disp,node_vel,node_acel,drift

def dinamicoIDA5(acceleration,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # acceleration es la lista de aceleraciones del registro. Se termina multiplicando por fact
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-values',*acceleration,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,Strains,cStress,sStress,node_disp,node_vel,node_acel,drift


#%% ========================== DINÁMICO CON REMOVAL ===========================
def dinamicoIDA4R(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,columns, beams, ele,der,elements,nodes_control,id_s,id_c,modes = [0,2],Kswitch = 1,Tol=1e-4):
    # recordName es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis (residual=nPts agergarle 2 seg dividos por el dt del registro)
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
    # columns es el tag de las columnas
    # beams es el tag de las vigas
    # ele es el tag de los puntales que definen el muro
    # der es la deriva en la que el muro falla y se rompe
    # ** elements son los elementos a los que se les calcula la fuerza **
    # nodes_control son los nodos de control
    
    maxNumIter = 10
    bandera = [0]*len(ele)
    
    # ------------------------- Creación del pattern --------------------------
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)

    # -------------------------------- Damping --------------------------------
    # Kswitch = 1 calcula amortiguamiento proporcional a la rigidez inicial, si no es 1 lo calcula 
    # proporcional a la rigidez tangencial
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
        
    # ------------------- Configuración básica del análisis -------------------
    
    wipeAnalysis()
    constraints('Plain')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # ---------------------- Otras opciones de análisis -----------------------
    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # -------------------------- Rutina del análisis --------------------------
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    ncol = len(columns)
    nnodos = len(nodes_control)
    f_puntalA = np.zeros((len(ele),Nsteps,6))                                   # Para grabar fuerzas en el puntal A
    f_puntalB = np.zeros((len(ele),Nsteps,6))                                   # Para grabar fuerzas en el puntal B
    Eds = np.zeros((nels, Nsteps+1, 6))                                         # Para grabar las fuerzas de los elementos
 
    node_disp = np.zeros((Nsteps + 1, nnodos))                                  # Para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos))                                   # Para grabar la velocidad de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos))                                  # Para grabar la aceleración de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1))                                  # Para grabar la deriva de entrepiso
    
    kfails = 1e9*np.ones(len(ele))                                              # Para identificar el punto donde hay fallas
  
    nodeI=[]
    nodeJ=[]
    for n in range(len(ele)):
        nodes= eleNodes(ele[n][1])
        nodeI.append(nodes[0])
        nodeJ.append(nodes[1])
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
        
        for el_i, ele_tag in enumerate(elements):
            
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]

            
        for i in range(len(ele)):
            if bandera[i] != 1:
                f1,a,flag1 = removalTH2(nodeI[i],nodeJ[i],ele[i],der[i])
                f_puntalA[i,k,:] = f1
                f_puntalB[i,k,:] = a
                bandera[i] = flag1
                if flag1 == 1:
                    kfails[i] = np.min((k,kfails[i]))
            else:
                f_puntalA[i,k,:] = [0]*6
                f_puntalB[i,k,:] = [0]*6
                # print('muerto el primero')
                if flag1 == 1:
                    kfails[i] = np.min((k,kfails[i]))
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
    
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    
    unicos = np.unique(kfails)                                                  # Para encontrar los valores únicos donde se producen fallas de la mamposteria
    unicos2 = unicos[unicos<1e8]                                                # Para quitar los valores de 1e9 que hicimos para el artificio
    
    # retornar tiempo, desplazamiento de recho, Fuerzas globales,fuerza en el 
    # puntal A, fuerza en el puntal B, Desplazamiento, velocidad, aceleración, 
    # deriva, instante de tiempo en el que el muro falla
    
    return tiempo,techo,Eds,f_puntalA,f_puntalB,node_disp,node_vel,node_acel,drift,unicos2, Prot, fiber_s, fiber_c

#%% ============================ FUNCIONES REMOVAL ============================

def removal(nodeI,nodeJ,ele,der):
    fuerza=0
    d1= nodeDisp(nodeI,1)
    d2= nodeDisp(nodeJ,1)
    y1= nodeCoord(nodeI)[1]
    y2= nodeCoord(nodeJ)[1]
    deriva= (d2-d1)/(y2-y1)
    if abs(deriva) > der:
        remove('ele',ele)
    else:
        fuerza = eleForce(ele,1)
    return fuerza

def removalTH(nodeI,nodeJ,ele,der):
    fuerza = 0
    flag = 0
    d1= nodeDisp(nodeI,1)
    d2= nodeDisp(nodeJ,1)
    y1= nodeCoord(nodeI)[1]
    y2= nodeCoord(nodeJ)[1]
    deriva= (d2-d1)/(y2-y1)
    if abs(deriva) > der:
        remove('ele',ele)
        flag = 1
    else:
        fuerzax = eleForce(ele,1,2,3)
       
    return fuerzax,flag
            
def removalTH2(nodeI,nodeJ,ele,der):
    fuerza1 = [0]*6
    fuerza2 = [0]*6
    flag = 0
    d1= nodeDisp(nodeI,1)
    d2= nodeDisp(nodeJ,1)
    y1= nodeCoord(nodeI)[1]
    y2= nodeCoord(nodeJ)[1]
    deriva= (d2-d1)/(y2-y1)
    if abs(deriva) > der:
        remove('ele',ele[0])
        remove('ele',ele[1])
        flag = 1
    else:
        fuerza1 = eleForce(ele[0])
        fuerza2 = eleForce(ele[1])

    return fuerza1, fuerza2, flag