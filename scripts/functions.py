#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 02:12:17 2022

@author: ginapantano
"""
import numpy as np
from collections import Counter
import pandas as pd
import os
import shutil


######  CREATE SUPERCELL INTERFACES & ASSOCIATED POSCAR FILES  #######
# =============================================================================
def supercell_interface(inf, pts, p, xl, xr, yl, yr, zl, zr, atoms, atom_labels, lattice, points, cgroup, POSCAR, rel):

    #Number of Cells in Each Direction for L and R Side of Interface.
    xl_size = xl
    xr_size = xr
    yl_size = yl
    yr_size = yr
    zl_size = zl
    zr_size = zr

    #Total Number of Cells in Each Direction
    xs = xl_size + xr_size
    ys = yl_size + yr_size
    zs = zl_size + zr_size

    #Number of Each Atom in Supercell 
    #Note: Array n will be in order atoms are ordered in input file from first to last coordinate
    #Example: [Sn, Sn, Co, Co, Co, Se, S] -> [2, 3, 1, 1] * (num of cells in x * in y * in z)
    n = np.array([])
    num_distinct_atoms = (Counter(atoms).keys()) #list of distinct atoms
    atom = atoms.tolist()
    for i in num_distinct_atoms:
        num = xs*ys*zs*atom.count(i)
        n = np.hstack([n, np.array([num])])
    
    #Atom Element Labels
    labels = str(atom_labels).replace('[', '').replace(']','').replace('\'','').replace(',',' ')

    #Changing primitive vectors to dimensions of supercell
    new_lat = np.array([[xs*lattice[0,0], ys*lattice[0,1], zs*lattice[0,2]],
                        [xs*lattice[1,0], ys*lattice[1,1], zs*lattice[1,2]], 
                        [xs*lattice[2,0], ys*lattice[2,1], zs*lattice[2,2]]])
    
    #Original Coordinates with Labels
    cell = np.hstack([points,atoms.reshape(len(atoms),1)])
 
    #Atomic Coordiantes in Left Cells 
    npl = np.array([]).reshape(0,4)

    nrel = 2 # num of atoms that can relax in nrel unit cells from the interface
    relce=np.array([]).reshape(0,3) #Array of labels for each atom in supercell. 
    trel=np.full((len(atoms), 3), 'T')  # T = atom is within unit cells allowed to relax
    frel=np.full((len(atoms), 3), 'F')  # F = atom is not within unit cells allowed to relax
    dis1=np.array([]) #coordinate vector for atom closest on left side of interface
    dis2=np.array([]) #coordinate vector for atom closest on right side of interface
    
    for i in range(xl_size):
        n_l = np.copy(cell)
        n_l[:,0]= (n_l[:,0]+i)/xs
        npl = np.vstack([npl,n_l])
        
        #Labels atoms in cells closet to interface allowed to relax as T
        if xl_size-i < nrel+1:
            relce = np.vstack([relce,trel])
        #Otherwise, label atoms as F
        else:
            relce = np.vstack([relce,frel])
        #Determines coordinate vector of atom closet to interface on left side
        if i==xl_size-1:
            rr=np.sqrt(np.sum(np.square(n_l),axis=1))
            ff=np.where(rr == max(rr)) 
            dis1=n_l[ff]
           
 
    #Sorted Transformed Coordinates with Labels
    sym_cell = cgroup[inf]

    #Atomic Coordinates in Right Cells
    npr = np.array([]).reshape(0,4)
    
    for i in range(xr_size):
        n_r = np.copy(sym_cell)
        n_r[:,0] = (n_r[:,0]+i+xl_size)/xs
        npr = np.vstack([npr,n_r])
        
        #Labels atoms in cells closet to interface allowed to relax as T
        if i < nrel:
            relce = np.vstack([relce,trel])
        #Otherwise, label atoms as F
        else:
            relce = np.vstack([relce,frel])
        #Determines coordinate vector of atom closet to interface on right side
        if i==1:
            rr=np.sqrt(np.sum(np.square(n_r),axis=1))
            ff=np.where(rr == min(rr))
            dis2=n_r[ff]
    
    x_max = np.max(npl[:,0])          
    x_min = np.min(npr[:,0])
    
    #Distance along x-axis between two closest atoms at interface
    x_dis = x_min - x_max
    #Total distance between two closet atoms at interface
    totd=np.sqrt(np.sum(np.square(dis2-dis1),axis=1))    

    dists = [x_dis, totd[0]]

    #Creates Supercell
    scell = np.vstack([npl,npr])
    
    
    if POSCAR == True:
            
        def POSCAR_file(scell, d, new_lat, n, labels, rel):
            s = scell
            rcell = np.hstack([scell,relce])
            vcell = scell[scell[:,3].argsort()] #organizes coordinates by atomic type (based on label)
            
            df = pd.DataFrame(vcell[:,:3])
            
            #CREATES SUPERCELL POSCAR FILE FOR EACH INTERFACE
            f4_name = os.path.join(p, "POSCAR")
            f4 = open(f4_name, "w+")
            f4.write(f"Interface {inf}\n")
            f4.write("1.000000 \n")
            f4.write(str(new_lat).replace(' [', '').replace('[', '').replace(']', '').replace('. ', '.0  '))
            f4.write("\n")

            f4.write(labels)
            f4.write("\n")
    
            f4.write(str(n).replace('.', ' ')[1:-1])
    
            f4.write("\n")
            f4.write("Direct \n")
            
            dfs = df.to_string(header=False, index=False)
            f4.write(str(dfs))
    
            f4.close()

            relcell = rcell[rcell[:,3].argsort()] #organizes coordinates by atomic type (based on label)
            relaxcell = np.hstack([relcell[:,:3],relcell[:,4:]]) #relax cell without atomic labels
            df = pd.DataFrame(relaxcell[:,:])
            
            if rel == True:
                #CREATES RELAX FOLDER
                path8 = f'{p}/relax'
                if not os.path.exists(path8):
                    os.makedirs(path8)
                else:
                    shutil.rmtree(path8) #Removes all the subdirectories! -> redo folder each time
                    os.makedirs(path8)
            
                f5_name = os.path.join(path8, "POSCAR")
                f5 = open(f5_name, "w+")
                f5.write(f"Relaxiation {xl}\{xr}")
                f5.write("1.000000 \n")
                f5.write(str(new_lat).replace(' [', '').replace('[', '').replace(']', '').replace('. ', '.0  ').replace("'", ""))
                f5.write("\n")
    
                f5.write(labels)
                f5.write("\n")
        
                f5.write(str(n).replace('.', ' ')[1:-1])
                  
                f5.write("\n")
                f5.write("Selective Dynamics \n")
                f5.write("Direct \n")
                
                dfs = df.to_string(header=False, index=False)
                f5.write(str(dfs))
        
                f5.close()
                
                return s,d
                
            else:
                path8 = f'{p}/relax'
                if not os.path.exists(path8):
                    pass
                else:
                    shutil.rmtree(path8) #Removes all the subdirectories!
                    
                return s,d
            

        supercell,dis = POSCAR_file(scell, dists , new_lat, n, labels, rel)
        
        return supercell,dis
    
        
    else:
        return scell,dists
        
# =============================================================================        
        




######  CREATE SUPERCELL VACUUM INTERFACES & ASSOCIATED POSCAR FILES  ######
# =============================================================================
def vac_interface(inf,p,xl,xr,yl,yr,zl,zr, atoms, atom_labels, lattice, points, POSCAR):

    #Total Number of Cells in Each Direction for Supercell
    xs = xl + xr
    ys = yl + yr
    zs = zl + zr

    #Number of Each Atom in Supercell 
    #Note: Array n will be in order atoms are ordered in input file from first to last coordinate
    #Example: [Sn, Sn, Co, Co, Co, Se, S] -> [2, 3, 1, 1] * (num of cells in x * in y * in z)
    n = np.array([])
    num_distinct_atoms = (Counter(atoms).keys()) #list of distinct atoms
    atom = atoms.tolist() #list of labeled atoms ordered from first to last original coordinate
    for i in num_distinct_atoms:
        num = xl*yl*zl*atom.count(i)
        n = np.hstack([n, np.array([num])])
    
    #Atom Element Labels
    labels = str(atom_labels).replace('[', '').replace(']','').replace('\'','').replace(',',' ')

    #Changing primitive vectors to dimensions of supercell
    new_lat = np.array([[xs*lattice[0,0], ys*lattice[0,1], zs*lattice[0,2]],
                        [xs*lattice[1,0], ys*lattice[1,1], zs*lattice[1,2]], 
                        [xs*lattice[2,0], ys*lattice[2,1], zs*lattice[2,2]]])
    
    #Original Coordinates with Labels
    cell = np.hstack([points,atoms.reshape(len(atoms),1)])

    #Creates Left Cells for Bulk
    npl = np.array([]).reshape(0,4)
    for i in range(xl):
        n_l = np.copy(cell)
        n_l[:,0]= (n_l[:,0]+i)/xs
        npl = np.vstack([npl,n_l])
        
        
    #Creates Vacuum Supercell
    vac_scell = npl
    
    if POSCAR == True:
            
        def POSCAR_file(vac_scell,new_lat,n,labels):
            vac_vcell = vac_scell[vac_scell[:,3].argsort()]
            
            df = pd.DataFrame(vac_vcell[:,:3])
            
            f6_name = os.path.join(p, "POSCAR")
            f6 = open(f6_name, "w+")
            f6.write(f"Bulk|Vacuum  {inf}x4\n")
            f6.write("1.000000 \n")
            f6.write(str(new_lat).replace(' [', '').replace('[', '').replace(']', '').replace('. ', '.0  '))
            f6.write("\n")

            f6.write(labels)
            f6.write("\n")
    
            f6.write(str(n).replace('.', ' ')[1:-1])
    
            f6.write("\n")
            f6.write("Direct \n")
            
            dfs = df.to_string(header=False, index=False)
            f6.write(str(dfs))
    
            f6.close()
            
            return(vac_scell)
        
        POSCAR_file(vac_scell,new_lat,n,labels)
        
    else: 
        return vac_scell
# =============================================================================





######  OBTAIN SYMMETRY OPERATION FROM DATABASE  ######
# =============================================================================

#Obtain symmetry name from database with symmetry operation matrix
def get_key(val, symmetry, rtol, atol):
    for key, value in symmetry.items():
        if np.allclose(val,value,rtol,atol):
            return key
    return "key doesn't exist"

# =============================================================================





######  CUSTOMIZE YOUR OUTPUT CRYSTAL SYSTEM FOLDER  ######
# =============================================================================



### CALCULATE DISTANCE BETWEEN CLOSEST ATOMS AT INTERFACE ###
def distances(plot_int, path, atoms, atom_labels, lattice, points, cgroup, size_range):
    
    x_distances = []
    tot_distances = []
    
    for i in range(len(plot_int)):
        
        l = size_range[0]
        
        xl_size = l
        xr_size = l
        yl_size = 1
        yr_size = 0 
        zl_size = 1
        zr_size = 0 
        
        s,D = supercell_interface(i, plot_int[i], path, xl_size, xr_size, yl_size, yr_size, 
                               zl_size, zr_size, atoms, atom_labels, lattice, points, cgroup, POSCAR=False, rel = False)
        
        x_distances.append(D[0])
        tot_distances.append(D[1])
        
    return x_distances, tot_distances




###  BULK FOLDER  ###
def bulk_folder(path, plot_int, comment, latconst, lattice, atom_labels, a, points, cgroup, symcells):

    #Creates bulk folder
    path2 = f'{path}/Bulk'
    if not os.path.exists(path2):
        os.makedirs(path2)
    else:
        shutil.rmtree(path2) 
        os.makedirs(path2)
    
    #Adds Primitive Cell POSCAR File to Bulk Folder
    f1_name = os.path.join(path2, "POSCAR")
    f1 = open(f1_name, "w+")
    f1.write(f"Primitive Unit Cell for {comment}")
    f1.write("\n")
    f1.write(str(latconst))
    f1.write("\n")
    f1.write(str(lattice).replace(' [', '').replace('[', '').replace(']', '').replace('. ', '.0  '))
    f1.write("\n")
    f1.write(str(atom_labels).replace('[', '').replace('\'','').replace(']',''))
    f1.write("\n")
    f1.write(str(a).replace('[', '').replace(']', ''))
    f1.write("\n")
    f1.write('Direct')
    f1.write("\n")
    f1.write(str(points).replace(' [', '').replace('[', '').replace(']', ''))
    f1.close()
    
    if symcells == True:
        for i in range(len(plot_int)):
            name = f'{i}'
            fill = name.zfill(5)
            path9 = f'{path}/Bulk/sym{fill}'
            if not os.path.exists(path9):
                os.makedirs(path9)
            else:
                shutil.rmtree(path9) #Removes all the subdirectories! -> redo folder each time
                os.makedirs(path9)
            
            cg = cgroup[i].copy()
            cg = cg[cg[:,3].argsort()]
         
            #Adds Transformed Unit Cell POSCAR File to Bulk Folder for Each Unique Interface
            f1_name = os.path.join(path9, "POSCAR")
            f1 = open(f1_name, "w+")
            f1.write(f"Sym {i} Unit Cell for {comment}")
            f1.write("\n")
            f1.write(str(latconst))
            f1.write("\n")
            f1.write(str(lattice).replace(' [', '').replace('[', '').replace(']', '').replace('. ', '.0  '))
            f1.write("\n")
            f1.write(str(atom_labels).replace('[', '').replace('\'','').replace(']',''))
            f1.write("\n")
            f1.write(str(a).replace('[', '').replace(']', ''))
            f1.write("\n")
            f1.write('Direct')
            f1.write("\n")
            f1.write(str(cg[:,:3]).replace(' [', '').replace('[', '').replace(']', ''))
            f1.close()
            
        return print("Bulk Folder with Sym Cells Complete")
        
    else:
        return print("Bulk Folder without Sym Cells Complete")




###  INTERFACE FOLDER  ###
def interface_folder(path, plot_int, atoms, atom_labels, lattice, points, cgroup, size_range, relax):

    path3 = f'{path}/Interface'
    if not os.path.exists(path3):
        os.makedirs(path3)
    else:
        shutil.rmtree(path3) 
        os.makedirs(path3)
     

    #Inside interface folder
    for i in range(len(plot_int)):
        name = f'{i}'
        fill = name.zfill(5)
        path4 = f'{path3}/sym{fill}'
        if not os.path.exists(path4):
            os.makedirs(path4)
        else:
            shutil.rmtree(path4)
            os.makedirs(path4)
            
        #Inside each symmetry folder
        for j in range(size_range[0],(size_range[1]+1)):
            path5 = f'{path4}/{j}x'
            if not os.path.exists(path5):
                os.makedirs(path5)
            else:
                shutil.rmtree(path5) 
                os.makedirs(path5)
            
            #Supercell Dimensions
            xl_size = j
            xr_size = j
            yl_size = 1
            yr_size = 0 
            zl_size = 1
            zr_size = 0 
    
            #Will add POSCAR file to each size supercell folder
            s,D = supercell_interface(i, plot_int[i], path5, xl_size, xr_size, yl_size, yr_size, 
                               zl_size, zr_size, atoms, atom_labels, lattice, points, cgroup, POSCAR = True, rel = relax)
            
    return print("Interface Folder Complete") 




###  VACUUM FOLDER  ###
def vacuum_folder(path, atoms, atom_labels, lattice, points, size_range):
    #Creates vacuum folder and different folders for each size supercell
    path6 = f'{path}/Vacuum'
    if not os.path.exists(path6):
        os.makedirs(path6)
    else:
        shutil.rmtree(path6) #Removes all the subdirectories! -> redo folder each time
        os.makedirs(path6)
        
    for j in range(size_range[0],(size_range[1]+1)):
        path7 = f'{path6}/{j}x'
        if not os.path.exists(path7):
            os.makedirs(path7)
        else:
            shutil.rmtree(path7) #Removes all the subdirectories!
            os.makedirs(path7)
            
        #Vacuum Supercell Input Parameters
        xl_size = j
        xr_size = 4
        yl_size = 1
        yr_size = 0 
        zl_size = 1
        zr_size = 0 
        
        vac_scell = vac_interface(j,path7, xl_size, xr_size, yl_size, yr_size, 
                                    zl_size, zr_size, atoms, atom_labels, lattice, points, POSCAR=True)
        
    return print("Vacuum Folder Complete")    




###  SYMMETRY OF CRYSTAL FILE  ###
def sym_of_cry_file(p, sym_of_cry, comment, symmetry, rtol, atol):
    
    #Numver of Crystal Symmetries
    num_sym_of_cry = len(sym_of_cry)
    
    f2_name = os.path.join(p, "Sym_of_Crystal.txt")
    f2 = open(f2_name, "w+")
    f2.write(f"Symmetries of {comment}")
    f2.write("\n")
    f2.write("\n")
    f2.write(f"Number of Crystal Symmetries: {num_sym_of_cry}")
    f2.write("\n")
    f2.write("\n")
    
    #Obtains name for each symmetry operation
    for i in sym_of_cry:
        f2.write(str(get_key(i, symmetry, rtol, atol)))
        f2.write("\n")
        f2.write(str(i))
        f2.write("\n")
        f2.write("\n")
    
    f2.close()
    
    return print("Symmetry of Crystal File Complete")




###  INFORMATION FILE  ###
def information_file(xd,td, symmetry, discarded, sym_of_cry, interface_sym, sgroup, path, comment, rtol, atol):
    num_operations = len(symmetry)      #Number of Symmetry Operations in Database
    num_dis_sym = len(discarded)        #Number of Discarded Symmetries
    num_sym_of_cry = len(sym_of_cry)    #Number of Crystal Symmetries (with translations)
    num_poss_int = len(interface_sym)   #Number of Symmetry Operations that Form an Interface
    num_unq_int = len(sgroup)           #Number of Unique Interfaces
    

    f3_name = os.path.join(path, "Information.txt")
    f3 = open(f3_name, "w+")
    f3.write(f"Information for {comment}")
    f3.write("\n")
    f3.write("\n")
    f3.write(f"Number of Applied Symmetry Operations: {num_operations}")
    f3.write("\n")
    f3.write(f"Number of Discarded Symmetries: {num_dis_sym}")
    f3.write("\n")
    f3.write(f"Number of Crystal Symmetries: {num_sym_of_cry}")
    f3.write("\n")
    f3.write(f"Number of Simple Interfaces: {num_poss_int}")
    f3.write("\n")
    f3.write(f"Number of Unique Simple Interfaces: {num_unq_int}")
    f3.write("\n")
    f3.write("\n")
    f3.write("Number of Symmetry Operations that Create Each Interface")
    f3.write("\n")
    f3.write("\n")
    
    #Obtains the number of operations that can form each unique interface
    n_syms = []
    for i in range(len(sgroup)):
        tup = sgroup[i].shape
        num_ops = tup[2]
        f3.write(f"{i}: {num_ops} symmetry operations, x distance = {xd[i]}, total distance = {td[i]}")
        f3.write("\n")
        n_syms.append(num_ops)
        
    f3.write("\n")
    f3.write("Symmetry Operations that Create Each Interface & Distances")
    f3.write("\n")
    
    #Lists the symmetry operations that form each unique interface
    for i in range(len(sgroup)):
        tup = sgroup[i].shape
        num_ops = tup[2]
        f3.write("\n")
        f3.write("\n")
        f3.write(f'Interface {i} : Num of Syms = {n_syms[i]}')     # Can you write the totd {line 412} here? it is the distance between interfaces
        f3.write("\n")
        f3.write(f'x distance = {xd[i]}')
        f3.write("\n")
        f3.write(f'total distance = {td[i]}')
        f3.write("\n")
        f3.write('-----------------------------------')
        f3.write("\n")
        
        #Obtains the name of each symmetry operation from database
        for j in range(num_ops):
            f3.write(str(get_key(sgroup[i][:,:,j], symmetry, rtol, atol)))
            f3.write("\n")
            f3.write(str(sgroup[i][:,:,j]))
            f3.write("\n")
            f3.write("\n")
    
    f3.close()
    
    return print("Information File is Complete")



