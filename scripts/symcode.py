#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 02:10:16 2022

@author: ginapantano
"""


#######  PACKAGE IMPORTS  #######
# =============================================================================
import numpy as np
import os
import shutil
import OPS as s #symmetry operation database
import functions as func

# =============================================================================





#######  CRYSTAL INPUT FOR PRIMITIVE UNIT CELL FROM FILE  #######
# =============================================================================
#NOTE: FILE CANNOT HAVE ANY SPACES BETWEEN LINES
#NOTE: INPUT FILE NEEDS TO BE IN SAME DIRECTORY OR YOU MUST SPECIFY PATH WITH POSCARFILE

poscarfile='aluminum.txt'   #path to initial POSCAR file


rf = open(poscarfile, 'r')
 
comment = rf.readline().strip()             #Name of Crystal
latconst = float(rf.readline().split()[0])  #Lattice Constant
lattice = np.zeros((3, 3))

#Primitive Vectors
lattice[0, :] = latconst * np.array([float(x) for x in rf.readline().split()[:3]])
lattice[1, :] = latconst * np.array([float(x) for x in rf.readline().split()[:3]])
lattice[2, :] = latconst * np.array([float(x) for x in rf.readline().split()[:3]])

#Atom Names in Unit Cell
atom_labels = np.array([x for x in rf.readline().split()])
#Number of Each Atom in Unit Cell
a = np.array([int(x) for x in rf.readline().split()])
#Number of Atoms in Unit Cell
na=np.sum(a)
#Number of Unique Atoms in Unit Cell
species=len(a)

#Atom Labels for Each Atomic Coordinate
atoms=[]
for x in range(species):
    atoms=np.append(atoms,np.ones(a[x])*int(1+x))
atoms=atoms.astype(int)

mode = rf.readline()

#Atomic Coordinates
points = np.zeros((na, 3))
for x in range(na):
    points[x, :] = np.array([float(x) for x in rf.readline().split()[:3]])


#Scale Coordinates to Single Unit Cell
for i in range(na):
    for j in range(3):
        if points[i,j] < 0:
            while points[i,j] < 0:
                points[i,j] = 1 + points[i,j]
        elif points[i,j] >= 1:
            while points[i,j] >= 1:
                points[i,j] = points[i,j] - 1


# =============================================================================





###### CREATE FOLDER FOR CRYSTAL SYSTEM ######
# =============================================================================

#Creates folder for subdirectories with POSCAR files
#Folder is labeled based on input files first comment line and placed in current directory
path = f'./{comment}'

if not os.path.exists(path):
    os.makedirs(path)
else:
    shutil.rmtree(path) #Removes all the subdirectories! -> redoes folder each time
    os.makedirs(path)

# =============================================================================





###### PREFORM SYMMETRY OPERATIONS ######
# =============================================================================

#Dictionary of crystallographic space group symmetry operations
symmetry = s.ops

discarded = []      #Symmetries where lattice vectors do not match
sym_of_cry = []     #Symmetries of the given crystal
interface_sym = []  #Symmetries that can form a possible interface

unique_coords = []  #Transformed sorted coordinates that do not match previous iterations of transformed coordinates.
plot_int = []       #Coordinates for each unique interface to plot in VESTA using first symmetry that creates it.
sym_num = []        #List of which interface is repeated for non-unique interface symmetry operations.
sgroup = {}         #Symmetry groups that form each interface.
cgroup = {}         #Transformed sorted coordinates for each interface symmetry grouped for each unique interface.
ind=0               #Starting index for dictionaries


number = [] #Labels which symmetry operation creates an interface from database.



for l, sym in enumerate(symmetry):
    
    #FOR TESTING SYMS WITHOUT TRANSLATION -> make own sym dictionary
    #t = np.array([0,0,0])
    #sym = np.hstack([sym, t.reshape(3,1)])
    
    #SYMMETRY DATABASE WITH TRANSLATION
    sym = symmetry[sym]
    
    #NOTE: sym, at this point, must have form 3x4
    #Preform Symmetry Operation on Lattice Vectors
    new_vecs = np.dot(sym[:,:3],lattice)
   
    #Preform Symmetry Operation on Atomic Coordinates
    tpoints = np.transpose(points[:,:3]) #transpose to get [[x1,x2,...], [y1,y2,...], [z1,z2,...]]
    point_sym = np.dot(sym[:,:3],tpoints) #point symmetry operation on coords

    #Creates Translational Symmetry Operation Matrix
    trans_sym = np.array([]).reshape(0,3)
    for i in range(na):
        #Adds new row to 2D array with translational symmetry operation until
        #the number of rows = number of atomic points
        trans_sym = np.vstack([trans_sym,sym[:,3]])
    
    ##Transformed Atomic Coordinates of Form [[x1,x2,...], [y1,y2,...], [z1,z2,...]] -> (point + translation symmetry)
    new_points = point_sym + np.transpose(trans_sym)
    

    #Scale Transformed Coordinates to Single Unit Cell
    for i in range(3):
        for j in range(na):
            if new_points[i,j] < 0:
                while new_points[i,j] < 0:
                    new_points[i,j] = 1 + new_points[i,j]
            elif new_points[i,j] >= 1:
                while new_points[i,j] >=1:
                    new_points[i,j] = new_points[i,j] - 1
    

    #New Atomic Coordinates
    #Transpose to get form [[x1,y1,z1], [x2,y2,z2]....]
    new_points = np.transpose(new_points)
    
    # =============================================================================
    
    
    ########  CHECK PRIMITIVE VECTORS MATCH  ########
    # =============================================================================
    #Note: If primitive vectors do not match, then it is not a symmetry of the crystal AND
    #not a symmetry capable of creating a perfect interface due to misaligned faces. The 
    #symmetry operation is then added to the list "discarded".

    if np.array_equal(np.sum(lattice, axis=1),abs(np.sum(new_vecs, axis=1))):
        #print("lattice vectors match", '\n')
        
    
    ######  CHECK IF SYMMETRY OF THE CRYSTAL OR FORMS INTERFACE  ######
    # =============================================================================
    #Note: If symmetry operation is a symmetry of the crystal, we add symmetry to list of
    #symmetries for the crystal (sym_of_cry). Else, symmetry operation forms an interface.
        
        #Sort Original Coordinates and Atomic Labels
        coords = points.copy()
        o_label = atoms.copy()
        o_label = o_label[coords[:,2].argsort()]
        
        coords = coords[coords[:,2].argsort()]                   #sort by "x" (col 1)
        o_label = o_label[coords[:,1].argsort(kind='mergesort')]
        coords = coords[coords[:,1].argsort(kind='mergesort')]   #sort by "y" (col 2)
        o_label = o_label[coords[:,0].argsort(kind='mergesort')]
        coords = coords[coords[:,0].argsort(kind='mergesort')]   #sort by "z" (col 3)
        
        #Sort Transformed Coordinates and Atomic Labels
        new_coords = new_points.copy()
        t_label = atoms.copy()
        t_label = t_label[new_coords[:,2].argsort()]
        
        new_coords = new_coords[new_coords[:,2].argsort()] 
        t_label = t_label[new_coords[:,1].argsort(kind='mergesort')]                #sort by "x" (col 1)
        new_coords = new_coords[new_coords[:,1].argsort(kind='mergesort')] #sort by "y" (col 2)
        t_label = t_label[new_coords[:,0].argsort(kind='mergesort')]
        new_coords = new_coords[new_coords[:,0].argsort(kind='mergesort')] #sort by "z" (col 3)
        
        #Check if coordinates are the same within a given tolerance (due to numerical error).
        #Note: if coordinates are the same after sorting, then the crystal remains
        #invariant under the symmetry operation and is thus a symmetry of the crystal.
        
        rtol = 1e-6  #relative tolerance
        atol = 1e-7  #absolute tolerance
        
        #Symmetry of the Crystal
        if np.allclose(coords,new_coords,rtol,atol) and np.array_equal(o_label,t_label):
            sym_of_cry.append(sym)


        #Symmetry that Can Form Interface
        else:
            interface_sym.append(sym)
            number.append(l)

            npts = new_points.copy()
            N = np.hstack([new_coords, t_label.reshape(len(t_label),1)])
            
            #If coordinates list is empty, append coordinates to list (only first iteration)
            if len(unique_coords) == 0:
                unique_coords.append(new_coords)
                plot_int.append(new_points)
                sgroup[ind] = sym
                cgroup[ind] = N
                result = -2
            
            #Check if sorted coordinates were created from previous symmetry
            #If yes, determine which interface coordinates match -> result = index of list
            #If no, append coordinates to list -> result = -1
            else:
                for n,c in enumerate(unique_coords):
                    if np.allclose(new_coords,c,rtol,atol):
                        result = n
                        break
                    
                    else:
                        result = -1
            
            if result == -2:
                pass
            
            
            elif result == -1:
                plot_int.append(new_points)
                unique_coords.append(new_coords)
                ind = ind+1
                sgroup[ind] = sym
                cgroup[ind] = N
              
            else:
                sym_num.append(result)
                sgroup[result] = np.dstack([sgroup[result],sym])
                #cgroup[result] = np.dstack([cgroup[result], N])
# =============================================================================   
        
        '''
        EDIT:
           - In else condition, add new condition to check mean or total energy. Is the
           interface worth exploring?
        '''        
        
    #Disgarded symmetries if lattice vectors do not match.       
    else:
        discarded.append(sym)
        # print("lattice vectors do not match", '\n')
        # print("original: ",'\n',lattice,'\n')
        # print("transformed:",'\n',new_vecs, '\n')
        # print('Cannot form interface with symmetry')
# =============================================================================   





######  CUSTOMIZE YOUR OUTPUT CRYSTAL SYSTEM FOLDER  ######
# =============================================================================

#Different Size Supercells You Want to Create
size_range = (4,12) #smallest size to largest size range along x-axis in terms of unit cells

#Calculates total distance and distance along x-axis between two closest atoms at interface
xd, td = func.distances(plot_int, path, atoms, atom_labels, 
                        lattice, points, cgroup, size_range) #Required to run

### OUTPUT OPTION ONE: Bulk Folder with or without Symmetry Unit Cells ###
#Creates bulk folder with primitive cell and with or without symmetry unit cells
func.bulk_folder(path, plot_int, comment, latconst, lattice, atom_labels, a, points, cgroup, symcells=True)

### OUTPUT OPTION TWO: Interface Folder with or without Relaxation Folder ###
func.interface_folder(path, plot_int, atoms, atom_labels, lattice, points, cgroup, size_range, relax = True)

### OUTPUT OPTION THREE: Vacuum Folder ###
func.vacuum_folder(path, atoms, atom_labels, lattice, points, size_range)

### OUTPUT OPTION FOUR: Symmetries of Crystal File ###
func.sym_of_cry_file(path, sym_of_cry, comment, symmetry, rtol, atol)

### OUTPUT OPTION FIVE: Information File ###
func.information_file(xd,td,symmetry, discarded, sym_of_cry, interface_sym, sgroup, path, comment, rtol, atol)

# =============================================================================
