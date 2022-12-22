# SymCode - Version 1.0
**SymCode** is a set of Python scripts for determining the symmetries and possible interfaces of various crystal systems. *WARNING*: This project is in its very early beta stages and has lots of errors and kinks to still work out. Any feedback is welcomed and strongly encouraged.  I included a list of known bugs below.

## Description 
### Purpose
A large portion of theoretical condensed matter physics research is focused on understanding the interesting emergent phenomena that arises at the interface between materials. For example, the two insulating perovskite oxides lanthanum aluminate (LaAlO3) and strontium titanate (SrTiO3) when interfaced together become highly conductive [1]. Or, one can think of the vast electronic devices that have come from interfacing different semiconducting and metallic materials together [2]. The exploration of new materials and interfacial effects is essential to overcome modern technological constraints with energy generation, conversion, storage, and overall device performance [3]. Theoretical studies require numerical computation and mathematical tools in order to understand these complex interfacial phenomena, and one important aspect is to understand the change in symmetry at the interface. The break in symmetry is often responsible for the physical and chemical changes that lead to such interesting properities.  
  
SymCode was originally created with the intended purpose of studying symmetries of nanoscale crystallographic domain boundaires, interfaces between different oriented variations of the same crystalline species, which are often developed during crystal growth processes. Specifically, our group is interested in studying a certain class of materials whose electronic properties could be enhanced rather than reduced by the presence of crystal imperfections (i.e. crystallographic domain boundaries). 

SymCode is currently limited to crystal systems with primitive vectors in a conventional cartesian basis of the form a = [1,0,0], b = [0,1,0], and c = [0,0,1] where a, b, c are the crystal axes. Additionally, SymCode only considers perfectly aligned crystallographic domain boundaries where the magnitude of the primitive vectors remains unchanged. The second version of the code is underdevelopment to be able to transform matrix operations to the users basis of choice and account for different types of interfaces such as heteorstructures, terminations, and misaligned planes (see more details below under "Future Edits"). 

### Features
* Access to an extensive database (OPS) of the point and translation symmetry operations for the seven crystal systems. Note: The operation matrices are written in a cartesian basis ($[1,0,0], [0,1,0], [0,0,1]$), and the database does not contain screw axis and glide plane symmetries with a sequential pure translation.
* Determines the symmetry operations of a given crystal and prints the matrix representations to a text file.
* Determines the number of symmetry operations that can form a perfect domain boundary interface. 
* Determines the number of unique interfaces that can form and groups the symmetry operations that form each unique interface into a easily accessible dictionary.
* Generates supercell POSCAR files spanning along the x-axis by n unit cells for each unique interface.
* Generates POSCAR files to perform structural relaxtion with VASP using the calculated total energy.
* Generates supercell POSCAR files spanning along the x-axis by n unit cells for bulk|vacuum interface. 
* Generates information text file.
    * Contains the subgroup of operations that form each unique interface.
    * The total distance and distance along the x-axis between the two closest atoms at the interface.
    * The number of symmetry operations that form each unique interface.
    * and more (see below).




### Operation Database - OPS
The matrix operations in *OPS.py* that is utilized by the Python script *symcode.py* were taken from the International Tables for Crystallography [4]. The script *OPS.py* has a dictionary called syms that contains all of the possible point symmetry operations for the seven crystal systems. The script then loops through the possible screw and glide translation vectors to create the possible screw and glide plane symmetries as well as the pure translation symmetry operations. However, the database does not account for screw axis or glide plane symmetries sequentially followed by a pure translation. This feature will be added to the second version of SymCode. The dictionary keys represent the symbol of the symmetry operation, the orientation of the symmetry element, and the transformed coordinates in terms of x, y, and z. For reference, the symmetry operations for each space group can be found in detail in the International Tables for Crystallography [4].

### Procedure of Symmetry Search
The symmetry search procedure was partly inspired by the open source code Spglib which is a C library for finding and handling crystal symmetries. Specifically, identifying the space group and magnetic space group type, and the associated symmetry operations of the crystal [5]. The procedure for finding the symmetries and possible interfaces of given crystal for SymCode is as follows:

1. Execute the Python script *symcode.py* after specifiying the input POSCAR file path assigning to the variable *poscarfile* and customizing the output folder by commenting/uncommenting desired functions to be called in the imported Python script *functions.py*.
2. Using the symmetry database imported from the Python script *OPS.py*, the matrix operations are applied to the atomic coordinates specified by the input POSCAR file. 
3. The code first checks if the new lattice vectors, after performing the symmetry operation, matches the original primitive vectors such that a perfect aligned boundary is formed. If they do not match, the symmetry operation is added to the list *discarded*.
4. If the new vectors match the original primitive vectors, the code then checks if the symmetry operation is a symmetry of the crystal. This is done by sorting the original and transformed atomic coordinates with their element labels by the smallest x value, then y, and then z while maintaining each previous sort. If the order of the atomic labels and sorted coordinates match, then its a symmetry of the crystal (crystal remained invariant under the operation) and is added to the text file *sym of cry*. If not, then it will form an interface.
5. If it forms an interface, the code then creates what I call the first "subgroup" for interface 0. This process gets repeated for all symmetry operations.
6. Each time a symmetry yields an interface it will check whether the interface was previously created or not. If yes, the symmetry operation is added to subgroup associated with that interface, if no, the symmetry operation is added to a new subgroup.

By the end, we have the number of unique interfaces and the symmetry operations that create each interface. Whichever interface has the most associated symmetries seems to be the most stable and experimentally viable interface based on our groups current DFT calculations. We still have a lot yet to discover and understand. Note: atoms are allowed to deviate from their ideal positions within some small tolerance to account for numerical error.


## INPUT
The only input file required by the user to execute the Python script *symcode.py* is the initial POSCAR file for the given crystal. The POSCAR file specifies the periodic simulation of the cell and information regarding the geometry of the system. The file must be structured in the following manner []:

~~~
Comment line  - (name of the output folder)
Universal scaling factor
Primitive lattice vectors
    a = 1.0  0.0  0.0 
    b = 0.0  1.0  0.0
    c = 0.0  0.0  1.0
Specify atom species (element names)
Specify number of each atom species
Type of coordinates: Cartesian or Direct
Position of atoms
~~~

Notice there is intentionally no spaces between the lines of the file. The comment line is often the crystals chemical formula and for our purpose will be the name of the output folder after executing the Python script *symcode.py*. The universal scaling factor is used to scale the given lattice vectors to the actual size of the primitive cell. Currently, the code only works for primitive vectors in a cartesian basis of the form shown above. In the second version of SymCode, the only restriction will be that the primitive vectors have to be linearly independent. The next line specifies the different type of atoms present in the cell. The order must match the order the positions of the atoms are given in. The next line specifies the number of each atom species present in the cell. The next line specifies the type of coordinates the positions are given in. In direct coordinates, the positions of the atoms are given as fractional coordinates of the primitive vectors. In cartesian coordinates, the positions of the atoms are given as global coordinates. The scaling factor would then scale the atomic coordinates. ***Please only write "Direct" for this line***. The rest of the file are the coordinates of each atom in the cell. 

Here is an example file for face center cubic (FCC) aluminum:
~~~
Al
1.0
4.0  0.0  0.0
0.0  4.0  0.0
0.0  0.0  4.0
Al
4
Direct
0.0  0.0  0.0
0.5  0.0  0.5
0.5  0.5  0.0
0.0  0.5  0.5
~~~
Please see the folder "tests" for other example POSCAR files.

## OUTPUT
The ouput folder will be added to your working directory and named the same as the first comment line in the initial POSCAR file. Please see the aluminum example output folder under the folder *example*. This has all the possible features SymCode currently has to offer.

### OPTION 1: Bulk Folder
The bulk folder will contain the initial POSCAR file representing the primitive unit cell of the crystal. Additionally, the user can also specify if they would like the POSCAR files for each transformed primitive unit cells for each of the unique symmetries that creates an interface with the original crystal. Please see the user instructions below how to turn this optional feature on and off. 

### OPTION 2: Interface Folder
The interface folder will contain a folder for each unique symmetry operation that creates a unique interface. Inside each symmetry folder are subfolders containing the various sized supercell POSCAR files. The size range of the supercells can be specified by the user in the Python script *symcode.py*. Please see the user instructions below. An additional option the user can turn on and off is to create the folder *relax* inside each supercell folder.  

#### Relax Folder
The relax folder contains the POSCAR file needed to perform structural relaxation methods to optimize the atomic coordinates and/or the cell with respect to the total energy calculated using VASP. The file indicates which atoms are allowed to move or relax by marking the coordinate position with logical flags. "T" marks which coordinates are allowed to change and "F" marks which coordinates are fixed. The user would then need to add the comment line "Selective Dynamics" above the comment line "Direct" in the POSCAR file to run with VASP. The user can specify which atoms get marked by indicating the number of unit cells away from the interface that atoms are allowed to move. This is currently set as 2 unit cells. The only way to change this is by editing the Python script *functions.py* by reassigning the variable *nrel* under the function *supercell_interface* to the new interger value. The next version of SymCode will add this feature to the function *interface_folder* input parameters. 

### OPTION 3: Vacuum Folder
The vacuum folder contains the POSCAR files for various sized supercells with bulk|vacuum interface. The right handside represents the vacuum which spans *n* number of empty unit cells (please see user instructions how to specify the size of supercells). The left handside from the interface are the bulk unit cells. The number of cells varies based on the size range specified by the user. 

### OPTION 4: Sym of Cry File
The symmetry of the crystal text file contains all of the symmetry operations that leave the crystal invariant. The file specifies the number of crystal symmetries and then sequentially lists the name of the operation, the orientation of the symmetry element, the transformed coordinates with translation vector, and its matrix representation. Here is an example of the first few lines of the symmetry file for aluminum. 

~~~
Symmetries of Al

Number of Crystal Symmetries: 192

Identity : x, y, z [[0 0 0]]
[[1 0 0 0]
 [0 1 0 0]
 [0 0 1 0]]

...
...
...
~~~

### OPTION 5: Information File
The information text file contains the number of applied symmetry operations, the number of discarded symmetries, the number of crystal symmetries, the number of operations that form an interface, the number of unique interfaces, and the number of symmetry operations that form each unique interface. The succeeding text are the matrix representations of the symmetry operations that form an interface, which are grouped together by which unique interface they create. Additionally, the total distance between the closest atoms at each unique interface and the distance between their x components is listed.  


## User Instructions

### Importing Python scripts and initial input file
First, make sure the Python scripts are all in the same working directory. If the Python scripts *functions.py* and *OPS.py* are not in the same directory as the Python script *symcode.py*, then the user must specify the path of the imported files in *symcode.py*. For example,
~~~
import sys
sys.path.append('/path/to/folder/with/file')
import OPS
~~~

If the chosen file path is incorrect, Python will raise the following error:
~~~
ModuleNotFoundError: No module named 'file name'
~~~
The path variable contains the directory the Python interpreter searches through to find the file, so if the file name does not exist within that directory, the above error is raised. Next, the user must specify the path where the input POSCAR file is located and assign the path to the variable *poscarfile* at the top of the Python script *symcode.py*:
~~~
poscarfile = 'path/to/input/file'
~~~

### Customizing the output folder
The user is able to control which features they would like to use by commenting/uncommenting the functions at the end of the *symcode.py* file. For example, if you only want the symmetries of the crystal, you would only need the function that outputs the *sym of cry* text file (option 4). Thus, you would comment out the functions for option 1, 2, 3, and 5.  

#### Specifying size range of supercells (optional)
If the user wants option 2 or 3, they are able to specify the size range of the supercell spanning along the x-direction. The size range is currently set as 4 unit cells to 12 unit cells. The range is represented as a tuple assigned to the variable *size_range* at the end of the *symcode.py* file. The first element of the tuple will specify the size of the vacuum (which stays constant for all bulk|vacuum cases), and the minimum size of the bulk/symmetry cells to right and left side of the interface respectively. The second element of the tuple specifies the maximum size the bulk/symmetry cells can be.

#### Plotting in VESTA 
VESTA is a free and user-friendly software for visualizing crystal structures [6]. The output POSCAR files can easily be plotted with VESTA. The user just needs to add the .vasp extension to the file name in order to open with the VESTA application.

### Repository folders & files
* **scripts** - this folder contains all three Python scripts: *functions.py, OPS.py, and symcode.py*.
    * *symcode.py* - the main Python file that is edited and executed by the user.
    * *functions.py* - this Python file contains all of the functions called by *symcode.py*. 
    * *OPS.py* - Python file containing the symmetry operation database called on by *symcode.py*. 
* **tests** - this folder contains several example input POSCAR files the user can test SymCode with. 
* **example** - this folder contains an example ouput folder for aluminum. This has all the possible features SymCode currently has to offer.
* **environment.yml** - contains a description of my conda environment I have been using for this project. This includes the versions of all the software and packages used (or will be used) for the project.

## Known Bugs
* Only works for systems whose primitive vectors are in a conventional cartesian basis. Thus, SymCode currently only works for cubic, tetragonal, and orthorhombic systems.
* SymCode only checks for perfectly aligned crystallographic domain boundary interfaces. 
* The symmetry database is incomplete. Does not account for screw axis or glide plane symmetries with sequential pure translations. An example case of this is with silicon (input POSCAR file under *tests*). The symmetry of the crystal text file is not a complete list of all the possible symmetries of silicon.
* SymCode input file can only have the VASP POSCAR file format.

## Future Edits
* Expand SymCode to accept any basis vectors.
* Produce heterostructures with two given input files and other types of interfaces.
* Determine the space group type for the crystal based on the given crystal symmetry operations.
* Expand SymCode to check magnetic symmetries and determine crystal magnetic space group type.
* Improve documentation.
* Create SymCode as a pip installable package to avoid users having to reference paths to imported files or keeping files in the same directory.
* Expand SymCode to accept various input file types such as .xsf structure file.

## References
[1] Ohtomo, A., & Hwang, H. Y. (2004). A High-Mobility Electron Gas at the LAALO3/SRTIO3 Heterointerface. Nature, 427(6973), 423–426. https://doi.org/10.1038/nature02308   

[2] Riordan, M., Hoddeson, L., & Herring, C. (1999). The Invention of the Transistor. Reviews of Modern Physics, 71(2). https://doi.org/10.1103/revmodphys.71.s336 

[3] Tsymbal, E. Y., & Dowben, P. A. (2013). Grand Challenges in Condensed Matter Physics: From Knowledge to Innovation. Frontiers in Physics, 1. https://doi.org/10.3389/fphy.2013.00032   

[4] Arnold, H., Bertaut, E. F., Billiet, Y., Buerger, M. J., Burzlaff, H., Fischer, W., Zimmermann, H., Wondratschek, H., de Wolff, P. M., Müller, U., Looijenga-Vos, A., Koch, E., Klapper, H., Hahn, T., & Gruber, B. (2006). International Tables For Crystallography. International Tables for Crystallography, A(5), 815–816. https://doi.org/10.1107/97809553602060000100   

[5] Togo, A., & Tanaka, I. (2018, August 5). Spglib: A Software Library for Crystal Symmetry Search. arXiv.org. Retrieved December 21, 2022, from https://arxiv.org/abs/1808.01590   

[6] Momma, K., & Izumi, F. (2011). VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data. Journal of Applied Crystallography, 44(6), 1272–1276. https://doi.org/10.1107/s0021889811038970 


## Contact Information
**Authors:** (1) Gina Pantano, (2) Dr. Jacob Gayles  
**Emails:** (1) gmpantano@usf.edu, (2) gayles@usf.edu  
**LinkedIn:** www.linkedin.com/in/gina-pantano  
**Group Website:** http://faculty.cas.usf.edu/jgayles/  


```python

```
