Calculates local (or position-resolved) density of states along z-axis.
This is usefull when invesigating interfaces between two semiconductors, see [here](https://arxiv.org/abs/1610.04119).
The program is a post-processing tool that reads output from the [Vienna Ab initio Simulation Package](http://vasp.at). 

The previous version of the program was based on the PROCAR file. This approach
led to a slightly deformed DOS since delocalized states were underestimated. The
present version also supports k-meshes.

# Installation
Just compile with `g++ -o dosloc dosloc.cc` and put the executable in your `bin` folder.

In oder to use the program one has to mofify the VASP source code. 


You will also need [Gnuplot](http://gnuplot.info) installed on your system.

# Usage

  -input: modified PARCHG files (summed in x & y direction) for all bands and k-points
  -the PARCHG files are calculated by a modified VASP

  -the dos in a slice at a certain y position (and at a certain eigenvalue) is calculated as:
   dos(y) = amount of charge in slice / volume of slice
  -the amount of charge corresponds to the number of states (2 electrons = 2 states = 1 band)
  -DOS units are 10^21 cm^-3 eV^-1

  -for some reason the VASP uses sigma=0.036 when sigma=0.05 is specified in the INCAR

  -we don't use symmetry, so all k-points have the same weight  

  -the program was tested: 
    -averaging the position-resolved DOS reproduces the total DOS
    -gives the same result as doslab3 (complete supercell is selected)
    -position-resolved DOS calculated with (Gamma+locdos2) and with (2x2x1 mesh + locdos3) looks similar

  -if the output file (or the corresponding .eps file) is too large you can change the for loops
   in the "write output" section to for (i = 0; i < nz; i=i+2)

  -use these settings in your INCAR:
   LPARD    = .TRUE.
   NBMOD    = 0
   LSEPB    = .TRUE.
   LSEPK    = .TRUE.


   

