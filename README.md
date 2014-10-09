Implicit solvation model for VASP: Monopole, Dipole/quadrapole corrections due to the solvent presence
==========================================

**Authors:** Kiran Mathew and Richard Hennig

http://theory.mse.cornell.edu/

PREREQUISITES
=============
[VASP](http://www.vasp.at/) version 5.3.5.

Compiler and library requirements are the same as that of VASP ([vasp wiki] (http://cms.mpi.univie.ac.at/wiki/index.php/Installing_VASP))

INSTALL
========

- Apply the interface patch to the original VASP source code.
```   
    cd <VASP src directory>
    patch -p1 < <path to the interface patch file>
```
- After applying the patch, copy the other 2 files to the VASP src directory:
- In the original VASP Makefile, put pot_lpcm_cav_k.o and pot_k.o object file names before pot.o in that order.
- ``` make clean ```
- ``` make ```
