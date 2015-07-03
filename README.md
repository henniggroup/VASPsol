Implicit solvation model for VASP: linear model with ionic-screening
==========================================

**Authors:** Kiran Mathew and Richard Hennig

http://vaspsol.mse.ufl.edu/

**For discussions and feedback:**

 https://groups.google.com/forum/#!forum/vaspsol

**References:**
If you use the VASPsol package in your research, please cite the following paper:
Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways.
K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, and R. G. Hennig. J. Chem. Phys. 140, 084106 (2014).
[doi:10.1063/1.4865107] (http://dx.doi.org/10.1063/1.4865107)

PREREQUISITES
=============
[VASP](http://www.vasp.at/) version 5.3.5.

Compiler and library requirements are the same as that of VASP ([vasp wiki] (http://cms.mpi.univie.ac.at/wiki/index.ph\
p/Installing_VASP))

INSTALL
========

- Apply the interface patch to the original VASP source code.
```
    cd <VASP src directory>
    patch -p1 < <path to the interface patch file>
```
- After applying the patch, copy the *_k.F files to the VASP src directory:
- In the original VASP Makefile, put pot_lpcm_k.o pot_k.o object file names before pot.o in that order.
- ``` make clean ```
- ``` make ```
