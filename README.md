Implicit solvation model for VASP: linear model with ionic-screening
==========================================

**Authors:** Kiran Mathew and Richard Hennig

http://vaspsol.mse.ufl.edu/

**For discussions and feedback:**

 https://groups.google.com/forum/#!forum/vaspsol

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
