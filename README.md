Implicit solvation model for VASP
==========================================

**Authors:** Kiran Mathew and Richard Hennig

**Collaborators:** Ravishankar Sundararaman, Kendra Letchworth-Weaver, Tomas Arias

http://vaspsol.mse.cornell.edu/

**For discussions and feedback:**

 https://groups.google.com/forum/#!forum/vaspsol 


**For the model details please see:** 

http://scitation.aip.org/content/aip/journal/jcp/140/8/10.1063/1.4865107

CONTENTS
=============
The two *.F files form the crux of the solvation model. 

The patch file, *interface_patch_535*, links *.F files to the original VASP code version 5.3.5.


PREREQUISITES
=============
[VASP](http://www.vasp.at/).5.3.5 source code.

Compiler and library requirements are the same as VASP.5.3.5. ([vasp wiki] (http://cms.mpi.univie.ac.at/wiki/index.php/Installing_VASP))

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


After installation, please see the test files available at http://vaspsol.mse.cornell.edu/

If you have used our code for your work, please cite our paper:
http://scitation.aip.org/content/aip/journal/jcp/140/8/10.1063/1.4865107