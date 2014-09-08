VASPsol: Implicit solvation model for VASP
==========================================

**Authors:** Kiran Mathew and Richard Hennig

http://vaspsol.mse.cornell.edu/

**For discussions and feedback:**

 https://groups.google.com/forum/#!forum/vaspsol 


**For the model details please see:** 

http://scitation.aip.org/content/aip/journal/jcp/140/8/10.1063/1.4865107

CONTENTS
=============
The 2 files, *_k.F, form the crux of the solvation model. 

The file interface_patch_535 is the patch file for linking *_k.F files to the original VASP code version 5.3.5.


PREREQUISITES
=============
VASP.5.3.5 source code.

Compiler and library requirements are the same as VASP.5.3.5.

INSTALL
========

- Apply the interface patch to the original VASP source code.
```   
    cd <VASP src directory>
    patch -p1 < <path to the interface patch file>
```
- After applying the patch, copy the other 2 files to the VASP src directory:
```
   	 git clone https://github.com/matk86/VASPsol <path to VASP src directory>
```
- In the original VASP Makefile, put pot_lpcm_cav_k.o and pot_k.o object file names before pot.o in that order.
- ``` make clean ```
- ``` make ```


After installation, please see the test files available at http://vaspsol.mse.cornell.edu/