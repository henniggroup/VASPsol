Implicit solvation model for VASP: nonlinear model
==========================================

**Authors:** Kiran Mathew and Richard Hennig

http://vaspsol.mse.ufl.edu/

**For discussions and feedback:**

 https://groups.google.com/forum/#!forum/vaspsol

**References:**

If you use the VASPsol package in your research, please cite the following paper:

Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways.
 K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, and R. G. Hennig. J. Chem. Phys. 140, 084106 (2014), [doi:10.1063/1.4865107] (http://dx.doi.org/10.1063/1.4865107).

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
- After applying the patch, copy the *.F files to the VASP src directory:
- In the original VASP Makefile, put pot_lpcm_k.o pot_nln_k.o pot_k.o solvation.o object file names before pot.o in that order.
- ``` make clean ```
- ``` make ```

```
For VASP6 versions >= 6.1.0:
----------------------------
For VASP versions > 6.1.0 there is currently a patch in /VASPsol/patches/VASPsol6.patch. This patch applies the FERMI_SHIFT patch to the pot.F as well as providing additional subroutines. This patch was developed by user shk11 in Issue #43. 

- copy solvation.F from path_to_VASPsol/src/solvation.F to path_to_VASP6_install/src/
- Navigate to path_to_VASP6_install/src/ and modify the .objects file by ensuring that “solvation.o” appears before “pot.o”
- Set the CPP option "-Dsol_compat" in the VASP makefile.include file. 
- Compile the code as described in the instructions on the VASP wiki.

```   
    cd vasp.6.1.0
    patch -p0 < <path to the pbz_patch_610 patch file>
```
