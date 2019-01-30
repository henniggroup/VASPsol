- For installation please refer to README.md

- Do a vacuum calculation for your system first and save the wavefunction file WAVECAR by 
  specifying LWAVE = .TRUE. in the INCAR file.

- Start the solvation calculation from the vacuum WAVECAR, specify ISTART = 1 in the INCAR file.

- The solvation parameters are read from the INCAR file.

- In the simplest case the only parameter that need to be set is the solvation flag LSOL = .TRUE.
```
    This corresponds to the default case where the solvent is water with a 
    relative permittivity of 80. Note, The cavitation energy in the solvation 
    model requires a higher resolution than neeed for vacuum calculations. 
    So please make sure that ENCUT is large enough and also set PREC=Accurate.
```

- The relative permittivity of the solvent can be changed from the default value of 78.4 for water 
  to any value using the INCAR parameter EB_k, e.g. EB_k = 20.

- You can also change the surface tension parameter tau in the INCAR file. For example, to neglect 
  the the cavitation energy contribution, set TAU = 0 in the INCAR file. Note that you may have to 
  increase the cutoff energy to converge the caviation energy because the grid must be fine enough 
  to resolve the cavity surface. 

- Set LRHOB = .TRUE. if you would like to write out the bound charge density in the CHGCAR format. The file is named RHOB.
```
    Note: Inorder to keep the interface between the vaspsol and the rest 
    of the vasp ecosystem simple, the RHOB fileio is done within the module. This 
    means the RHOB file will be written in each scf iteration. This could slow down 
    the calculations for bigger systems. So it is reccomended that you do the solvation 
    calculations without this parameter and do another static calculation(starting 
    from the converged WAVECAR)with RHOB=.TRUE. if you wish to visualize the bound 
    charge. Sorry for the inconvenience. I will try to isolate the file io without 
    creating too much of a mess in the rest of the VASP code.
```

- For examples, please see the Examples folder

- Set LAMBDA_D_K (the Debye length in Angstroms) parameter in the INCAR file to use the linearized Poisson-Boltzmann model(electrolyte model)
```
    The VASP choice of the electrostatic reference potential sets the average potential
    in the simulation cell to zero, not the potential in the electrolyte region. To
    correct the reference electrostatic reference potential, the constant, FERMI_SHIFT,
    which is printed out to stdout, needs to be added to the Fermi level. The sign of
    the shift becomes apparent once you plot the x-y averaged local potential in
    the z-direction. Without the shift, the electrostatic potential would not go to zero
    in the bulk of the electrolyte as expected from the Poisson-Boltzmann solution.
    
    Furthermore, if you shift the electrostatic potential to align the potential in the
    electrolyte region to zero, or any other value, you also need to add a correction to
    the energy of your system. This correction is given by Q V, where  Q is the net
    charge of the simulation cell and V the shift in reference potential, e.g.,
    V = FERMI_SHIFT. 
```
