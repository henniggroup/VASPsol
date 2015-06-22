INCAR PARAMETERS: Example
=========================
- LSOL = .TRUE., Solvation on/off
```
	 NOTE :
	 	 for linear model just set LSOL and leave LNLSOL
		 again default is water. For other solvents set EB_K and TAU params
		
		 to include ionic screening just set the debye length, lambda_d_k
		  default lambda_d_k corresponds to 1M ionic conc at 298K --> 3.04 Angstroms

		  CAUTION: linear ionic screening model is applicable only in small potential
		  regime --> do not use the nonlinear cut off charge density which is 
		  almost 3 times that of the linear model cutoff ==> larger potentials 
		  will be present in system
```
- LAMBDA_D_K = 3.04 # debye length in Angstroms 1/lambda_d^2 = kappa^2  = sum(n_i * z_i**2) / (epsilon_r * epsilon_0) /kT, 

- LNLSOL = .TRUE., Nonlinear on/off switch


 Solvent info example: water @ 298K (the default)
----------------------------------------------------
- EB_K = 78.4, dimensionless static relative permittivity of the solvent
 
- NSOL_K = 55.35, solvent conc in mol/L

- EINFTY_K = 1.78, dimensionless high frequency relative permittivity of the solvent

- P0_K = 2.3892, permenant dipole moment of the solvent molecule in Debye

- TAU =  9.23E-4_q, surface tension parmameter in eV/Ang^2


 Electrolyte info example: KClO4
----------------------------------------------------
- NION_COUNT_K = 2, number of ions:  K+ and ClO4-

- NION_K = 1 1, conc of ionic species in mol/L

- ZION_K = 1 -1, charge on each ions : K+, ClO4-

- RION_K = 1.52 2.26, ionic radii in Angstroms
```
		if lnsol is set to true, always set the ionic radii else the akappa2 contribution
		to the energy will blow up because of the large value of phi(reminder: nc is pretty large for the nonlinear sol)
		default is some meaningless impossibly large number , around 500 Angstrom
```

- LRHOB = .TRUE.
```
	write out the bound charge density to the file RHOB
	only for nonlinear
```
