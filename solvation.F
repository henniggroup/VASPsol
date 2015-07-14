!**********************************************************************
! Copyright 2013-2015
! Kiran Mathew and Richard Hennig
!
! This file is part of VASPsol.
!
! VASPsol is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! VASPsol is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with VASPsol.  If not, see <http://www.gnu.org/licenses/>.
!**********************************************************************

#include "symbol.inc"

!******************** MODULE SOLVATION *********************************
!
!
! interfaces the solvation engine with the rest of vasp
!
!
!***********************************************************************
MODULE solvation

  USE prec
  USE base
  USE pot_k

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SOL_READER, SOL_WRITER, XML_WRITE_SOL, SOL_VCORRECTION
  ! vaclevel
#ifdef debugsol_pbz
  PUBLIC :: GET_FERMISHIFT
#endif
  ! vaclevel

  LOGICAL, SAVE :: LSOL=.FALSE.
  LOGICAL, SAVE :: LNLSOL=.FALSE.
  LOGICAL, SAVE :: LRHOB=.FALSE.
  LOGICAL, SAVE :: LRHOION=.FALSE.

  !energy and force
  REAL(q), PUBLIC, SAVE :: Ediel_SOL = 0._q
  REAL(q), PUBLIC, ALLOCATABLE, SAVE :: EIFOR_SOL(:,:)
  
  REAL(q), ALLOCATABLE, PRIVATE, SAVE :: NION_K(:), ZION_K(:), RION_K(:)
  ! shape function parameters for linear solvation
  ! linearpcm params, 
  !if LNLSOL = .TRUE. then nc_k and tau will be reset
  REAL(q), PRIVATE, SAVE  :: sigma_k = 0.6_q
  REAL(q), PRIVATE, SAVE  :: nc_k = 0.0025_q !! 0.00473_q
  REAL(q), PRIVATE, SAVE  :: tau = 5.25E-4_q

  !default values for the paramters that describe the solvent environment
  ! water @ 25C
  REAL(q), PRIVATE, SAVE  :: eb_k = 78.4_q !80.0_q
  REAL(q), PRIVATE, SAVE  :: p0_k = 0.94/0.393430_q ! in Debye, 1 D = 0.393430 ea0
  REAL(q), PRIVATE, SAVE  :: einfty_k = 1.78_q !1.33 * 1.33 ! square of water refractive index 
  REAL(q), PRIVATE, SAVE  :: Nsol_k = 55.35_q !@ 25C, 997.075/18.015 = 55.35!!!55.5 !mol/L, for water
  REAL(q), PRIVATE, SAVE  :: T_k = 298.15_q !!!25C, temperature in kelvin
  !debye length, set to infinity
  REAL(q), PRIVATE, SAVE  :: lambda_d_k = 1E100_q
  !kappa2b_k =  (1/debye screeing lenth in ang)**2
  ! to convert to au, mulitply by 0.5292**2
  !debye screening length = 3 Ang ==> (1./3^2) * (0.5292)^2 = 0.03
  ! To override set lambda_d_k in the INCAR file
  REAL(q), PRIVATE, SAVE  :: kappa2b_k = 0._q 

  !default ionic species parameters
  ! default electrolyte: 1M NaF
  !number of ionic species, used to 
  !allocate space for the other ionic params
   INTEGER, PRIVATE, SAVE :: NION_COUNT_K = 2 

  REAL(q), PRIVATE, SAVE  :: EDIFFSOL ! set in module pot_k

!local variables
  TYPE (in_struct), PRIVATE, SAVE :: IO_1
  REAL(q), PRIVATE, SAVE :: EDIFF_1

CONTAINS

!******************** SUBROUTINE SOL_READER ****************************
!
!
! Reads in the solvation model parameters
!
!      
!***********************************************************************
  SUBROUTINE SOL_READER(NIONS,EDIFF,IO)

    USE base
    USE vaspxml
    
    IMPLICIT NONE

    TYPE (in_struct), INTENT(in) :: IO
    REAL(q), INTENT(in) :: EDIFF
    INTEGER, INTENT(in) :: NIONS
    INTEGER :: IDUM, N, IERR
    REAL(q) :: RDUM
    COMPLEX(q) CDUM
    LOGICAL :: LOPEN,LDUM,LIONODE_FLAG,LFLUEX
    CHARACTER (1) CHARAC
    
    ! this has to be done ALWAYS
    IF (ALLOCATED(EIFOR_SOL)) DEALLOCATE(EIFOR_SOL)
    ALLOCATE(EIFOR_SOL(3,NIONS))

    !ionic concentrations in mol/L
    !default is no electrolyte
    IF (ALLOCATED(NION_K)) DEALLOCATE(NION_K)
    ALLOCATE(NION_K(NION_COUNT_K))
    NION_K = 0._q
    !charges on ionic species (dimensionless)
    IF (ALLOCATED(ZION_K)) DEALLOCATE(ZION_K)
    ALLOCATE(ZION_K(NION_COUNT_K))
    ZION_K(1) = 0._q
    ZION_K(2) = 0._q
    !ionic radii (Angstrom)
    !default set to some impossible value
    IF (ALLOCATED(RION_K)) DEALLOCATE(RION_K)
    ALLOCATE(RION_K(NION_COUNT_K))
    RION_K(1) = 1000._q * AUTOA
    RION_K(2) = 1000._q * AUTOA

    EIFOR_SOL = 0._q
    Ediel_SOL = 0._q
    IO_1 = IO
    EDIFF_1 = EDIFF
    EDIFFSOL = EDIFF/100._q

    LOPEN = .FALSE.
   
    OPEN(UNIT=IO%IU5,FILE='INCAR',STATUS='OLD')      

  !---------------------------------------------------------------
  !LSOL
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'LSOL', 'L', 1, IDUM, &
       & RDUM, CDUM, LSOL, 'switch for solvation')

  IF (LSOL) THEN

  !---------------------------------------------------------------
  !sigma_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'SIGMA_K', 'F', 1, IDUM, &
       & SIGMA_K, CDUM, LDUM, 'the width of dielectric cavity (dimensionless)')

  !---------------------------------------------------------------
  !nc_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'NC_K', 'F', 1, IDUM, &
       & NC_K, CDUM, LDUM, 'cavity turn-on charge density (Angstrom^-3)')

  !---------------------------------------------------------------
  !eb_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'EB_K', 'F', 1, IDUM, &
       & EB_K, CDUM, LDUM, 'relative permittivity of the bulk solvent (dimensionless)')

  !---------------------------------------------------------------
  !tau
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'TAU', 'F', 1, IDUM, &
       & TAU, CDUM, LDUM, 'cavity surface tension (eV/Angstrom^2)')

  !---------------------------------------------------------------
  !EDIFFSOL
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'EDIFFSOL', 'F', 1, IDUM, &
       & EDIFFSOL, CDUM, LDUM, 'solvation convergence tolerance')

  !---------------------------------------------------------------
  !LRHOB
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'LRHOB', 'L', 1, IDUM, &
       & RDUM, CDUM, LRHOB, 'write the bound charge file')

  !---------------------------------------------------------------
  !LAMBDA_D_K
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'LAMBDA_D_K', 'F', 1, IDUM, &
       & LAMBDA_D_K, CDUM, LDUM, 'Debye length (Angstroms)')

  !kappa^2 from debye length
  IF (LAMBDA_D_K > 1E30_q) THEN
     KAPPA2B_K = 0._q
  ELSE
     KAPPA2B_K = 1._q / (LAMBDA_D_K ** 2) * AUTOA**2  ! in 1/a.u^2
  ENDIF

  !---------------------------------------------------------------
  !LNLSOL
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(IO, LOPEN, 'LNLSOL', 'L', 1, IDUM, &
       & RDUM, CDUM, LNLSOL, 'switch for nonlinear solvation')

  ! read the nonlinear solvation parameters
  IF (LNLSOL) THEN
     !reset the default solvation params, nc_k and tau
     nc_k    = 0.00675_q 
     tau     = 9.23E-4_q ! for water
     !---------------------------------------------------------------
     !P0_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(IO, LOPEN, 'P0_K', 'F', 1, IDUM, &
          & P0_K, CDUM, LDUM, 'dipole moment of the molecule (Debye)')

     !---------------------------------------------------------------
     !EINFTY_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(IO, LOPEN, 'EINFTY_K', 'F', 1, IDUM, &
          & EINFTY_K, CDUM, LDUM, '\epsilon_{\infty} (dimensionless)')

     !---------------------------------------------------------------
     !NSOL_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(IO, LOPEN, 'NSOL_K', 'F', 1, IDUM, &
          & NSOL_K, CDUM, LDUM, 'Solvent concentration (mol/L)')

     !---------------------------------------------------------------
     !T_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(IO, LOPEN, 'T_K', 'F', 1, IDUM, &
          & T_K, CDUM, LDUM, 'Temperature (in K)')

     !---------------------------------------------------------------
     !NION_COUNT_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(IO, LOPEN, 'NION_COUNT_K', 'I', 1, NION_COUNT_K, &
          & RDUM, CDUM, LDUM, 'Number of types of ionic species in the solvent')
     IF (NION_COUNT_K > 2) THEN
        IF (ALLOCATED(NION_K)) DEALLOCATE(NION_K)
        IF (ALLOCATED(ZION_K)) DEALLOCATE(ZION_K)
        IF (ALLOCATED(RION_K)) DEALLOCATE(RION_K)
     ENDIF

     !---------------------------------------------------------------
     !NION_K
     !---------------------------------------------------------------
     IF (.NOT. ALLOCATED(NION_K)) ALLOCATE(NION_K(NION_COUNT_K))
     CALL READ_FROM_INCAR(IO, LOPEN, 'NION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Ionic concentrations (mol/L)', NION_K)

     !---------------------------------------------------------------
     !ZION_K
     !---------------------------------------------------------------
     IF (.NOT. ALLOCATED(ZION_K)) ALLOCATE(ZION_K(NION_COUNT_K))
     CALL READ_FROM_INCAR(IO, LOPEN, 'ZION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Charge on each ionic species (dimensionless)', ZION_K)

     !---------------------------------------------------------------
     !RION_K
     !---------------------------------------------------------------
     IF (.NOT. ALLOCATED(RION_K)) ALLOCATE(RION_K(NION_COUNT_K))
     CALL READ_FROM_INCAR(IO, LOPEN, 'RION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Ionic radius of  each ionic species (Angstroms)', RION_K)

  ! override kappa^2 for nonlinear sol
  ! in 1/a.u^2
     KAPPA2B_K = 4._q * PI * SUM(NION_K * 6.022 * 1E-4_q * (AUTOA**3) * ZION_K * ZION_K) / (BOLKEV * T_k/(2 * RYTOEV))/ EB_K
  
  ENDIF !!! end of if (lnlnsol)

ENDIF !!! end of lsol


  RETURN

  END SUBROUTINE SOL_READER



!******************** SUBROUTINE SOL_WRITER ***************************
!      
!
! writes the solvation model parameters to the OUTCAR file
!
!
!***********************************************************************
  SUBROUTINE SOL_WRITER(IO)

    USE base

    IMPLICIT NONE

    TYPE (in_struct), INTENT(in) :: IO

    ! early exit if possible
    IF (.NOT.LSOL) RETURN

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,100) LNLSOL, sigma_k, nc_k, T_k, Nsol_k, eb_k, einfty_k, tau
    ENDIF

100   FORMAT( &
           ' Solvation parameters'/ &
           '   LNLSOL   =',L10,  '    Nonlinear solvation (T or F)' /&
           '   SIGMA_K  =',F10.6,  '  width of the dielectric cavity' /&
           '   NC_K     =',F10.6,  '  cutoff charge density (Angstrom^-3)' /&
           ' Parameters that describe the solvent environement'/ &
           '   T_K      =',F10.6,  '  Temperature (K)' /&
           '   Nsol_K   =',F10.6,  '  molarity (mol/L)' /&
           '   EB_K     =',F10.6,  '  relative permittivity of the bulk solvent' /&
           '   EINFTY_K =',F10.6,  '  relative permittivity(optical) of the bulk solvent' /&
           '   P0_K     =',F10.6,  '  dipole moment of the independent dipole, computed from eb_k and einfty_k' /&
           '   TAU      =',F10.6,  '  cavity surface tension (eV/Angstrom^2)' /)

    RETURN

  END SUBROUTINE SOL_WRITER



!******************** SUBROUTINE XML_WRITE_SOL *************************
!      
!
! writes the solvation model parameters to vasprun.xml
!
!
!***********************************************************************
  SUBROUTINE XML_WRITE_SOL

    USE vaspxml

    IMPLICIT NONE

    INTEGER IDUM
    REAL(q) RDUM
    COMPLEX(q) CDUM
    LOGICAL LDUM
    CHARACTER (1) CHARAC

    CALL XML_TAG("separator","solvation model")

    CALL XML_INCAR('LSOL','L',IDUM,RDUM,CDUM,LSOL,CHARAC,1)

    IF (LSOL) THEN
       CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,1)
       CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,1)
       CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,1)
       CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,1)
       CALL XML_INCAR('LAMBDA_D_K','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,1)
       CALL XML_INCAR('LRHOB','L',IDUM,RDUM,CDUM,LRHOB,CHARAC,1)
       CALL XML_INCAR('LRHOION','L',IDUM,RDUM,CDUM,LRHOION,CHARAC,1)
       CALL XML_INCAR('LNLSOL','L',IDUM,RDUM,CDUM,LNLSOL,CHARAC,1)
       IF (LNLSOL) THEN
          CALL XML_INCAR('P0_K','F',IDUM,p0_k,CDUM,LDUM,CHARAC,1)
          CALL XML_INCAR('EINFTY_K','F',IDUM,einfty_k,CDUM,LDUM,CHARAC,1)
          CALL XML_INCAR('NSOL_K','F',IDUM,Nsol_k,CDUM,LDUM,CHARAC,1)
          CALL XML_INCAR('T_K','F',IDUM,T_k,CDUM,LDUM,CHARAC,1)
       ENDIF
    ENDIF

    CALL XML_CLOSE_TAG

    RETURN

  END SUBROUTINE XML_WRITE_SOL



!******************** SUBROUTINE SOL_VCORRECTION ***********************
!
!
! Computes the potential, energy and force corrections due to solvation
!
!
!***********************************************************************
  SUBROUTINE SOL_VCORRECTION(INFO, T_INFO, LATT_CUR, P, WDES, GRIDC, &
       & CHTOT, CVTOT)

    USE base
    USE poscar
    USE lattice
    USE pseudo
    USE mgrid
    USE wave
    USE mdipol

    IMPLICIT NONE

    TYPE (info_struct), INTENT(in) :: INFO
    TYPE (type_info), INTENT(in) :: T_INFO
    TYPE (latt), INTENT(IN) :: LATT_CUR
    TYPE (potcar), INTENT(IN) :: P(T_INFO%NTYP)
    TYPE (wavedes), INTENT(IN) :: WDES
    TYPE (grid_3d), INTENT(IN) :: GRIDC
    
    COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ)
    COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
    
! local variables
    REAL(q) :: Ecorr, EDIFF
    COMPLEX(q), ALLOCATABLE :: Vcorr(:)
    INTEGER :: I

! early exit if possible
    IF (.NOT.LSOL) RETURN

    ALLOCATE(Vcorr(GRIDC%MPLWV))
    
    CALL SOL_INIT(IO_1, EDIFF_1, INFO, LATT_CUR, GRIDC, T_INFO, P, &
      & sigma_k, nc_k, eb_k, tau, lambda_d_k, kappa2b_k, &
      & LNLSOL, p0_k, einfty_k, Nsol_k, T_k, &
      & NION_COUNT_K, NION_K, ZION_K, RION_K, &
      & EDIFFSOL, LRHOB, LRHOION)
         
    CALL VCORRECTION(GRIDC, P, LATT_CUR, T_INFO, WDES, &
         & CHTOT, Vcorr, Ecorr, EIFOR_SOL) 

    DO I = 1, GRIDC%RC%NP
       CVTOT(I,1) = CVTOT(I,1) + Vcorr(I)
    ENDDO

    Ediel_sol = Ecorr

!setting the value of epsilon used in computing the dipole and monopole corrections            
!to EB_K
    IF (((DIP%IDIPCO>0) .OR. (DIP%LMONO)) ) THEN
#ifdef debugsol
       IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'monopole or dipole/quadrapole energy corrections set in INCAR'
       IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Setting DIP%EPSILON = EB_K, the bulk dielectric constant of the solvent, since solvent is included in the vacuum region'
       IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Current value of DIP%EPSILON = ', DIP%EPSILON
#endif

       DIP%EPSILON = EB_K

#ifdef debugsol
       IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'New value of DIP%EPSILON = ', DIP%EPSILON
#endif

    ENDIF

    DEALLOCATE(Vcorr)
      
    RETURN

  END SUBROUTINE SOL_VCORRECTION



END MODULE solvation
