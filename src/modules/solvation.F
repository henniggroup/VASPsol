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
  LOGICAL, SAVE :: LRHOB=.FALSE.
  LOGICAL, SAVE :: LRHOION=.FALSE.

!energy and force
  REAL(q), PUBLIC, SAVE :: Ediel_SOL = 0._q
  REAL(q), PUBLIC, ALLOCATABLE, SAVE :: EIFOR_SOL(:,:)
!solvation params
  REAL(q), PRIVATE, SAVE  :: sigma_k = 0.6_q
  REAL(q), PRIVATE, SAVE  :: nc_k = 0.0025_q !! 0.00473_q
  REAL(q), PRIVATE, SAVE  :: eb_k = 78.4_q
  REAL(q), PRIVATE, SAVE  :: tau = 5.25E-4_q
  !kappa2b_k =  (1/debye screeing lenth in ang)**2
  ! to convert to au, mulitply by 0.5292**2
  !debye screening length = 3 Ang ==> (1./3^2) * (0.5292)^2 = 0.03
  ! To override set lambda_d_k in the INCAR file
  REAL(q), PRIVATE, SAVE  :: lambda_d_k = 1E100_q ! some large value
  REAL(q), PRIVATE, SAVE  :: kappa2b_k = 0._q 
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

    EIFOR_SOL = 0._q
    Ediel_SOL = 0._q
    IO_1 = IO
    EDIFF_1 = EDIFF
    EDIFFSOL = EDIFF/100._q

    LOPEN = .FALSE.
   
    OPEN(UNIT=IO%IU5,FILE='INCAR',STATUS='OLD')      

!-----------------------------------------------------------------
! LSOL, solvation on/off
!-----------------------------------------------------------------
    CALL RDATAB(LOPEN,'INCAR',IO%IU5,'LSOL','=','#',';','L', &
         &            IDUM,RDUM,CDUM,LSOL,CHARAC,N,1,IERR)
    IF (((IERR/=0).AND.(IERR/=3)).OR. &
         &                    ((IERR==0).AND.(N<1))) THEN
       IF (IO%IU0>=0) &
            WRITE(IO%IU0,*)'Error reading item ''LSOL'' from file INCAR.'
       STOP
    ENDIF
    CALL XML_INCAR('LSOL','L',IDUM,RDUM,CDUM,LSOL,CHARAC,N)

    IF (LSOL) THEN

!-----------------------------------------------------------------
!sigma_k
!-----------------------------------------------------------------
       CALL RDATAB(LOPEN,'INCAR',IO%IU5,'SIGMA_K','=','#',';','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
       IF (IERR == 3 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'''SIGMA_K'',the width of dielectric cavity, not specified, using the default(fit parameter for solvent = water)'
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'sigma_k = ',sigma_k          
       ELSE IF (IERR == 0 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Found : sigma_k(width of dielectric cavity) = ',sigma_k 
         
       ENDIF
#endif

       IF ((IERR/=0).AND.(IERR/=3)) THEN
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Error reading ''SIGMA_K''from INCAR.'
          IF (IO%IU0>=0) WRITE(IO%IU0,*) ' Error code = ', IERR
          STOP
       ENDIF
       CALL XML_INCAR('SIGMA_K','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N)

!----------------------------------------------------------------------
!nc_k
!----------------------------------------------------------------------
       CALL RDATAB(LOPEN,'INCAR',IO%IU5,'NC_K','=','#',';','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
       IF (IERR == 3 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'''NC_K'', cavity turn-on charge density, not specified, using the default(fit parameter for solvent = water)'

          IF (IO%IU0>=0) WRITE(IO%IU0,*)'nc_k = ', nc_k          
       ELSE IF (IERR == 0 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Found : nc_k (cavity turn-on charge density) = ', nc_k 
         
       ENDIF
#endif

       IF ((IERR/=0).AND.(IERR/=3)) THEN
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Error reading ''NC_K''from INCAR.'
          IF (IO%IU0>=0) WRITE(IO%IU0,*) ' Error code = ', IERR
          STOP
       ENDIF
       CALL XML_INCAR('NC_K','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N)

!-----------------------------------------------------------------------
!eb_k
!-----------------------------------------------------------------------
       CALL RDATAB(LOPEN,'INCAR',IO%IU5,'EB_K','=','#',';','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
       IF (IERR == 3 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'''EB_K'', relative permittivity of the bulk solvent, not specified, using the default(for solvent = water)'

          IF (IO%IU0>=0) WRITE(IO%IU0,*)'eb_k = ', eb_k          
       ELSE IF (IERR == 0 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Found : eb_k (relative permittivity of the bulk solvent) = ', eb_k 
         
       ENDIF
#endif

       IF ((IERR/=0).AND.(IERR/=3)) THEN
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Error reading ''EB_K''from INCAR.'
          IF (IO%IU0>=0) WRITE(IO%IU0,*) ' Error code = ', IERR
          STOP
       ENDIF
       CALL XML_INCAR('EB_K','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N)

!-----------------------------------------------------------------------
!tau
!-----------------------------------------------------------------------
       CALL RDATAB(LOPEN,'INCAR',IO%IU5,'TAU','=','#',';','F',IDUM,tau,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
       IF (IERR == 3 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'''TAU'',cavity surface tension, not specified, using the default'
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'tau = ',tau
         
       ELSE IF (IERR == 0 ) THEN
      
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Found : tau(cavity surface tension) = ',tau
         
       ENDIF
#endif

       IF ((IERR/=0).AND.(IERR/=3)) THEN
          IF (IO%IU0>=0) WRITE(IO%IU0,*)'Error reading ''TAU''from INCAR.'
          IF (IO%IU0>=0) WRITE(IO%IU0,*) ' Error code = ', IERR
          STOP
       ENDIF
       CALL XML_INCAR('TAU','F',IDUM,tau,CDUM,LDUM,CHARAC,N)

    ENDIF

  !-------------------------------------------------------------------
  !lambda_d_k = debye screening length in angstrom : the user input from INCAR.
  !kappa2b_k = (1/debye screening length)^2 
  ! if lambd_d_k not set in INCAR, kappa2b_k will be set to zero(default value)
  !--------------------------------------------------------------------
  CALL RDATAB(LOPEN,'INCAR',IO%IU5,'LAMBDA_D_K','=','#',';','F',IDUM,lambda_d_k,CDUM,LDUM,CHARAC,N,1,IERR)

  IF (IERR == 3 ) THEN
#ifdef debugsol      
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) '''LAMBDA_D_K'', debye screening length not specified'

     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'kappa2b_k  = (1/ debye screening length^2) = ', kappa2b_k         
#endif         
  ELSE IF (IERR == 0 ) THEN
#ifdef debugsol
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Found : lambda_d_k (debye screening length) = ', lambda_d_k
     
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'lambda_d_k  = (debye screening length) = ', lambda_d_k
#endif

     kappa2b_k = 1._q / (lambda_d_k ** 2)

#ifdef debugsol
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'kappa2b_k  = (1/ debye screening length^2) = ', kappa2b_k
#endif         
  ENDIF

  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Error reading ''LAMBDA_D_K''from INCAR.'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF

  !---------------------------------------------------------------
  !EDIFFSOL
  !---------------------------------------------------------------
  CALL RDATAB(LOPEN,'INCAR',IO%IU5,'EDIFFSOL','=','#',';','F',IDUM,EDIFFSOL,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
  IF (IERR == 3 ) THEN
      
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) '''EDIFFSOL'', solvation convergence tolerance, not specified, using the default(EDIFF/100)'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'EDIFFSOL = ', EDIFFSOL          
         
  ELSE IF (IERR == 0 ) THEN
      
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Found : EDIFFSOL = ', EDIFFSOL 
         
  ENDIF
#endif

  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Error reading ''EDIFFSOL''from INCAR.'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF

  !---------------------------------------------------------------
  !LRHOB
  !---------------------------------------------------------------
  CALL RDATAB(LOPEN,'INCAR',IO%IU5,'LRHOB','=','#',';','L',IDUM, &
       &  RDUM,CDUM, LRHOB, CHARAC,N,1,IERR)
#ifdef debugsol
  IF (IERR == 3 ) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) '''LRHOB'', write the bound charge density or not'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'LRHOB = ', LRHOB                  
  ELSE IF (IERR == 0 ) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Found : LRHOB = ', LRHOB
  ENDIF
#endif
  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Error reading ''LRHOB''from INCAR.'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF

  !---------------------------------------------------------------
  !LRHOION
  !---------------------------------------------------------------
  CALL RDATAB(LOPEN,'INCAR',IO%IU5,'LRHOION','=','#',';','L',IDUM, &
       &  RDUM,CDUM, LRHOION, CHARAC,N,1,IERR)
#ifdef debugsol
  IF (IERR == 3 ) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) '''LRHOION'', write the ionic charge density or not'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'LRHOION = ', LRHOION                 
  ELSE IF (IERR == 0 ) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Found : LRHOION = ', LRHOION
  ENDIF
#endif
  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) 'Error reading ''LRHOION''from INCAR.'
     IF (IO_1%IU0>=0) WRITE(IO_1%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF
      
  CLOSE(IO%IU5)

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
       WRITE(IO%IU6,100) LSOL, eb_k, sigma_k, nc_k, tau, &
            & lambda_d_k, LRHOB, LRHOION
    ENDIF

100 FORMAT( &
             ' Solvation model:'/ &
             '   LSOL    =',L6, '         switch on solvation model'/ &
             '   EB_K    =',F10.6,  '     relative permittivity of the bulk solvent' /&
             '   SIGMA_K =',F10.6,  '     width of the dielectric cavity' /&
             '   NC_K    =',F10.6,  '     cutoff charge density' /&
             '   TAU     =',F10.6,  '     cavity surface tension' /&
             '   LAMBDA_D_K    =',F10.6,  '     Debye length in Angstroms' /&
             '   LRHOB     =',F10.6,  '     write boundcharge  density' /&
             '   LRHOION     =',F10.6,  '     write ioniccharge density' /)


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
