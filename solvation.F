#include "symbol.inc"

MODULE solvation
  USE prec
  USE base

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SOL_READER, SOL_WRITER, SOL_Vcorrection

  LOGICAL, SAVE :: LSOL=.FALSE.

!energy and force
  REAL(q), PUBLIC, SAVE :: Ediel_SOL=0._q
  REAL(q), PUBLIC, ALLOCATABLE, SAVE :: EIFOR_SOL(:,:)
!solvation params
  REAL(q), PRIVATE, SAVE  :: sigma_k = 0.6_q
  REAL(q), PRIVATE, SAVE  :: nc_k = 0.0025_q !0.00473_q
  REAL(q), PRIVATE, SAVE  :: eb_k = 80.0_q
  REAL(q), PRIVATE, SAVE  :: tau = 5.25E-4_q

  TYPE (in_struct), PRIVATE, SAVE :: IO_1
  REAL(q), PRIVATE, SAVE :: EDIFF_1


CONTAINS



!******************** SUBROUTINE SOL_READER ***************************************
!
!
! Reads in the solvation model parameters
!
!      
!**********************************************************************************

  SUBROUTINE SOL_READER(NIONS,EDIFF,IO)

    USE base
    USE pot_k, only : SOL_INIT
    
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

    EIFOR_SOL=0._q; Ediel_SOL=0._q
    IO_1 = IO
    EDIFF_1 = EDIFF

      LOPEN = .FALSE.
   
      OPEN(UNIT=IO%IU5,FILE='INCAR',STATUS='OLD')      

! read in flag LSOL
    CALL RDATAB(LOPEN,'INCAR',IO%IU5,'LSOL','=','#',';','L', &
         &            IDUM,RDUM,CDUM,LSOL,CHARAC,N,1,IERR)
    IF (((IERR/=0).AND.(IERR/=3)).OR. &
         &                    ((IERR==0).AND.(N<1))) THEN
       IF (IO%IU0>=0) &
            WRITE(IO%IU0,*)'Error reading item ''LSOL'' from file INCAR.'
       STOP
    ENDIF

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

!---------------------------------------------------------------------------
!nc_k
!--------------------------------------------------------------------------
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

!---------------------------------------------------------------------------
!eb_k
!---------------------------------------------------------------------------
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

!---------------------------------------------------------------------------
!tau
!--------------------------------------------------------------------------
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

    ENDIF
      
    RETURN

  END SUBROUTINE SOL_READER



!******************** SUBROUTINE SOL_WRITER ***************************************
!      
!
! writes the solvation model paramters to the OUTCAR file
!
!
!**********************************************************************************
  SUBROUTINE SOL_WRITER(IO)

    USE base

    IMPLICIT NONE

    TYPE (in_struct), INTENT(in) :: IO

    ! early exit if possible
    IF (.NOT.LSOL) RETURN

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,100) eb_k, sigma_k, nc_k, tau
    ENDIF

100 FORMAT( &
             ' Solvation parameters'/ &
             '   EB_K    =',F10.6,  '     relative permittivity of the bulk solvent' /&
             '   SIGMA_K =',F10.6,  '     width of the dielectric cavity' /&
             '   NC_K    =',F10.6,  '     cutoff charge density' /&
             '   TAU     =',F10.6,  '     cavity surface tension' /)


    RETURN

  END SUBROUTINE SOL_WRITER



!******************** SUBROUTINE SOL_Vcorrection *********************************
!
!
! Computes the potential, energy and force corrections due to solvation
!
!
!********************************************************************************

  SUBROUTINE SOL_Vcorrection(INFO, T_INFO, LATT_CUR, P, WDES, GRIDC, CHTOT, CVTOT)

    USE base
    USE poscar
    USE lattice
    USE pseudo
    USE mgrid
    USE wave
    USE mdipol
    USE pot_k, ONLY : SOL_INIT, VCORRECTION

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
    
    CALL SOL_INIT(IO_1, EDIFF_1, INFO, LATT_CUR, GRIDC, T_INFO, P, sigma_k, nc_k, eb_k, tau)
    CALL VCORRECTION(GRIDC, P, LATT_CUR, T_INFO, WDES, CHTOT, Vcorr, Ecorr, EIFOR_SOL) 

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

  END SUBROUTINE SOL_Vcorrection
  

END MODULE solvation
