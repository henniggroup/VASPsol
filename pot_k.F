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
#define sol_io  IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)

!>==========================MODULE POT_K===============================
!!
!! This module glues together the solvation engine module POT_LPCM_K
!!(defined in pot_lpcm_cav_k.F) with the rest of the VASP code
!!
!! Contains subroutines that perform solvation I/O, pre-process 
!!the charge densities, call the main subroutine, Vcorrection_lpcm,  
!!from the module POT_LPCM_K and compute the force corrections
!!
!!=====================================================================
MODULE POT_K

 USE prec
 USE base
 USE mpimy 
 USE pot_lpcm_k, ONLY : Vcorrection_lpcm, MY_D_PROD, SP_IO, SimplePoisson
 USE pot_lpcm_k, ONLY : GET_VAC_LEVEL, VACPOT_POIS, UNPACK_C2R
#ifdef debugsolwf
 USE pot_lpcm_k, ONLY : WRITE_TO_FILE
#endif

#ifdef jdftx
 EXTERNAL initjdftx, getionsigma, minimizefluid
#endif

 PUBLIC :: SOL_INIT, Vcorrection
 ! vaclevel
#ifdef debugsol_pbz
  PUBLIC :: GET_FERMISHIFT
#endif
  ! vaclevel

 PRIVATE :: IONIC_CHARGE, GETRCORE, Fcorrection, JDFTX_INTERFACE

 INTEGER, PRIVATE, SAVE :: PK_counter = 1
 ! numerical params
 REAL(q), PRIVATE, SAVE :: sigma_nc_k
 REAL(q), ALLOCATABLE, PRIVATE, SAVE :: sigma_rc_k(:), CORE_C(:)
 ! empirical params
 REAL(q), PRIVATE, SAVE :: sigma_k, nc_k, eb_k
 REAL(q), PRIVATE, SAVE :: lambda_d_k, kappa2b_k, tau, EDIFFSOL

 LOGICAL, PRIVATE, SAVE :: LSOL = .FALSE.
 LOGICAL, PRIVATE, SAVE :: LRHOB = .FALSE., LRHOION = .FALSE.

 LOGICAL, SAVE :: LDEFAULTPCM = .FALSE.
 LOGICAL, SAVE :: LJDFTX = .FALSE.

 CONTAINS

!>==========================SUBROUTINE SOL_INIT======================
!!
!! Initialize the module parameters
!! the numerical parameters(sigma_nc_k, sigma_rc_k, core_c) are read 
!! from INCAR if specified(for testing prpose)
!! and initiaize the grid for jdftx(if LJDFTX flag is set in INCAR)
!! 
!!======================================================================
 SUBROUTINE SOL_INIT(IO, EDIFF, INFO, LATT_CUR, GRIDC, T_INFO, P, &
      & sigma_k_1, nc_k_1, eb_k_1, tau_1, lambda_d_k_1, kappa2b_k_1, &
      & EDIFFSOL_1, LRHOB_1, LRHOION_1)

   USE poscar, ONLY : type_info 
   USE lattice, ONLY : latt, DIRKAR 
   USE mgrid, ONLY : grid_3d
   USE pseudo, ONLY : potcar

   IMPLICIT NONE

   TYPE(type_info),INTENT(IN) :: T_INFO
   TYPE(potcar),INTENT(IN) :: P(T_INFO%NTYP)
   TYPE (info_struct) :: INFO
   TYPE (latt),INTENT(IN) :: LATT_CUR
   TYPE(grid_3d),INTENT(IN) :: GRIDC
   TYPE (in_struct),INTENT(IN) :: IO

   REAL(q) :: EDIFF, EDIFFSOL_1
   REAL(q) :: sigma_gauss, RCS(T_INFO%NTYP)
   REAL(q) :: sigma_k_1, nc_k_1, eb_k_1, tau_1
   REAL(q) :: lambda_d_k_1, kappa2b_k_1
   INTEGER :: IDUM, N, IERR, K, I, II, ATOMIC_NUM
   REAL(q) :: RDUM,delta_grd,sigma_rc_k0(T_INFO%NTYP)
   COMPLEX(q) CDUM
   LOGICAL :: LRHOB_1, LRHOION_1
   LOGICAL :: LOPEN,LDUM,LIONODE_FLAG,LFLUEX
   CHARACTER (1) CHARAC
   
   LIONODE_FLAG = .FALSE.
 
   SP_IO = IO
   LSOL = .TRUE. 

   ! set the module params
   sigma_k = sigma_k_1
   nc_k = nc_k_1
   eb_k = eb_k_1
   tau = tau_1
   lambda_d_k = lambda_d_k_1
   kappa2b_k = kappa2b_k_1
   LRHOB = LRHOB_1
   LRHOION = LRHOION_1
   EDIFFSOL = EDIFFSOL_1
   
#ifdef MPI
   IF (GRIDC%COMM%NODE_ME==GRIDC%COMM%IONODE) THEN
#endif
      LIONODE_FLAG = .TRUE.
#ifdef MPI
   ENDIF
#endif   

 IF (PK_COUNTER == 1) THEN
    ALLOCATE(sigma_rc_k(T_INFO%NTYP))
    ALLOCATE(CORE_C(T_INFO%NTYP))
#ifdef debugsol
    sol_io 'Setting up the paramters.'
    sol_io 'Unless  explicity specified in the INCAR file these values will be used as the defaults'
#endif
    ! delta_grd = ((LATT_CUR%OMEGA)/(GRIDC%NPLWV))**(1._q/3._q)
    delta_grd = MAXVAL(LATT_CUR%ANORM/GRIDC%NGPTAR)
#ifdef debugsol      
    sol_io 'delta_grd = ', delta_grd      
#endif

    !---------------------------------------------------------------
    !Defaults values of sigma_nc_k(sigma for the gaussian nuclei), 
    !sigma_rc_k(sigma for the pseudo core)
    !and CORE_C
    !---------------------------------------------------------------
    sigma_nc_k = 1.6_q * delta_grd !sigma for point charges

    CALL GETRCORE(T_INFO%NTYP,RCS)

    !PSDMAX, maximal r for nonlocal contribution 
    !(in fact rmax=PSRMAX/NPSNL*(NPSNL-1))
    ! RCS = (/(P(K)%PSRMAX,K=1,T_INFO%NTYP)/) 
    ! sigma_rc_k = MINVAL(RCS)*2.0/4.0_q 
#ifdef debugsol
    sol_io 'RCOREs = ', RCS      
#endif
    sigma_rc_k = RCS/2.5_q   ! sigma for the pseudo core
    DO II = 1, T_INFO%NTYP
       ATOMIC_NUM = INT(P(II)%ZVALF_ORIG+P(II)%ZCORE)
       IF (ATOMIC_NUM <= 10) THEN
#ifdef debugsol
          IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'ATOMIC NUMBER, SYMBOL, VALENCE, CORE :',ATOMIC_NUM, P(II)%ELEMENT, P(II)%ZVALF_ORIG, P(II)%ZCORE
          IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Setting CORE_C(the magnitude of the pseudo core charge) to zero.'
#endif
          CORE_C(II) = 0._q
          IF (ATOMIC_NUM <= 4 .AND. P(II)%ZVALF_ORIG /= ATOMIC_NUM) THEN
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'!------------------------------------------------------------------------'
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'WARNING: '
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Use all-electron pseudopotential or set CORE_C explicitly '
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'for atomic-species with atomic number less than or equal to 4.'
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Problem species ATOMIC NUMBER, SYMBOL & VALENCE :',ATOMIC_NUM, P(II)%ELEMENT, P(II)%ZVALF_ORIG
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'I HOPE YOU KNOW WHAT YOU ARE DOING !'
             IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'!------------------------------------------------------------------------'
          ENDIF
       ELSE
          CORE_C(II) = 1._q
       ENDIF
    ENDDO

    !------------------------------------------------------------------
    !the afore-mentioned three paramters are essential for all pcm models
    ! Try to read those from the INCAR file
    !-------------------------------------------------------------------
    LOPEN = .FALSE.
   
    OPEN(UNIT=SP_IO%IU5,FILE='INCAR',STATUS='OLD')      

    !---------------------------------------------------------------
    !sigma_nc_k
    !---------------------------------------------------------------
    CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'SIGMA_NC_K','=','#',';','F',IDUM,sigma_nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol      
      IF (IERR == 3 ) THEN
      
         sol_io '''SIGMA_NC_K'',the width of ionic gaussians, not specified, using the default'
         sol_io 'sigma_nc_k = ',sigma_nc_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         sol_io 'Found : sigma_nc_k(width of ionic gaussians) = ',sigma_nc_k 

      ENDIF
#endif
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         sol_io 'Error reading ''SIGMA_NC_K''from INCAR.'
         sol_io ' Error code = ', IERR
         STOP
      ENDIF

      !---------------------------------------------------------------
      !sigma_rc_k
      !---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'SIGMA_RC_K','=','#',';','F',IDUM,sigma_rc_k0,CDUM,LDUM,CHARAC,N,T_INFO%NTYP,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         sol_io '''SIGMA_RC_K'',the width of pseudo-core gaussians, not specified, using the default'
         sol_io 'sigma_rc_k = ',sigma_rc_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         sol_io 'Found : sigma_rc_k(width of pseudo-core gaussians) = ',sigma_rc_k0
         DO I=1,T_INFO%NTYP
            sigma_rc_k(I) = sigma_rc_k0(I)
         ENDDO
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         sol_io 'Error reading ''SIGMA_RC_K''from INCAR.'
         sol_io  ' Error code = ', IERR
         STOP
      ENDIF

      !---------------------------------------------------------------
      !CORE_C
      !---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'CORE_C','=','#',';','F',IDUM,CORE_C,CDUM,LDUM,CHARAC,N,T_INFO%NTYP,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         sol_io '''CORE_C'', magnitude of pseudo-core charge, not specified, using the default'

         sol_io 'CORE_C = ', CORE_C          
         
      ELSE IF (IERR == 0 ) THEN
      
         sol_io 'Found : CORE_C (magnitude of pseudo-core charge) = ', CORE_C 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         sol_io 'Error reading ''CORE_C''from INCAR.'
         sol_io ' Error code = ', IERR
         STOP
      ENDIF

      !---------------------------------------------------------------
      !LJDFTX
      !---------------------------------------------------------------
#ifdef jdftx   
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'LJDFTX','=','#',';','L',IDUM,RDUM,CDUM,LJDFTX,CHARAC,N,1,IERR)
#endif
!error opening INCAR      
      IF ( IERR == 1 .OR. IERR == 2) THEN
         IF (LIONODE_FLAG) WRITE(SP_IO%IU0,*)'Error reading INCAR'
         STOP
      ENDIF
#ifdef jdftx
      !check for flucar
#ifdef debugsol
      IF (LIONODE_FLAG .AND. LJDFTX) WRITE(*,*)'Checking for FLUCAR file'       
#endif
      INQUIRE(FILE='FLUCAR',EXIST=LFLUEX)

!if LJDFTX is true and flucar exists initialize jdftx
      IF (LJDFTX .AND. LFLUEX) THEN
         
         IF (LIONODE_FLAG) THEN
#ifdef debugsol
            PRINT *, "Initializing jdftx"
#endif
            CALL initjdftx(LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3),GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ)
#ifdef debugsol
            PRINT *,"Setting ionic gaussian width, sigma_nc_k"
#endif
            CALL getionsigma(sigma_gauss)
         ENDIF
         CALLMPI( M_bcast_d(GRIDC%COMM, sigma_gauss, 1)) 
         sigma_nc_k = sigma_gauss
#ifdef debugsol
         IF (LIONODE_FLAG) WRITE(SP_IO%IU0,*)'sigma_nc_k set to ',sigma_nc_k
#endif
!otherwise use the default pcm
      ELSE
#ifdef debugsol
         IF (LIONODE_FLAG) WRITE(SP_IO%IU0,*)'Either LJDFTX = .FALSE. or FLUCAR file not present'
         IF (LIONODE_FLAG) WRITE(SP_IO%IU0,*)'using the default linear pcm in vasp'
#endif
#endif
         LDEFAULTPCM = .TRUE.
#ifdef debugsol
         IF (LIONODE_FLAG) WRITE(SP_IO%IU0,'(/A/)') 'All solvation model paramters are set'
#endif
#ifdef jdftx        
      ENDIF
#endif
#ifdef debugsol
      IF (LIONODE_FLAG) WRITE(SP_IO%IU0,*)'Solvation model intialization DONE'
#endif

   ENDIF

   CLOSE(SP_IO%IU5)

   RETURN

 END SUBROUTINE SOL_INIT

 
!>==========================SUBROUTINE Vcorrection=====================
!!
!! Creates PS_TOTN(n_val+pseudo_core) and TOTN(n_val+N_gauss)
!! Passes them on to vcorrection_lpcm in pot_lpcm_cav_k.F or to 
!! the jdftx_interface
!! Uses the so obtained Vcorr, Vdiel and Atot to compute the net 
!! correction to the local potential, Vcorr(not the same as above) 
!! and the net correction to the total energy, Ecorr
!!
!! called from pot.F
!!
!!=====================================================================
 SUBROUTINE Vcorrection(SP_GRIDC, SP_P, SP_LATT_CUR, SP_T_INFO, WDES,&
                      & SP_CHTOT, Vcorr, Ecorr, EIFOR_SOL)

  USE constant 
  USE poscar, ONLY : type_info 
  USE lattice, ONLY : latt, DIRKAR 
  USE mgrid, ONLY : grid_3d
  USE pseudo, ONLY : potcar
  USE wave, ONLY : wavedes

  IMPLICIT NONE

  TYPE(type_info),INTENT(IN) :: SP_T_INFO
  TYPE(potcar),INTENT(IN) :: SP_P (SP_T_INFO%NTYP)
  TYPE (latt),INTENT(IN) :: SP_LATT_CUR
  TYPE(grid_3d),INTENT(IN) :: SP_GRIDC
  TYPE (wavedes) ::    WDES

  REAL(q) :: Atot
  COMPLEX(q), INTENT(OUT) ::  Vcorr(SP_GRIDC%MPLWV)
  REAL(q) :: EIFOR_SOL(3, SP_T_INFO%NIONS)
  COMPLEX(q) :: SP_CHTOT(SP_GRIDC%MPLWV, WDES%NCDIJ)
  COMPLEX(q),ALLOCATABLE :: SP_N(:),TOTN(:),SP_N_ION(:)
  COMPLEX(q),ALLOCATABLE :: Vdiel(:)
  COMPLEX(q),ALLOCATABLE :: PS_TOTN(:)
  ! electrost. pot. diff * n   and Vdiel * n respectively
  REAL(q) :: Ecorr3, Ecorr4 
  REAL(q),INTENT(OUT) :: Ecorr ! Atot - Ecorr3 - Ecorr4
  INTEGER(q) :: I,J,K
  REAL(q)  :: dx,dy,dz,dv

  IF(LSOL) THEN

    ALLOCATE(SP_N(SP_GRIDC%MPLWV),TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(Vdiel(SP_GRIDC%MPLWV))
    ALLOCATE(PS_TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(SP_N_ION(SP_GRIDC%MPLWV))
  
    Vdiel = 0._q
    Ecorr3 = 0._q
    Ecorr4 = 0._q
    Atot = 0._q
    EIFOR_SOL = 0._q

#ifdef debugsol
    sol_io 'Begin Solvation calculations'
    sol_io 'Setting up the charge densities'
#endif
!-------------------------------------------------------------
!Gaussian ionic charge in fourier space
!-------------------------------------------------------------
    CALL IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,1)

    SP_N_ION = SP_N

!-------------------------------------------------------------
!total charge in fourier space
!-------------------------------------------------------------
    CALL RC_ADD(SP_CHTOT(1,1),1.0_q,SP_N_ION(1),1.0_q,TOTN(1),SP_GRIDC)
    CALL SETUNB(TOTN,SP_GRIDC)
!    TOTN = SP_CHTOT + SP_N_ION

!-------------------------------------------------------------
!CHARGE BALL fourier, flag=3
!-------------------------------------------------------------
    CALL IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,3)

!-------------------------------------------------------------
!total charge in k space,valence+pseudo core
!-------------------------------------------------------------
   CALL RC_ADD(SP_CHTOT(1,1),1.0_q,SP_N(1),1.0_q,PS_TOTN(1),SP_GRIDC)
   CALL SETUNB(PS_TOTN,SP_GRIDC)

!   PS_TOTN = SP_CHTOT + SP_N

#ifdef debugsol
   sol_io 'Begin minimization'
#endif

!-------------------------------------------------------------
!
!              BEGIN Solvation minimization
!
!           Compute Vcorr, the correction potential
!           and Atot, the total solvation energy
!
!-------------------------------------------------------------    
   IF (LDEFAULTPCM) THEN

#ifdef debugsol
      sol_io 'Calling vcorr_lpcm minimization'
#endif

!-------------------------------------------------------------
!
!                   VCORRECTION_LPCM
!
!-------------------------------------------------------------
      CALL Vcorrection_lpcm(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO,&
           & SP_CHTOT, TOTN, PS_TOTN, Vcorr, Vdiel, Atot, sigma_nc_k,&
           & sigma_k, nc_k, eb_k, kappa2b_k, CORE_C, sigma_rc_k,&
           & tau, EDIFFSOL, LRHOB, LRHOION)

#ifdef debugsol
      sol_io 'Return from vcorr_lpcm minimization'
#endif

!-------------------------------------------------------------
!JDFTx interface
!-------------------------------------------------------------
#ifdef jdftx
   ELSE

      CALL JDFTX_INTERFACE(SP_GRIDC, SP_LATT_CUR, TOTN, PS_TOTN, &
           & Vdiel, Vcorr, Atot)
#endif

   ENDIF
!-------------------------------------------------------------
!
!              END Solvation minimization
!
!-------------------------------------------------------------    

#ifdef debugsol
   sol_io 'End minimization'
#endif

   CALL SETUNB(Vcorr,SP_GRIDC)
   CALL SETUNB(Vdiel,SP_GRIDC)

!-------------------------------------------------------------
!Ecorr3 = Vcorr*n(just electrostatic correction)
!-------------------------------------------------------------
!test
!   sol_io 'Vcorr(1), electrostatic corr. = ', Vcorr(1)
!test
   CALL MY_D_PROD(Ecorr3, SP_CHTOT(1,1), Vcorr, SP_GRIDC)
!   CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr3,1))
   CALLMPI( M_sum_s(SP_GRIDC%COMM,1,Ecorr3,0,0,0))

!-------------------------------------------------------------
!Ecorr4 = Vdiel*n, 
!where Vdiel = epsilon and kappa2 derivative terms + cavity pot
!-------------------------------------------------------------
!test
!    sol_io 'Vdiel(1), the rest = ', Vdiel(1)
!test
   CALL MY_D_PROD(Ecorr4, SP_CHTOT(1,1), Vdiel, SP_GRIDC)
!   CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr4,1))
   CALLMPI( M_sum_s(SP_GRIDC%COMM,1,Ecorr4,0,0,0))

!-------------------------------------------------------------
!calculate force correction,flag=5(not used)
!-------------------------------------------------------------
   CALL Fcorrection(SP_GRIDC,SP_P,SP_T_INFO,SP_LATT_CUR, &
        & Vcorr,Vdiel, EIFOR_SOL)

!-------------------------------------------------------------
!net correction potential
!Vcorr = electrostatic correction + cavity correction(Vdiel)
!-------------------------------------------------------------
   Vcorr = Vcorr + Vdiel

#ifdef debugsol
   sol_io 'Atot, Ecorr3, Ecorr4 ', Atot, Ecorr3, Ecorr4
#endif

!-------------------------------------------------------------
!net correction to the total energy, E%Ediel = Ecorr
!-------------------------------------------------------------
   Ecorr = Atot - Ecorr3 - Ecorr4

#ifdef debugsol
   IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,'(/a, (es15.6))')"E%Ediel = Ecorr", Ecorr
#endif

!-------------------------------------------------------------
!keeps count of the number of times Vcorrection is called
!-------------------------------------------------------------
   PK_COUNTER = PK_COUNTER + 1

   DEALLOCATE(SP_N,TOTN)
   DEALLOCATE(Vdiel)
   DEALLOCATE(PS_TOTN)
   DEALLOCATE(SP_N_ION)

!-------------------------------------------------------------
!if LSOL = .FALSE., sets Vcorr and Ecorr to 0._q
!-------------------------------------------------------------
  ELSE
   IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0, *) 'LSOL not set. No solvation for you'
   Vcorr = 0._q
   Ecorr = 0._q

  ENDIF

#ifdef debugsol
  IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0, *) 'End Solvation calculations'
#endif

  RETURN

 END SUBROUTINE Vcorrection


!>=====================SUBROUTINE IONIC_CHARGE=========================
!!
!! evaluates the gaussian ionic charge(flag=1) and pseudo-core charge
!!(flag=3)in k space
!!also calculates force corrections(flag=5) in real space
!!=====================================================================
SUBROUTINE IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N, &
     & flag,Vcorr,Vdiel)

 USE constant
 USE mgrid, ONLY : grid_3d
 USE lattice, ONLY : latt,DIRKAR
 USE poscar, ONLY : type_info
 USE pseudo, ONLY : potcar

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: SP_GRIDC
 TYPE (latt),INTENT(IN) :: SP_LATT_CUR
 TYPE (type_info),INTENT(IN) :: SP_T_INFO
 TYPE(potcar),INTENT(IN) :: SP_P (SP_T_INFO%NTYP)

 !ionic charge in k space 
 COMPLEX(q),INTENT(OUT) :: SP_N(SP_GRIDC%MPLWV)
 ! local variables   
 INTEGER,INTENT(IN) :: flag
 ! REAL(q),INTENT(IN) :: sigma_nc_k,sigma_rc_k,CORE_C
 INTEGER :: SP_NC,SP_N1,SP_N2,SP_N3,SP_NT,SP_NIS,SP_NI,SP_NG
 REAL(q) :: dx, dy,dz,dv
 REAL(q) :: SP_G1, SP_G2, SP_G3,SP_GX, SP_GY, SP_GZ, SP_GSQU
 COMPLEX(q) :: SP_CEXPF
 COMPLEX(q) :: SP_FOR_gion, SP_FOR1_gion, SP_FOR2_gion, SP_FOR3_gion
 COMPLEX(q) :: SP_FOR_cb, SP_FOR1_cb, SP_FOR2_cb, SP_FOR3_cb
 COMPLEX(q) :: imag_i
 COMPLEX(q), OPTIONAL, INTENT(IN) :: Vcorr(SP_GRIDC%MPLWV)
 COMPLEX(q), OPTIONAL, INTENT(IN) :: Vdiel(SP_GRIDC%MPLWV)
 INTEGER :: I, J
 
 imag_i = cmplx(0._q,1._q)

 SP_N = 0._q

 SP_NIS = 1  

 !-------------------------------------------------------------
 !loop over all ions
 !-------------------------------------------------------------
 type:DO SP_NT=1,SP_T_INFO%NTYP
    ions: DO SP_NI=SP_NIS,SP_T_INFO%NITYP(SP_NT)+SP_NIS-1

       SP_FOR1_gion=0._q
       SP_FOR2_gion=0._q
       SP_FOR3_gion=0._q

       SP_FOR1_cb=0._q
       SP_FOR2_cb=0._q
       SP_FOR3_cb=0._q

       DO SP_NG=1,SP_GRIDC%RC%NP
          SP_N1= MOD((SP_NG-1),SP_GRIDC%RC%NROW) +1
          SP_NC= (SP_NG-1)/SP_GRIDC%RC%NROW+1
          SP_N2= SP_GRIDC%RC%I2(SP_NC)
          SP_N3= SP_GRIDC%RC%I3(SP_NC)

          SP_G1 = SP_T_INFO%POSION(1,SP_NI)*SP_GRIDC%LPCTX(SP_N1)
          SP_G2 = SP_T_INFO%POSION(2,SP_NI)*SP_GRIDC%LPCTY(SP_N2)
          SP_G3 = SP_T_INFO%POSION(3,SP_NI)*SP_GRIDC%LPCTZ(SP_N3)
          !------------------------------------------------------
          !phase factor
          !------------------------------------------------------
          SP_CEXPF=EXP(-CITPI*(SP_G3+SP_G2+SP_G1))*SP_T_INFO%VCA(SP_NT)
          
          SP_GX = (SP_GRIDC%LPCTX(SP_N1) * SP_LATT_CUR%B(1,1) + &
               & SP_GRIDC%LPCTY(SP_N2) * SP_LATT_CUR%B(1,2) + &
               & SP_GRIDC%LPCTZ(SP_N3) * SP_LATT_CUR%B(1,3))
          SP_GY = (SP_GRIDC%LPCTX(SP_N1) * SP_LATT_CUR%B(2,1) + &
               & SP_GRIDC%LPCTY(SP_N2) * SP_LATT_CUR%B(2,2) + &
               & SP_GRIDC%LPCTZ(SP_N3) * SP_LATT_CUR%B(2,3))
          SP_GZ = (SP_GRIDC%LPCTX(SP_N1) * SP_LATT_CUR%B(3,1) + &
               & SP_GRIDC%LPCTY(SP_N2) * SP_LATT_CUR%B(3,2) + &
               & SP_GRIDC%LPCTZ(SP_N3) * SP_LATT_CUR%B(3,3))
          !------------------------------------------------------
          !G**2
          !------------------------------------------------------
          SP_GSQU = (SP_GX*SP_GX + SP_GY*SP_GY + SP_GZ*SP_GZ) * (TPI**2)
          !------------------------------------------------------
          !flag=1, gaussian ionic charge
          !------------------------------------------------------
          IF ( flag == 1) THEN

             SP_N(SP_NG) = SP_N(SP_NG) - (SP_P(SP_NT)%ZVALF) * &
                  & SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_nc_k**2))
             !-------------------------------------------------------
             !flag=3, pseudo core charge density to repel fluid 
             ! from core
             !-------------------------------------------------------
          ELSE IF ( flag == 3 ) THEN
             
             SP_N(SP_NG) = SP_N(SP_NG) + CORE_C(SP_NT) * SP_CEXPF * &
                  & EXP(-0.5_q*SP_GSQU*(sigma_rc_k(SP_NT)**2))

          ENDIF

       ENDDO !end of loop over reciproc vectors

    ENDDO ions
    
    SP_NIS = SP_NIS + SP_T_INFO%NITYP(SP_NT)

 END DO type

 RETURN

 END SUBROUTINE IONIC_CHARGE


!>==========================SUBROUTINE Fcorrection=====================
!!
!!Force Correction
!!
!!=====================================================================
SUBROUTINE Fcorrection(GRIDC,P,T_INFO,LATT_CUR,Vcorr,Vdiel, EIFOR_SOL)
  USE constant 
  USE poscar, ONLY : type_info 
  USE lattice, ONLY : latt, DIRKAR 
  USE mgrid, ONLY : grid_3d
  USE pseudo, ONLY : potcar

  IMPLICIT NONE

  TYPE (grid_3d)     GRIDC
  TYPE (type_info)   T_INFO
  TYPE (potcar)      P (T_INFO%NTYP)
  TYPE (latt)        LATT_CUR
  
  COMPLEX(q) Vcorr(GRIDC%RC%NP)
  COMPLEX(q) Vdiel(GRIDC%RC%NP)
  INTEGER :: NIS,NT,NIADD,NI,N,N1,NC,N2,N3,NG,NGP,N1P,I,J
  REAL(q) ::G,GX,GY,GZ,G1,G2,G3,FOR1,FOR2,FOR3,FOR,FACTM
  REAL(q) :: EIFOR_SOL(3,T_INFO%NIONS)
  COMPLEX(q) :: CE,CX,CEXPF

  REAL(q), ALLOCATABLE :: GWORK(:),PCWORK(:)
   
  ALLOCATE(GWORK(GRIDC%RC%NP),PCWORK(GRIDC%RC%NP))

  NIS=1

  typ: DO NT=1,T_INFO%NTYP
     
     NIADD=T_INFO%NITYP(NT)

     DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
         
        SETFACT1
        SETFACT
        
        GX = GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*&
             & LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY = GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*&
             & LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ = GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*&
             & LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        GWORK (N) = -(P(NT)%ZVALF)*EXP(-0.5_q*(G**2)*(sigma_nc_k**2))
        PCWORK (N) = (CORE_C(NT))*EXP(-0.5_q*(G**2)*(sigma_rc_k(NT)**2))

     ENDDO

     ion: DO NI=NIS,NIADD+NIS-1
        !-------------------------------------------------------------
        ! initialise the force on the ion to zero
        !-------------------------------------------------------------
        FOR1=0
        FOR2=0
        FOR3=0

        !-------------------------------------------------------------
        !CGXDX,Y,Z = I* the changes in the phase factor g.r on 
        !moving one reciprocal lattice vector in the x,y,z 
        !directions, respectively
        !-------------------------------------------------------------

        !-------------------------------------------------------------
        ! calculate the total force on the ions by summing over 
        !reciprocal lattice vectors
        ! first calculate phase factor:
        ! there are two version for calculating the phase factor
        ! on vector machines you might try the first version
        ! (see stufak.F)
        !-------------------------------------------------------------
#ifdef vector
        DO NG=1,GRIDC%RC%NP
           N1= MOD((NG-1),GRIDC%RC%NROW) +1
           NC= (NG-1)/GRIDC%RC%NROW+1
           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)

           G1=T_INFO%POSION(1,NI)*GRIDC%LPCTX(N1)
           G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
           G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
           SETFACT1
           SETFACT
           CEXPF=EXP(-CITPI*(G3+G2+G1))*T_INFO%VCA(NT)
#else
           CX =EXP(-CITPI*T_INFO%POSION(1,NI))
           G1 =T_INFO%POSION(1,NI)*xmin(GRIDC%NGX)

           DO NC=1,GRIDC%RC%NCOL
              NGP=(NC-1)*GRIDC%RC%NROW+1

              N2= GRIDC%RC%I2(NC)
              N3= GRIDC%RC%I3(NC)
              G2=T_INFO%POSION(2,NI)*GRIDC%LPCTY(N2)
              G3=T_INFO%POSION(3,NI)*GRIDC%LPCTZ(N3)
              CE=EXP(-CITPI*(G3+G2+G1))*T_INFO%VCA(NT)

              DO N1P=0,GRIDC%RC%NROW-1
                 SETN1fromN1P
                 NG=NGP+N1
                 N1=N1+1
                  
                 SETFACT1
                 SETFACT
                 CEXPF=CE
                 CE=CE*CX
#endif

                 !---------------------------------------------------
                 ! add the contribution to the force from the present
                 !  reciprocal lattice
                 ! vector  and multiply by i (ie take imaginary part)
                 !---------------------------------------------------
                 FOR = GWORK(NG)*MULFACT AIMAG(CONJG(Vcorr(NG))*CEXPF)
                 FOR = FOR + PCWORK(NG)*MULFACT AIMAG(CONJG(Vdiel(NG))*CEXPF)
                 FOR1=FOR1-GRIDC%LPCTX_(N1)*FOR
                 FOR2=FOR2-GRIDC%LPCTY_(N2)*FOR
                 FOR3=FOR3-GRIDC%LPCTZ_(N3)*FOR
              ENDDO
#ifndef vector
           ENDDO
#endif

           !---------------------------------------------------------
           ! multiply forces by 2*Pi
           !----------------------------------------------------------
           EIFOR_SOL(1,NI)=FOR1*TPI
           EIFOR_SOL(2,NI)=FOR2*TPI
           EIFOR_SOL(3,NI)=FOR3*TPI
           
        ENDDO ion
        NIS=NIS+NIADD
     ENDDO typ
     !----------------------------------------------------------------
     ! forces are now in the reciprocal lattice transform it to
     ! cartesian coordinates
     !----------------------------------------------------------------
     CALLMPI( M_sum_d(GRIDC%COMM, EIFOR_SOL(1,1),T_INFO%NIONS*3))
     
     CALL  DIRKAR(T_INFO%NIONS,EIFOR_SOL,LATT_CUR%B)

     DEALLOCATE(GWORK,PCWORK)

#ifdef debugsol
     sol_io 'FORCE CORRECTIONS'
     IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,'(1X,3F12.8)')((EIFOR_SOL(I,J),I=1,3),J=1,T_INFO%NIONS)
#endif

     RETURN

   END SUBROUTINE Fcorrection


!>==========================SUBROUTINE GETRCORE=======================
!!
!!get the core radii from potcar
!!
!!=====================================================================
SUBROUTINE GETRCORE(NTYP,RCS)
 USE constant 
 USE mpimy

 implicit none

 INTEGER IERR, I, ist,L
 INTEGER, INTENT(IN) :: NTYP
 INTEGER,PARAMETER :: ISDIM=100
 CHARACTER (80) STRING(ISDIM)
 CHARACTER (80) STRING1,STRING2,VALUE
 CHARACTER (6)  TAG
 REAL(q) :: RDUM 
 REAL(q),INTENT(OUT) :: RCS(NTYP)

 OPEN(UNIT=10,FILE='POTCAR',STATUS='OLD',IOSTAT=IERR)
 REWIND(10)
 IF (IERR/=0) THEN
    OPEN(UNIT=10,FILE='POTCAR',STATUS='OLD')
 ENDIF
 
 I = 1
 READ(10,'(A)',IOSTAT=ist) STRING1
 DO WHILE(ist == 0)
   READ(10,'(A)',IOSTAT=ist) STRING1
   TAG = STRING1(4:9)
   L = LEN(TAG)
   IF (TAG(1:L)=='RCORE') THEN
      VALUE = STRING1(13:80) 
!      print *, "VALUE",VALUE
      READ(VALUE,'(G10.3)',IOSTAT=IERR) RDUM
!      print *, "RDUM", RDUM
      RCS(I) = RDUM
      I = I + 1
      CYCLE
   ENDIF
 ENDDO

 CLOSE(10) 
 
 RCS = RCS * AUTOA

 RETURN

END SUBROUTINE GETRCORE


!>==========================SUBROUTINE JDFTX_INTERFACE=================
!!
!!interface to the jdftx library
!!
!! returns Vcorr, Vdiel and Atot
!!
!! Vcorr = the electrostatic contrib.
!! Vdiel = epsilon and kappa derivative terms + pot. due to cavitation
!! Atot  = electrostatic contribution to the solv. energy + cavitation 
!!
!!=====================================================================
SUBROUTINE JDFTX_INTERFACE(SP_GRIDC, SP_LATT_CUR, TOTN, PS_TOTN, &
     & Vdiel, Vcorr, Atot)

 USE lattice, ONLY : latt
 USE mgrid, ONLY : grid_3d

 IMPLICIT NONE

 TYPE (latt),INTENT(IN) :: SP_LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: SP_GRIDC

 COMPLEX(q) ::  Vcorr(SP_GRIDC%MPLWV), Vdiel(SP_GRIDC%MPLWV)
 COMPLEX(q), INTENT(IN) :: TOTN(SP_GRIDC%MPLWV)
 COMPLEX(q), INTENT(IN) :: PS_TOTN(SP_GRIDC%MPLWV)
 REAL(q) :: Atot
  
 REAL(q), ALLOCATABLE ::  A_n_r(:),A_rho_r(:)
!  REAL(q), ALLOCATABLE ::  A_n_r_tmp(:),A_rho_r_tmp(:)
 REAL(q), ALLOCATABLE ::  n_r(:),rho_r(:)
 COMPLEX(q), ALLOCATABLE ::  A_n_c(:),A_rho_c(:)

#ifdef jdftx
 ALLOCATE(A_n_r(SP_GRIDC%NPLWV),A_rho_r(SP_GRIDC%NPLWV))
 ALLOCATE(n_r(SP_GRIDC%NPLWV),rho_r(SP_GRIDC%NPLWV))
 ALLOCATE(A_n_c(SP_GRIDC%MPLWV),A_rho_c(SP_GRIDC%MPLWV))

!TOTN to real space
 CALL FFT3D(TOTN,SP_GRIDC,1)

!PS_TOTN to real space
 CALL FFT3D(PS_TOTN,SP_GRIDC,1)

!rho_explicit, real in real,old
!  rho_r = REAL(TOTN)

!merge TOTN to rho_r
 CALL MRG_GRID_RL(SP_GRIDC, rho_r,TOTN)

!merge PS_TOTN to n_r
 CALL MRG_GRID_RL(SP_GRIDC,n_r,PS_TOTN)

!n_cavity, real in real,old
! n_r = REAL(PS_TOTN)

!remove scaling
 rho_r = rho_r/SP_LATT_CUR%OMEGA
 n_r = n_r/SP_LATT_CUR%OMEGA
!  ALLOCATE(A_n_r_tmp(SP_GRIDC%NPLWV),A_rho_r_tmp(SP_GRIDC%NPLWV))

#ifdef MPI
 IF (SP_GRIDC%COMM%NODE_ME==SP_GRIDC%COMM%IONODE) THEN
#endif

    print *, 'Calling minimizefluid'
    call minimizefluid(Atot, n_r, rho_r, A_n_r, A_rho_r)
    print *, "Atot = ", Atot

#ifdef MPI
 ENDIF
#endif

 CALLMPI( M_bcast_d(SP_GRIDC%COMM, Atot, 1)) 

!derivatives wrt free energy 
! A_n_c = CMPLX(A_n_r)
! A_rho_c = CMPLX(A_rho_r)

!ditribute A_n_c and A_rhol_c
 CALL DIS_GRID_RL(SP_GRIDC, A_n_r, A_n_c, .TRUE.)
 CALL DIS_GRID_RL(SP_GRIDC, A_rho_r, A_rho_c, .TRUE.)

!to fourier space
 CALL FFT3D(A_n_c,SP_GRIDC,-1)
 CALL FFT3D(A_rho_c,SP_GRIDC,-1)
!remove the scaling from FFT
 A_n_c = A_n_c / SP_GRIDC%NPLWV
 A_rho_c = A_rho_c / SP_GRIDC%NPLWV

 Vdiel = A_n_c
 Vcorr = A_rho_c

 DEALLOCATE(A_n_r,A_rho_r)
! DEALLOCATE(A_n_r_tmp,A_rho_r_tmp)
 DEALLOCATE(A_n_c,A_rho_c)
 DEALLOCATE(n_r,rho_r)

#endif

 RETURN

END SUBROUTINE JDFTX_INTERFACE


! vaclevel
#ifdef debugsol_pbz
!>==========================SUBROUTINE GET_FERMISHIFT=================
!!
!! compute the fermi shift
!!
!!=====================================================================
SUBROUTINE GET_FERMISHIFT(LATT_CUR, T_INFO, GRIDC, WDES, CVPSP, CHTOT)

  USE poscar, ONLY : type_info 
  USE lattice, ONLY : latt, DIRKAR 
  USE mgrid, ONLY : grid_3d
  USE pseudo, ONLY : potcar
  USE wave, ONLY : wavedes

  IMPLICIT NONE

  TYPE(type_info),INTENT(IN) :: T_INFO
  TYPE (latt),INTENT(IN) :: LATT_CUR
  TYPE (grid_3d) :: GRIDC
  TYPE (wavedes) :: WDES
   
  INTEGER(q) :: I
  COMPLEX(q) :: CVPSP(GRIDC%MPLWV,WDES%NCDIJ)
  COMPLEX(q) :: CHTOT(GRIDC%MPLWV, WDES%NCDIJ)

  RGRID, ALLOCATABLE :: RCHTOT(:), RV(:)
  COMPLEX(q), ALLOCATABLE::  CWORK(:), CWORK_V(:), CVHAR(:), CV(:) 
  REAL(q) :: VACPOT_PSP, VACPOT

  INTEGER :: NODE_ME, IONODE

  ALLOCATE(CWORK(GRIDC%MPLWV), CWORK_V(GRIDC%MPLWV) )
  ALLOCATE(CVHAR(GRIDC%MPLWV), CV(GRIDC%MPLWV))
  ALLOCATE(RCHTOT(DIMREAL(GRIDC%MPLWV)), RV(DIMREAL(GRIDC%MPLWV)))

  CWORK = 0._q
  CWORK_V = 0._q
  CVHAR = 0._q
  CV = 0._q
  RV = 0._q
  RCHTOT = 0._q

  !copy CHTOT to CWORK                                         
  CALL RC_ADD(CHTOT, 1.0_q, CHTOT, 0.0_q, CWORK, GRIDC)
  !hartree pot, CVHAR from CHTOT
  CALL SimplePoisson(GRIDC, LATT_CUR, CWORK, CVHAR, 1)
  !CVHAR --> real
  CALL RC_ADD(CVHAR, 1.0_q, CVHAR, 0.0_q, CV, GRIDC)
  CALL FFT3D(CV, GRIDC, 1)
  CALL UNPACK_C2R(CV, 1.0_q, CV, 0.0_q, RV, GRIDC)
#ifdef debugsolwf
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing PHI_HAR'
  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'PHI_HAR', RV)
#endif

  CV = 0._q
  RV = 0._q
  ! CVPSP --> real
  CALL RC_ADD(CVPSP, 1.0_q, CVPSP, 0.0_q, CV, GRIDC)
  CALL FFT3D(CV, GRIDC, 1)
  CALL UNPACK_C2R(CV, 1.0_q, CV, 0.0_q, RV, GRIDC)
#ifdef debugsolwf
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing PHI_PSP'
  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'PHI_PSP', RV)
#endif

  CV = 0._q
  RV = 0._q
  ! CVPSP + CVHAR --> real
  CALL RC_ADD(CVPSP, 1.0_q, CVHAR, 1.0_q, CV, GRIDC)
  CALL FFT3D(CV, GRIDC, 1)
  CALL UNPACK_C2R(CV, 1.0_q, CV, 0.0_q, RV, GRIDC)
#ifdef debugsolwf
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing PHI_PSP_HAR'
  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'PHI_PSP_HAR', RV)
#endif

  !charge density to real space, CWORK --> RCHTOT
  CALL FFT3D(CWORK, GRIDC, 1)
  CALL UNPACK_C2R(CWORK, 1.0_q, CWORK, 0.0_q, RCHTOT, GRIDC)

  ! vacuum level
  CALL GET_VAC_LEVEL(GRIDC, RCHTOT, RV, VACPOT_PSP)
#ifdef debugsol
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'VACPOT_PSP = ', VACPOT_PSP
#endif
  VACPOT = VACPOT_POIS - VACPOT_PSP
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) ' FERMI_SHIFT = ', VACPOT
!#ifdef debugsolwf
!  !shift by vacpot                                                
!  RV = RV + VACPOT
!  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing PHI_PSP_HAR (shifted by VACPOT)'
!  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'PHI_PSP_HAR_SHIFTED', RV)
!#endif

  DEALLOCATE(RCHTOT, RV)
  DEALLOCATE(CWORK_V, CWORK, CVHAR, CV)

  RETURN 

END SUBROUTINE GET_FERMISHIFT
#endif
! vaclevel

END MODULE POT_K
