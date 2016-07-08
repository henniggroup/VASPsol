!**********************************************************************
! Copyright 2013-2015 Kiran Mathew and Richard Hennig
!
!Licensed under the Apache License, Version 2.0 (the "License");
!you may not use this file except in compliance with the License.
!You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
!Unless required by applicable law or agreed to in writing, software
!distributed under the License is distributed on an "AS IS" BASIS,
!WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!See the License for the specific language governing permissions and
!limitations under the License.
!**********************************************************************

#include "symbol.inc"

MODULE POT_K

 USE prec
 USE constant
 USE base
 USE mpimy 
 USE pot_nln_k, ONLY : SET_PARAMS_LPCM, VCORRECTION_LPCM, MY_D_PROD
 USE pot_nln_k, ONLY : CONSTRUCTOR_NLN, VCORRECTION_NLN_PCM, DESTRUCTOR_NLN
 USE pot_nln_k, ONLY : SIMPLEPOISSON, GET_VAC_LEVEL, VACPOT_POIS

 EXTERNAL initjdftx, getionsigma, minimizefluid

 PUBLIC :: READ_FROM_INCAR, SOL_INIT, VCORRECTION
! vaclevel
 PUBLIC :: GET_FERMISHIFT
! vaclevel

 PRIVATE :: SOL_READER, IONIC_CHARGE, GETRCORE
 PRIVATE :: Fcorrection

 INTEGER, PRIVATE, SAVE :: PK_counter = 1

 LOGICAL, PRIVATE, SAVE :: LSOL = .FALSE.
 REAL(q), ALLOCATABLE, PRIVATE, SAVE :: sigma_rc_k(:), CORE_C(:)
 REAL(q), PRIVATE, SAVE :: sigma_nc_k,sigma_k,nc_k,eb_k
 REAL(q), PRIVATE, SAVE :: tau,EDIFFSOL
 REAL(q), PRIVATE, SAVE :: kappa2b_k, lambda_d_k
 ! nonlinear
 LOGICAL, PRIVATE, SAVE :: LNLSOL = .FALSE.
 REAL(q), PRIVATE, SAVE :: p0_k, einfty_k, Nsol_k, T_k
 INTEGER, PRIVATE, SAVE :: NION_COUNT_K
 REAL(q), ALLOCATABLE, PRIVATE, SAVE ::NION_K(:), ZION_K(:), RION_K(:)
 ! io
 LOGICAL, PRIVATE, SAVE :: LRHOB = .FALSE.
 LOGICAL, PRIVATE, SAVE :: LRHOION = .FALSE.
 TYPE (in_struct), PRIVATE, SAVE :: SP_IO

 LOGICAL, PRIVATE, SAVE :: LDEFAULTPCM = .FALSE.
 LOGICAL, PRIVATE, SAVE :: LJDFTX = .FALSE.

 CONTAINS
 
!>==========================SUBROUTINE VCORRECTION====================
!!
!! passes n_cav(n_val+pseudo_core) and rho(n_val+N_gauss) to jdftx 
!! which returns the derivative of free energy wrt to n_cav and rho,
!! returns Vcorr and Ecorr
!!
!!====================================================================
 SUBROUTINE VCORRECTION(SP_GRIDC, SP_P, SP_LATT_CUR, SP_T_INFO, WDES, &
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

  REAL(q), ALLOCATABLE ::  A_n_r(:),A_rho_r(:)
  !  REAL(q), ALLOCATABLE ::  A_n_r_tmp(:),A_rho_r_tmp(:)
  REAL(q), ALLOCATABLE ::  n_r(:),rho_r(:)
  COMPLEX(q), ALLOCATABLE ::  A_n_c(:),A_rho_c(:)
  REAL(q) :: Adiel

  COMPLEX(q), INTENT(OUT) ::  Vcorr(SP_GRIDC%MPLWV)
  COMPLEX(q) :: SP_CHTOT(SP_GRIDC%MPLWV, WDES%NCDIJ)
  COMPLEX(q),ALLOCATABLE :: SP_N(:),TOTN(:),SP_N_ION(:)
  COMPLEX(q),ALLOCATABLE :: Vdiel(:)
  COMPLEX(q),ALLOCATABLE :: PS_TOTN(:)
  REAL(q),INTENT(OUT) :: Ecorr
  REAL(q) :: EIFOR_SOL(3, SP_T_INFO%NIONS)
  REAL(q) :: Ecorr3,Ecorr4,Ecorr5,eps_corr,Ecorr6,Ehart
  INTEGER :: I,J,K
  REAL(q)  :: dx,dy,dz,dv

  IF(LSOL) THEN

    ALLOCATE(SP_N(SP_GRIDC%MPLWV),TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(Vdiel(SP_GRIDC%MPLWV))
    ALLOCATE(PS_TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(SP_N_ION(SP_GRIDC%MPLWV))
  
    Vdiel = 0._q
    Ecorr3 = 0._q
    Ecorr4 = 0._q
    Adiel = 0._q
    EIFOR_SOL = 0._q

#ifdef debugsol
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Begin Solvation calculations'
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Setting up the charge densities'
#endif
    !-------------------------------------------------------------
    !Gaussian ionic charge in fourier space
    !-------------------------------------------------------------
    CALL IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,1)
    SP_N_ION = SP_N

    !-------------------------------------------------------------
    !total charge in fourier space
    !-------------------------------------------------------------
    CALL RC_ADD(SP_CHTOT(1,1),1.0_q,SP_N_ION(1),1.0_q,TOTN(1), &
         & SP_GRIDC)
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
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Begin minimization'
#endif
    !-------------------------------------------------------------
    !
    !              BEGIN Solvation minimization
    !
    !-------------------------------------------------------------
    
    !-------------------------------------------------------------
    !LINEAR PCM
    !-------------------------------------------------------------
    IF (LDEFAULTPCM .AND. (.NOT. LNLSOL)) THEN

#ifdef debugsol
       IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Calling vcorr_lpcm minimization'
#endif
       CALL SET_PARAMS_LPCM(sigma_nc_k, sigma_k, nc_k, eb_k, CORE_C, &
            & sigma_rc_k, tau, ediffsol, SP_T_INFO%NTYP, kappa2b_k)
       CALL VCORRECTION_LPCM(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, &
            & SP_T_INFO, TOTN, PS_TOTN, Vcorr, Vdiel, Adiel)
#ifdef debugsol
       IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Return from vcorr_lpcm minimization'
#endif

       !-------------------------------------------------------------
       !NON-LINEAR PCM
       !-------------------------------------------------------------
    ELSEIF (LDEFAULTPCM .AND. LNLSOL) THEN

#ifdef debugsol
       IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Calling vcorr_nln_pcm minimization'
#endif
       CALL CONSTRUCTOR_NLN(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, &
            & SP_T_INFO, sigma_nc_k, sigma_k, nc_k, eb_k, CORE_C, &
            & sigma_rc_k, tau, ediffsol, p0_k, einfty_k, Nsol_k, T_k,&
            & NION_K, ZION_K, RION_K, NION_COUNT_K, LRHOB, kappa2b_k)

       !VCORR= delta phi, VDIEL= the rest
       CALL VCORRECTION_NLN_PCM(TOTN, PS_TOTN, VDIEL, VCORR, ADIEL )

       CALL DESTRUCTOR_NLN()

#ifdef debugsol
       IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Return from vcorr_nln_pcm minimization'
#endif

       !-------------------------------------------------------------
       !BEGIN JDFTx
       !-------------------------------------------------------------
#ifdef jdftx
    ELSE
       ALLOCATE(A_n_r(SP_GRIDC%NPLWV),A_rho_r(SP_GRIDC%NPLWV))
       ALLOCATE(n_r(SP_GRIDC%NPLWV),rho_r(SP_GRIDC%NPLWV))
       ALLOCATE(A_n_c(SP_GRIDC%MPLWV),A_rho_c(SP_GRIDC%MPLWV))
       !TOTN to real space
       CALL FFT3D(TOTN,SP_GRIDC,1)
       !PS_TOTN to real space
       CALL FFT3D(PS_TOTN,SP_GRIDC,1)

       !rho_explicit, real in real,old
       !        rho_r = REAL(TOTN)
       
       !merge TOTN to rho_r
       CALL MRG_GRID_RL(SP_GRIDC, rho_r,TOTN)

       !merge PS_TOTN to n_r
       CALL MRG_GRID_RL(SP_GRIDC,n_r,PS_TOTN)
       
       !n_cavity, real in real,old
       !         n_r = REAL(PS_TOTN)

       !remove scaling
       rho_r = rho_r/SP_LATT_CUR%OMEGA
       n_r = n_r/SP_LATT_CUR%OMEGA
       !ALLOCATE(A_n_r_tmp(SP_GRIDC%NPLWV),A_rho_r_tmp(SP_GRIDC%NPLWV))
#ifdef MPI
       IF (SP_GRIDC%COMM%NODE_ME==SP_GRIDC%COMM%IONODE) THEN
#endif
          print *, 'Calling minimizefluid'
          call minimizefluid(Adiel, n_r, rho_r, A_n_r, A_rho_r)
          print *, "Adiel=", Adiel
          !do i = 1,SP_GRIDC%NPLWV
          !   A_n_r_tmp(i) = A_n_r(i)
          !   A_rho_r_tmp(i) = A_rho_r(i)
          !end do
#ifdef MPI
       ENDIF
#endif
       CALLMPI( M_bcast_d(SP_GRIDC%COMM, Adiel, 1)) 

       !derivatives wrt free energy 
       !          A_n_c = CMPLX(A_n_r)
       !          A_rho_c = CMPLX(A_rho_r)
       
       !distribute A_n_c and A_rhol_c
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
    ENDIF
    !-------------------------------------------------------------
    !END JDFTx
    !-------------------------------------------------------------     

    !-------------------------------------------------------------
    !
    !              END Solvation minimization
    !
    !-------------------------------------------------------------    
#ifdef debugsol
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'End minimization'
#endif
    CALL SETUNB(Vcorr,SP_GRIDC)
    CALL SETUNB(Vdiel,SP_GRIDC)

    !Ecorr3 = Vcorr*n(just electrostatic correction)
    CALL MY_D_PROD(Ecorr3, SP_CHTOT(1,1), Vcorr, SP_GRIDC)
    !    CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr3,1))
    CALLMPI( M_sum_s(SP_GRIDC%COMM,1,Ecorr3,0,0,0))

    !Ecorr4 = Vdiel*n
    CALL MY_D_PROD(Ecorr4, SP_CHTOT(1,1), Vdiel, SP_GRIDC)
    !    CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr4,1))
    CALLMPI( M_sum_s(SP_GRIDC%COMM,1,Ecorr4,0,0,0))

    !calculate force correction,flag=5(not used)
    CALL Fcorrection(SP_GRIDC, SP_P, SP_T_INFO, SP_LATT_CUR, &
         & Vcorr, Vdiel, EIFOR_SOL)
    !CALL IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,5,Vcorr,Vdiel)

    !-------------------------------------------------------------
    !Vcorr = electrostatic correction + cavity correction(Vdiel)
    !-------------------------------------------------------------
    Vcorr = Vcorr + Vdiel

    !-------------------------------------------------------------
    !net correction to the total energy, TOTEN
    !-------------------------------------------------------------
#ifdef debugsol
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Adiel, Ecorr3, Ecorr4 ', Adiel, Ecorr3, Ecorr4
#endif

    Ecorr = Adiel - Ecorr3 - Ecorr4

#ifdef debugsol
    IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,'(/a, (es15.6))')"E%Ediel = Ecorr", Ecorr
#endif

    !-------------------------------------------------------------
    !keeps count of the number of times Vcorrection is called
    !-------------------------------------------------------------
    PK_COUNTER = PK_COUNTER + 1
    
    DEALLOCATE(SP_N,TOTN)
    DEALLOCATE(Vdiel)
    DEALLOCATE(PS_TOTN)
    DEALLOCATE(SP_N_ION)

    !if LSOL = .FALSE., sets Vcorr and Ecorr to 0._q
 ELSE
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'LSOL not set. No solvation for you'
    Vcorr = 0._q
    Ecorr = 0._q

 ENDIF
#ifdef debugsol
 IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'End Solvation calculations'
#endif
 RETURN
 
END SUBROUTINE VCORRECTION


!>=====================SUBROUTINE IONIC_CHARGE========================
!!
!! evaluates the gaussian ionic charge(flag=1) and pseudo-core charge
!!(flag=3)in k space
!!
!!====================================================================
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
 INTEGER :: I,J
 
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
          SP_GX=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(1,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
               &_CUR%B(1,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(1,3))
          SP_GY=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(2,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
               &_CUR%B(2,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(2,3))
          SP_GZ=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(3,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
               &_CUR%B(3,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(3,3))
          !----------------------------------------------------------
          !G**2
          !----------------------------------------------------------
          SP_GSQU=(SP_GX*SP_GX+SP_GY*SP_GY+SP_GZ*SP_GZ)*(TPI**2)
          !---------------------------------------------------------
          !flag=1, gaussian ionic charge
          !--------------------------------------------------------
          IF ( flag == 1) THEN
             SP_N(SP_NG) = SP_N(SP_NG)-(SP_P(SP_NT)%ZVALF)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_nc_k**2))
             !--------------------------------------------------------
             !flag=3, pseudo core cHarge density to repel fluid from core
             !--------------------------------------------------------
          ELSE IF ( flag == 3 ) THEN
             SP_N(SP_NG) = SP_N(SP_NG)+CORE_C(SP_NT)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_rc_k(SP_NT)**2))

          ENDIF

       ENDDO !end of loop over reciproc vectors

    ENDDO ions

    SP_NIS = SP_NIS + SP_T_INFO%NITYP(SP_NT)

 END DO type

 RETURN

 END SUBROUTINE IONIC_CHARGE

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
      & LNLSOL_1, p0_k_1, einfty_k_1, Nsol_k_1, T_k_1, &
      & NION_COUNT_K_1, NION_K_1, ZION_K_1, RION_K_1, &
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
   ! nonlinear
   LOGICAL :: LNLSOL_1
   REAL(q) :: p0_k_1, einfty_k_1, Nsol_k_1, T_k_1
   INTEGER :: NION_COUNT_K_1
   REAL(q) :: NION_K_1(NION_COUNT_K_1), ZION_K_1(NION_COUNT_K_1)
   REAL(q) :: RION_K_1(NION_COUNT_K_1) 
   ! io
   LOGICAL :: LRHOB_1, LRHOION_1

   INTEGER :: IDUM, N, IERR, K, I, II, ATOMIC_NUM
   REAL(q) :: RDUM,delta_grd,sigma_rc_k0(T_INFO%NTYP)
   COMPLEX(q) CDUM
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
   ! nonlinear
   LNLSOL = LNLSOL_1
   p0_k = p0_k_1
   einfty_k = einfty_k_1
   Nsol_k = Nsol_k_1
   T_k = T_k_1
   ! io
   LRHOB = LRHOB_1
   LRHOION = LRHOION_1
   EDIFFSOL = EDIFFSOL_1
   NION_COUNT_K = NION_COUNT_K_1
   IF (.NOT. ALLOCATED(NION_K)) ALLOCATE(NION_K(NION_COUNT_K))
   NION_K = NION_K_1
   IF (.NOT. ALLOCATED(ZION_K)) ALLOCATE(ZION_K(NION_COUNT_K))
   ZION_K = ZION_K_1
   IF (.NOT. ALLOCATED(RION_K)) ALLOCATE(RION_K(NION_COUNT_K))
   RION_K = RION_K_1
   
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
      IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Setting up the paramters.Unless explicity specified in the INCAR file these values will be used as the defaults'
#endif
      ! delta_grd = ((LATT_CUR%OMEGA)/(GRIDC%NPLWV))**(1._q/3._q)
      delta_grd = MAXVAL(LATT_CUR%ANORM/GRIDC%NGPTAR)
#ifdef debugsol      
      IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'delta_grd = ', delta_grd 
#endif

      sigma_nc_k = 1.6_q * delta_grd !sigma for point charges

      CALL GETRCORE(T_INFO%NTYP,RCS)

      !RCS = (/(P(K)%PSRMAX,K=1,T_INFO%NTYP)/) 
      !PSDMAX, maximal r for nonlocal contribution 
      !(in fact rmax=PSRMAX/NPSNL*(NPSNL-1))
      ! sigma for the pseudo core
      ! sigma_rc_k = MINVAL(RCS)*2.0/4.0_q   

      sigma_rc_k = RCS/2.5_q   ! sigma for the pseudo core
#ifdef debugsol
      IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'RCOREs = ', RCS      
#endif

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

      LOPEN = .FALSE.
   
      OPEN(UNIT=SP_IO%IU5,FILE='INCAR',STATUS='OLD')      

      !the following three paramters are essential for all pcm models
      !---------------------------------------------------------------
      !sigma_nc_k
      !---------------------------------------------------------------
      CALL READ_FROM_INCAR(SP_IO, LOPEN, 'SIGMA_NC_K', 'F', 1, IDUM, &
           & SIGMA_NC_K, CDUM, LDUM, 'the width of ionic gaussians (Angstrom)')

      !---------------------------------------------------------------
      !sigma_rc_k
      !---------------------------------------------------------------
      CALL READ_FROM_INCAR(SP_IO, LOPEN, 'SIGMA_RC_K', 'F', T_INFO%NTYP, IDUM, &
          & RDUM, CDUM, LDUM, 'the width of pseudo-core gaussians (Angstrom)', &
          & sigma_rc_k)
      !   CALL READ_FROM_INCAR(LOPEN, 'SIGMA_RC_K', 'F', 1, IDUM, &
      ! & SIGMA_RC_K0, CDUM, LDUM, 'the width of pseudo-core gaussians (Angstrom)')
      !   DO I = 1, T_INFO%NTYP
      !      sigma_rc_k(I) = sigma_rc_k0
      !   ENDDO

      !---------------------------------------------------------------
      !CORE_C
      !---------------------------------------------------------------
      CALL READ_FROM_INCAR(SP_IO, LOPEN, 'CORE_C', 'F', T_INFO%NTYP, IDUM, &
          & RDUM, CDUM, LDUM, 'magnitude of pseudo-core charge (dimensionless)', &
          & CORE_C)
      !CALL READ_FROM_INCAR(LOPEN, 'CORE_C', 'F', 1, IDUM, &
      !     & CORE_C, CDUM, LDUM, 'magnitude of pseudo-core charge (dimensionless)')

      !---------------------------------------------------------------
      !LJDFTX
      !---------------------------------------------------------------
#ifdef jdftx   
      CALL READ_FROM_INCAR(SP_IO, LOPEN, 'LJDFTX', 'L', 1, IDUM, &
           & RDUM, CDUM, LJDFTX, 'Use jdftx library')
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
            CALL initjdftx(LATT_CUR%A(1,1),LATT_CUR%A(1,2), &
                 & LATT_CUR%A(1,3),GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ)
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

 END SUBROUTINE SOL_INIT


!>==========================SUBROUTINE Fcorrection====================
!!
!!Force Correction
!!
!!====================================================================
SUBROUTINE Fcorrection(GRIDC, P, T_INFO, LATT_CUR, Vcorr, Vdiel, EIFOR_SOL)

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

 ! work arrays
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

       GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
       GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
       GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

       G=SQRT(GX**2+GY**2+GZ**2)*TPI

       GWORK (N) = -(P(NT)%ZVALF)*EXP(-0.5_q*(G**2)*(sigma_nc_k**2))
       PCWORK (N) = (CORE_C(NT))*EXP(-0.5_q*(G**2)*(sigma_rc_k(NT)**2))

    ENDDO

    ion: DO NI=NIS,NIADD+NIS-1
       !------------------------------------------------------------
       ! initialise the force on the ion to zero
       !------------------------------------------------------------
       FOR1=0
       FOR2=0
       FOR3=0
       !------------------------------------------------------------
       ! CGXDX,Y,Z = I* the changes in the phase factor g.r on 
       ! moving one reciprocal lattice vector in the x,y,z 
       !directions, respectively
       !
       ! calculate the total force on the ions by summing over 
       ! reciprocal lattice vectors
       ! first calculate phase factor:
       ! there are two version for calculating the phase factor
       ! on vector machines you might try the first version
       ! (see stufak.F)
       !------------------------------------------------------------
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
                ! add the contribution to the force from the 
                ! present reciprocal lattice
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
          !---------------------------------------------------
          ! multiply forces by 2*Pi
          !---------------------------------------------------
          EIFOR_SOL(1,NI)=FOR1*TPI
          EIFOR_SOL(2,NI)=FOR2*TPI
          EIFOR_SOL(3,NI)=FOR3*TPI
          
       ENDDO ion

       NIS=NIS+NIADD

    ENDDO typ

    !---------------------------------------------------
    ! forces are now in the reciprocal lattice transform it to
    ! cartesian coordinates
    !---------------------------------------------------
    CALLMPI( M_sum_d(GRIDC%COMM, EIFOR_SOL(1,1),T_INFO%NIONS*3))
    
    CALL  DIRKAR(T_INFO%NIONS,EIFOR_SOL,LATT_CUR%B)

    DEALLOCATE(GWORK,PCWORK)
#ifdef debugsol
    IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*)'FORCE CORRECTIONS'
    IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,'(1X,3F12.8)')((EIFOR_SOL(I,J),I=1,3),J=1,T_INFO%NIONS)
#endif

    RETURN

  END SUBROUTINE Fcorrection


!>==========================SUBROUTINE GETRCORE=======================
!!
!!get the core radii from potcar
!!
!!===================================================================
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

!>========================== SUBROUTINE READ_FROM_INCAR ===========
!!
!!
!! The following comment copied from drdatab.F file in vasp.5.lib
!!
! There may be specified five types (given in variable TYPE) for the  
! data to be extracted:                                               
!                                                                     
!    -- 'S'    read some string       --> result goes to variable STR 
!    -- 'I'    read integer data list --> result goes to array INTRES 
!    -- 'F'    read real    data list --> result goes to array FLTRES 
!    -- 'C'    read complex data list --> result goes to array CMPRES 
!    -- 'L'    read logical data list --> result goes to array LOGRES 
!                                                                     
! for TYPE=I,F,C and L the variable N returns the number of elements  
! found in the data list, if an error occurred (e.g. wrong format) N  
! will be zero (and IERR non-zero). For TYPE=S variable N returns the 
! position of the last non-blank character in variable STR ... .      
!                                                                     
! IERR returns an error code if any problem occured (normal: IERR=0)  
!    -- IERR=1   no free I/O-unit found to open file                  
!    -- IERR=2   OPEN error (e.g. file not found, ...)                
!    -- IERR=3   'keyword' (WHAT) not found on specified file         
!    -- IERR=4   invalid data type (TYPE)                             
!    -- IERR=5   error reading/parsing data list (check format!)      
!    -- IERR=6   cannot open scratch file for conversion of data      
!                                                                     
! RDUM_ARRAY : reading in an array of real values, optional
!!====================================================================
SUBROUTINE READ_FROM_INCAR(IO, LOPEN, VAR_NAME, VAR_TYPE, VAR_SIZE, IDUM,&
     & RDUM, CDUM, LDUM, MSG, RDUM_ARRAY)

  USE prec
  USE mgrid
  
  IMPLICIT NONE

  TYPE (in_struct) :: IO

  INTEGER :: IDUM, N, IERR, VAR_SIZE
  !INTEGER :: IDUM_ARRAY(VAR_SIZE)
  REAL(q) :: RDUM
  REAL(q),OPTIONAL :: RDUM_ARRAY(VAR_SIZE)
  COMPLEX(q) :: CDUM !, CDUM_ARRAY(VAR_SIZE)
  LOGICAL :: LOPEN, LDUM !, LDUM_ARRAY(VAR_SIZE)
  CHARACTER (1) CHARAC
  CHARACTER (1) VAR_TYPE
  CHARACTER (*) VAR_NAME
  CHARACTER (*) MSG

  IF ((VAR_SIZE == 1) .AND. (.NOT. PRESENT(RDUM_ARRAY))) THEN
     CALL RDATAB(LOPEN, 'INCAR', IO%IU5, VAR_NAME, '=', '#', ';', VAR_TYPE, &
          & IDUM, RDUM, CDUM, LDUM, CHARAC, N, VAR_SIZE, IERR)
  ELSE
     CALL RDATAB(LOPEN, 'INCAR', IO%IU5, VAR_NAME, '=', '#', ';', VAR_TYPE, &
          & IDUM, RDUM_ARRAY, CDUM, LDUM, CHARAC, N, VAR_SIZE, IERR)
  ENDIF

#ifdef debugsol      
  IF (IERR == 3 ) THEN      
     IF (IO%IU0>=0) THEN 
        WRITE(IO%IU0,*) VAR_NAME, ' not found in INCAR. Default value will be used'
        WRITE(IO%IU0,*) MSG
        IF (VAR_TYPE == 'I') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', IDUM
        ELSEIF (VAR_TYPE == 'C') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', CDUM
        ELSEIF (VAR_TYPE == 'L') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', LDUM
        ELSEIF ( (VAR_TYPE == 'F') .AND. (VAR_SIZE >= 1) .AND. PRESENT(RDUM_ARRAY) ) THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', RDUM_ARRAY
        ELSEIF (VAR_TYPE == 'F') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', RDUM
        ENDIF
     ENDIF
  ELSE IF (IERR == 0 ) THEN
     IF (IO%IU0>=0) THEN
        WRITE(IO%IU0,*)'Found : ', VAR_NAME, ', ', MSG
        IF (VAR_TYPE == 'I') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', IDUM
        ELSEIF (VAR_TYPE == 'C') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', CDUM
        ELSEIF (VAR_TYPE == 'L') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', LDUM
        ELSEIF ( (VAR_TYPE == 'F') .AND. (VAR_SIZE >= 1) .AND. PRESENT(RDUM_ARRAY) ) THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', RDUM_ARRAY
        ELSEIF (VAR_TYPE == 'F') THEN
           WRITE(IO%IU0,*) VAR_NAME, ' = ', RDUM
        ENDIF
     ENDIF
  ELSEIF ( IERR == 1 .OR. IERR == 2) THEN
     IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Error reading INCAR'
     STOP
  ENDIF
#endif

  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Error reading ', VAR_NAME, 'from INCAR.'
     IF (IO%IU0>=0) WRITE(IO%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF

END SUBROUTINE READ_FROM_INCAR


! vaclevel
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
  ! CVPSP + CVHAR --> real
  CALL RC_ADD(CVPSP, 1.0_q, CVHAR, 1.0_q, CV, GRIDC)
  ! vacuum level
  CALL GET_VAC_LEVEL(GRIDC, LATT_CUR, T_INFO, CWORK, CV,& 
       & VACPOT_PSP)
#ifdef debugsol
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'VACPOT_PSP = ', VACPOT_PSP
#endif
  VACPOT = VACPOT_POIS - VACPOT_PSP
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) ' FERMI_SHIFT = ', VACPOT

  DEALLOCATE(RCHTOT, RV)
  DEALLOCATE(CWORK_V, CWORK, CVHAR, CV)

  RETURN 

END SUBROUTINE GET_FERMISHIFT
! vaclevel

END MODULE POT_K
