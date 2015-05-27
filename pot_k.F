!*********************************************************************
! Copyright 2013
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
!*********************************************************************

#include "symbol.inc"

MODULE POT_K

 USE prec
 USE constant
 USE base
 USE mpimy 
 USE pot_nln_k, ONLY : SET_PARAMS_LPCM, Vcorrection_lpcm, MY_D_PROD
 USE pot_nln_k, ONLY : CONSTRUCTOR_NLN, Vcorrection_nln_pcm, DESTRUCTOR_NLN
 USE pot_nln_k, ONLY : GET_VAC_LEVEL, VACPOT_POIS

 EXTERNAL initjdftx, getionsigma, minimizefluid

 PUBLIC :: Vcorrection, SET_SOL_IO, SOL_INIT, SOL_WRITER
 PUBLIC :: LSOL, EIFOR_SOL
 PRIVATE :: SOL_READER, IONIC_CHARGE, GETRCORE
 PRIVATE :: Fcorrection, READ_FROM_INCAR

 INTEGER, PRIVATE, SAVE :: PK_counter = 1

 LOGICAL, SAVE :: LSOL = .FALSE.
 REAL(q), ALLOCATABLE, SAVE :: EIFOR_SOL(:,:)
 REAL(q), ALLOCATABLE, PRIVATE, SAVE :: sigma_rc_k(:), CORE_C(:)
 REAL(q), PRIVATE, SAVE :: sigma_nc_k,sigma_k,nc_k,eb_k
 REAL(q), PRIVATE, SAVE :: tau,EDIFFSOL
 REAL(q), PRIVATE, SAVE :: kappa2b_k, lambda_d_k

 LOGICAL, SAVE :: LNLSOL = .FALSE.
 REAL(q), PRIVATE, SAVE :: p0_k, einfty_k, Nsol_k, T_k

 INTEGER, PRIVATE, SAVE :: NION_COUNT_K
 REAL(q), ALLOCATABLE, PRIVATE, SAVE ::NION_K(:), ZION_K(:), RION_K(:)

 TYPE (in_struct), PUBLIC, SAVE :: SP_IO
 LOGICAL, SAVE :: LBCHG = .FALSE.

 LOGICAL, SAVE :: LDEFAULTPCM = .FALSE.
 LOGICAL, SAVE :: LJDFTX = .FALSE.

 CONTAINS
 
!>==========================SUBROUTINE Vcorrection====================
!!
!! passes n_cav(n_val+pseudo_core) and rho(n_val+N_gauss) to jdftx 
!! which returns the derivative of free energy wrt to n_cav and rho,
!! returns Vcorr and Ecorr
!!
!!====================================================================
 SUBROUTINE Vcorrection(SP_GRIDC,SP_P,SP_LATT_CUR,SP_T_INFO,WDES, &
      & SP_CHTOT,Vcorr,Ecorr)

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
  REAL(q) :: Ecorr3,Ecorr4,Ecorr5,eps_corr,Ecorr6,Ehart
  INTEGER :: I,J,K
  REAL(q)  :: dx,dy,dz,dv

  IF(LSOL) THEN

    ALLOCATE(SP_N(SP_GRIDC%MPLWV),TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(Vdiel(SP_GRIDC%MPLWV))
    ALLOCATE(PS_TOTN(SP_GRIDC%MPLWV))
    ALLOCATE(SP_N_ION(SP_GRIDC%MPLWV))
  
    IF (PK_COUNTER == 1) ALLOCATE(EIFOR_SOL(3,SP_T_INFO%NIONS))

    Vdiel = 0._q
    Ecorr3 = 0._q
    Ecorr4 = 0._q
    Adiel = 0._q
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
       CALL Vcorrection_lpcm(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, &
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
            & NION_K, ZION_K, RION_K, NION_COUNT_K, LBCHG, kappa2b_k)

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
    CALL Fcorrection(SP_GRIDC,SP_P,SP_T_INFO,SP_LATT_CUR,Vcorr,Vdiel)
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
 
END SUBROUTINE Vcorrection


!>=====================SUBROUTINE IONIC_CHARGE========================
!!
!! evaluates the gaussian ionic charge(flag=1) and pseudo-core charge
!!(flag=3)in k space
!!also calculates force corrections(flag=5) in real space
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
          !----------------------------------------------------------
          !flag=5(not used), force correction due the gaussian ions 
          !and chrge ball
          !----------------------------------------------------------
          IF (PRESENT(Vcorr) .AND. PRESENT(Vdiel) ) THEN
             
             SP_CEXPF=EXP(CITPI*(SP_G3+SP_G2+SP_G1))*SP_T_INFO%VCA(SP_NT)
             SP_FOR_gion = (Vcorr(SP_NG))*imag_i*((SP_P(SP_NT)%ZVALF)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_nc_k**2)))
             SP_FOR_cb = -(Vdiel(SP_NG))*imag_i*(CORE_C(SP_NT)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_rc_k(SP_NT)**2)))

             SP_FOR1_gion = SP_FOR1_gion + SP_GRIDC%LPCTX_(SP_N1)*SP_FOR_gion
             SP_FOR2_gion = SP_FOR2_gion + SP_GRIDC%LPCTY_(SP_N2)*SP_FOR_gion
             SP_FOR3_gion = SP_FOR3_gion + SP_GRIDC%LPCTZ_(SP_N3)*SP_FOR_gion

             SP_FOR1_cb = SP_FOR1_cb + SP_GRIDC%LPCTX_(SP_N1)*SP_FOR_cb
             SP_FOR2_cb = SP_FOR2_cb + SP_GRIDC%LPCTY_(SP_N2)*SP_FOR_cb
             SP_FOR3_cb = SP_FOR3_cb + SP_GRIDC%LPCTZ_(SP_N3)*SP_FOR_cb
             
          ENDIF

       ENDDO !end of loop over reciproc vectors

       IF (PRESENT(Vcorr) .AND. PRESENT(Vdiel) ) THEN
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR1_gion,1))
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR2_gion,1))
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR3_gion,1))
          
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR1_cb,1))
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR2_cb,1))
          CALLMPI( M_sum_d(SP_GRIDC%COMM, SP_FOR3_cb,1))

          SP_FOR1_gion =  SP_FOR1_gion * TPI
          SP_FOR2_gion =  SP_FOR2_gion * TPI
          SP_FOR3_gion =  SP_FOR3_gion * TPI

          SP_FOR1_cb =  SP_FOR1_cb * TPI
          SP_FOR2_cb =  SP_FOR2_cb * TPI
          SP_FOR3_cb =  SP_FOR3_cb * TPI

          EIFOR_SOL(1,SP_NI) =  SP_FOR1_gion +  SP_FOR1_cb
          EIFOR_SOL(2,SP_NI) =  SP_FOR2_gion +  SP_FOR2_cb
          EIFOR_SOL(3,SP_NI) =  SP_FOR3_gion +  SP_FOR3_cb

       ENDIF

    ENDDO ions

    SP_NIS = SP_NIS + SP_T_INFO%NITYP(SP_NT)

 END DO type

 !-------------------------------------------------------------------
 !force corrections to real space
 !--------------------------------------------------------------------
 IF (PRESENT(Vcorr) .AND. PRESENT(Vdiel) ) THEN
    !CALLMPI(M_sum_d(SP_GRIDC%COMM, EIFOR_SOL(1,1),SP_T_INFO%NIONS*3))
    CALL  DIRKAR(SP_T_INFO%NIONS,EIFOR_SOL,SP_LATT_CUR%B)
#ifdef debugsol
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'FORCE CORRECTIONS'
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,'(1X,3F12.8)')((EIFOR_SOL(I,J),I=1,3),J=1,SP_T_INFO%NIONS)
#endif
    !STOP
 ENDIF

 RETURN

 END SUBROUTINE IONIC_CHARGE


!>==========================SUBROUTINE SOL_INIT======================
!!
!! set the numerical parameters(sigma_nc_k,sigma_rc_k,core_c)
!! and initiaize the grid for jdftx(if LJDFTX flag is set in INCAR)
!! called in potlok
!! 
!!===================================================================
 SUBROUTINE SOL_INIT(INFO,LATT_CUR,GRIDC,T_INFO,P)

   USE constant
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
   
   REAL(q) :: sigma_gauss,RCS(T_INFO%NTYP)
   
   INTEGER :: IDUM, N, IERR, K, I, II, ATOMIC_NUM
   REAL(q) :: RDUM,delta_grd,sigma_rc_k0
   COMPLEX(q) CDUM
   LOGICAL :: LOPEN,LDUM,LIONODE_FLAG,LFLUEX
   CHARACTER (1) CHARAC
   
   LIONODE_FLAG = .FALSE.
   
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
      CALL READ_FROM_INCAR(LOPEN, 'SIGMA_NC_K', 'F', 1, IDUM, &
           & SIGMA_NC_K, CDUM, LDUM, 'the width of ionic gaussians (Angstrom)')

      !---------------------------------------------------------------
      !sigma_rc_k
      !---------------------------------------------------------------
      CALL READ_FROM_INCAR(LOPEN, 'SIGMA_RC_K', 'F', T_INFO%NTYP, IDUM, &
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
      CALL READ_FROM_INCAR(LOPEN, 'CORE_C', 'F', T_INFO%NTYP, IDUM, &
          & RDUM, CDUM, LDUM, 'magnitude of pseudo-core charge (dimensionless)', &
          & CORE_C)
      !CALL READ_FROM_INCAR(LOPEN, 'CORE_C', 'F', 1, IDUM, &
      !     & CORE_C, CDUM, LDUM, 'magnitude of pseudo-core charge (dimensionless)')

      !---------------------------------------------------------------
      !LJDFTX
      !---------------------------------------------------------------
#ifdef jdftx   
      CALL READ_FROM_INCAR(LOPEN, 'LJDFTX', 'L', 1, IDUM, &
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


!>==========================SUBROUTINE SET_SOL_IO====================
!!
!!gets the IO object and ediff from main, sets LSOL to .TRUE. 
!!and initializes the solvent parameters(eb_k,nc_k,sigma_k,tau, by callingsol_reader)
!!called from main.F(called only if INFO%LSOL is true)
!!
!!===================================================================
SUBROUTINE SET_SOL_IO(IO,EDIFF)

  IMPLICIT NONE

  TYPE (in_struct),INTENT(IN) :: IO
  REAL(q) :: EDIFF

  SP_IO = IO
  LSOL = .TRUE. 
  EDIFFSOL = EDIFF/100._q

  !default values
  ! shape function parameters for linear solvation
  ! linearpcm params, 
  !if LNLSOL = .TRUE. then nc_k and tau will be reset
  sigma_k = 0.6_q !!!determines the width of the interface
  nc_k    = 0.0025_q !!!cut-off charge density,old value: 0.00473_q
  tau     = 5.25E-4_q

  !default values for the paramters that describe the solvent environment
  ! water @ 25C
  eb_k     = 78.4_q !80.0_q
  p0_k     = 0.94/0.393430_q ! in Debye, 1 D = 0.393430 ea0
  einfty_k = 1.78_q !1.33 * 1.33 ! square of water refractive index 
  Nsol_k   = 55.35_q !@ 25C, 997.075/18.015 = 55.35!!!55.5 !mol/L, for water
  T_k      = 298.15_q !!!25C, temperature in kelvin
  !debye length, set to infinity
  lambda_d_k = 1E100_q

  !default ionic species parameters
  ! default electrolyte: 1M NaF
  !number of ionic species, used to 
  !allocate space for the other ionic params
  NION_COUNT_K = 2 

  !ionic concentrations in mol/L
  !default is no electrolyte
  ALLOCATE(NION_K(NION_COUNT_K))
  NION_K = 0._q

  !charges on ionic species (dimensionless)
  ALLOCATE(ZION_K(NION_COUNT_K))
  ZION_K(1) = 0._q
  ZION_K(2) = 0._q

  !ionic radii (Angstrom)
  !default set to some impossible value
  ALLOCATE(RION_K(NION_COUNT_K))
  RION_K(1) = 1000._q * AUTOA
  RION_K(2) = 1000._q * AUTOA

#ifdef debugsol
  IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,'(/A/)')'Reading in solvation model parameters from INCAR.&
       & If not present, default values will be used'
#endif

  !read in the solvation parameters
  CALL SOL_READER()

  RETURN

END SUBROUTINE SET_SOL_IO


!>==========================SUBROUTINE SOL_READER=====================
!!
!!Reads in the linear pcm model parameters
!!(eb_k,nc_k,sigma_k,tau,ediffsol))
!!for the solvent from INCAR
!!
!!===================================================================
SUBROUTINE SOL_READER()

  USE constant
  USE setexm
  USE vaspxml
  
  IMPLICIT NONE
  
  INTEGER :: IDUM, N, IERR
  REAL(q) :: RDUM
  COMPLEX(q) CDUM
  LOGICAL LOPEN,LDUM
  CHARACTER (1) CHARAC
  
  LOPEN = .FALSE.
      
  OPEN(UNIT=SP_IO%IU5,FILE='INCAR',STATUS='OLD')      

  !---------------------------------------------------------------
  !LNLSOL
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'LNLSOL', 'L', 1, IDUM, &
       & RDUM, CDUM, LNLSOL, 'switch for nonlinear solvation')

  ! read the nonlinear solvation parameters
  IF (LNLSOL) THEN
     !reset the default solvation params, nc_k and tau
     nc_k    = 0.00675_q 
     tau     = 9.23E-4_q ! for water
     !---------------------------------------------------------------
     !P0_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(LOPEN, 'P0_K', 'F', 1, IDUM, &
          & P0_K, CDUM, LDUM, 'dipole moment of the molecule (Debye)')

     !---------------------------------------------------------------
     !EINFTY_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(LOPEN, 'EINFTY_K', 'F', 1, IDUM, &
          & EINFTY_K, CDUM, LDUM, '\epsilon_{\infty} (dimensionless)')

     !---------------------------------------------------------------
     !NSOL_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(LOPEN, 'NSOL_K', 'F', 1, IDUM, &
          & NSOL_K, CDUM, LDUM, 'Solvent concentration (mol/L)')

     !---------------------------------------------------------------
     !T_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(LOPEN, 'T_K', 'F', 1, IDUM, &
          & T_K, CDUM, LDUM, 'Temperature (in K)')

     !---------------------------------------------------------------
     !NION_COUNT_K
     !---------------------------------------------------------------
     CALL READ_FROM_INCAR(LOPEN, 'NION_COUNT_K', 'I', 1, NION_COUNT_K, &
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
     CALL READ_FROM_INCAR(LOPEN, 'NION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Ionic concentrations (mol/L)', NION_K)

     !---------------------------------------------------------------
     !ZION_K
     !---------------------------------------------------------------
     IF (.NOT. ALLOCATED(ZION_K)) ALLOCATE(ZION_K(NION_COUNT_K))
     CALL READ_FROM_INCAR(LOPEN, 'ZION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Charge on each ionic species (dimensionless)', ZION_K)

     !---------------------------------------------------------------
     !RION_K
     !---------------------------------------------------------------
     IF (.NOT. ALLOCATED(RION_K)) ALLOCATE(RION_K(NION_COUNT_K))
     CALL READ_FROM_INCAR(LOPEN, 'RION_K', 'F', NION_COUNT_K, IDUM, &
          & RDUM, CDUM, LDUM, 'Ionic radius of  each ionic species (Angstroms)', RION_K)
  
  ENDIF !!! end of if (lnlnsol)

  !---------------------------------------------------------------
  !sigma_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'SIGMA_K', 'F', 1, IDUM, &
       & SIGMA_K, CDUM, LDUM, 'the width of dielectric cavity (dimensionless)')

  !---------------------------------------------------------------
  !nc_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'NC_K', 'F', 1, IDUM, &
       & NC_K, CDUM, LDUM, 'cavity turn-on charge density (Angstrom^-3)')

  !---------------------------------------------------------------
  !eb_k
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'EB_K', 'F', 1, IDUM, &
       & EB_K, CDUM, LDUM, 'relative permittivity of the bulk solvent (dimensionless)')

  !---------------------------------------------------------------
  !tau
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'TAU', 'F', 1, IDUM, &
       & TAU, CDUM, LDUM, 'cavity surface tension (eV/Angstrom^2)')

  !---------------------------------------------------------------
  !EDIFFSOL
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'EDIFFSOL', 'F', 1, IDUM, &
       & EDIFFSOL, CDUM, LDUM, 'solvation convergence tolerance')

  !---------------------------------------------------------------
  !LBCHG
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'LBCHG', 'L', 1, IDUM, &
       & RDUM, CDUM, LBCHG, 'write the bound charge file')

  !---------------------------------------------------------------
  !LAMBDA_D_K
  !---------------------------------------------------------------
  CALL READ_FROM_INCAR(LOPEN, 'LAMBDA_D_K', 'F', 1, IDUM, &
       & LAMBDA_D_K, CDUM, LDUM, 'Debye length (Angstroms)')

  !kappa^2 from debye length
  IF (LAMBDA_D_K > 1E30_q) THEN
     KAPPA2B_K = 0._q
  ELSE
     KAPPA2B_K = 1._q / (LAMBDA_D_K ** 2) * AUTOA**2  ! in 1/a.u^2
  ENDIF

  ! override kappa^2 for nonlinear sol
  ! in 1/a.u^2
  IF (LNLSOL) KAPPA2B_K = 4._q * PI * SUM(NION_K * 6.022 * 1E-4_q * (AUTOA**3) * ZION_K * ZION_K) / (BOLKEV * T_k/(2 * RYTOEV))/ EB_K

  RETURN

END SUBROUTINE SOL_READER


!>==========================SUBROUTINE SOL_WRITER=====================
!!
!!Writes the linear pcm model parameters(eb_k,nc_k,sigma_k,tau))
!!for the solvent to outcar
!!
!!====================================================================
SUBROUTINE SOL_WRITER()

 IMPLICIT NONE
 
 IF (SP_IO%IU6>=0) THEN
   WRITE(SP_IO%IU6,100) LNLSOL, sigma_k, nc_k, T_k, Nsol_k, eb_k, einfty_k, tau
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

END SUBROUTINE SOL_WRITER 


!>==========================SUBROUTINE Fcorrection====================
!!
!!Force Correction
!!
!!====================================================================
SUBROUTINE Fcorrection(GRIDC,P,T_INFO,LATT_CUR,Vcorr,Vdiel)

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
SUBROUTINE READ_FROM_INCAR(LOPEN, VAR_NAME, VAR_TYPE, VAR_SIZE, IDUM,&
     & RDUM, CDUM, LDUM, MSG, RDUM_ARRAY)

  USE prec
  USE mgrid
  
  IMPLICIT NONE

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
     CALL RDATAB(LOPEN, 'INCAR', SP_IO%IU5, VAR_NAME, '=', '#', ';', VAR_TYPE, &
          & IDUM, RDUM, CDUM, LDUM, CHARAC, N, VAR_SIZE, IERR)
  ELSE
     CALL RDATAB(LOPEN, 'INCAR', SP_IO%IU5, VAR_NAME, '=', '#', ';', VAR_TYPE, &
          & IDUM, RDUM_ARRAY, CDUM, LDUM, CHARAC, N, VAR_SIZE, IERR)
  ENDIF

#ifdef debugsol      
  IF (IERR == 3 ) THEN      
     IF (SP_IO%IU0>=0) THEN 
        WRITE(SP_IO%IU0,*) VAR_NAME, ' not found in INCAR. Default value will be used'
        WRITE(SP_IO%IU0,*) MSG
        IF (VAR_TYPE == 'I') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', IDUM
        ELSEIF (VAR_TYPE == 'C') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', CDUM
        ELSEIF (VAR_TYPE == 'L') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', LDUM
        ELSEIF ( (VAR_TYPE == 'F') .AND. (VAR_SIZE >= 1) .AND. PRESENT(RDUM_ARRAY) ) THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', RDUM_ARRAY
        ELSEIF (VAR_TYPE == 'F') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', RDUM
        ENDIF
     ENDIF
  ELSE IF (IERR == 0 ) THEN
     IF (SP_IO%IU0>=0) THEN
        WRITE(SP_IO%IU0,*)'Found : ', VAR_NAME, ', ', MSG
        IF (VAR_TYPE == 'I') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', IDUM
        ELSEIF (VAR_TYPE == 'C') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', CDUM
        ELSEIF (VAR_TYPE == 'L') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', LDUM
        ELSEIF ( (VAR_TYPE == 'F') .AND. (VAR_SIZE >= 1) .AND. PRESENT(RDUM_ARRAY) ) THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', RDUM_ARRAY
        ELSEIF (VAR_TYPE == 'F') THEN
           WRITE(SP_IO%IU0,*) VAR_NAME, ' = ', RDUM
        ENDIF
     ENDIF
  ELSEIF ( IERR == 1 .OR. IERR == 2) THEN
     IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) 'Error reading INCAR'
     STOP
  ENDIF
#endif

  IF ((IERR/=0).AND.(IERR/=3)) THEN
     IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) 'Error reading ', VAR_NAME, 'from INCAR.'
     IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
     STOP
  ENDIF

END SUBROUTINE READ_FROM_INCAR


END MODULE POT_K

