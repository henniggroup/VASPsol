!***********************************************************************
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
!***********************************************************************

#include "symbol.inc"

MODULE POT_K

 USE prec
 USE base
 USE mpimy 

 EXTERNAL initjdftx, getionsigma, minimizefluid

 PUBLIC :: Vcorrection, SET_SOL_IO, SOL_INIT, SOL_WRITER, LSOL, EIFOR_SOL
 PRIVATE :: SOL_READER, IONIC_CHARGE, GETRCORE, Fcorrection

 INTEGER, PRIVATE, SAVE                 :: PK_counter = 1
 REAL(q), ALLOCATABLE, SAVE             :: EIFOR_SOL(:,:)
 REAL(q), ALLOCATABLE, PRIVATE, SAVE    :: sigma_rc_k(:)
 REAL(q), PRIVATE, SAVE                 :: sigma_nc_k,sigma_k,nc_k,eb_k,CORE_C,tau,EDIFFSOL
 REAL(q), PRIVATE, SAVE                 :: p0_k, einfty_k, Nsol_k, T_k
 TYPE (in_struct), PRIVATE, SAVE        :: SP_IO
 LOGICAL, SAVE                          :: LSOL = .FALSE.
 LOGICAL, SAVE                          :: LNLSOL = .FALSE.
 LOGICAL, SAVE                          :: LDEFAULTPCM = .FALSE.
 LOGICAL, SAVE                          :: LJDFTX = .FALSE.

 CONTAINS
 
!>==========================SUBROUTINE Vcorrection======================
!!
!! passes n_cav(n_val+pseudo_core) and rho(n_val+N_gauss) to jdftx which
!! returns the derivative of free energy wrt to n_cav and rho,
!! returns Vcorr and Ecorr
!!
!!======================================================================

 SUBROUTINE Vcorrection(SP_GRIDC,SP_P,SP_LATT_CUR,SP_T_INFO,WDES,SP_CHTOT,Vcorr,Ecorr)

  USE constant 
  USE poscar, ONLY : type_info 
  USE lattice, ONLY : latt, DIRKAR 
  USE mgrid, ONLY : grid_3d
  USE pseudo, ONLY : potcar
  USE wave, ONLY : wavedes
!  USE pot_lpcm_k, ONLY : Vcorrection_lpcm, MY_D_PROD
  USE pot_nln_k, ONLY : SET_PARAMS_LPCM, Vcorrection_lpcm, MY_D_PROD
  USE pot_nln_k, ONLY : CONSTRUCTOR_NLN, Vcorrection_nln_pcm, DESTRUCTOR_NLN

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
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Begin minimization'
#endif
!-------------------------------------------------------------
!
!
!              BEGIN Solvation minimization
!
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
                            & sigma_rc_k, tau, ediffsol, SP_T_INFO%NTYP)

     CALL Vcorrection_lpcm(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, TOTN, &
                           & PS_TOTN, Vcorr, Vdiel, Adiel)

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

     CALL CONSTRUCTOR_NLN(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, &
                          & sigma_nc_k, sigma_k, nc_k, eb_k, CORE_C, sigma_rc_k, tau, ediffsol, &
                           & p0_k, einfty_k, Nsol_k, T_k)

!VCORR= just the poisson difference, VDIEL= the rest
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
!        ALLOCATE(A_n_r_tmp(SP_GRIDC%NPLWV),A_rho_r_tmp(SP_GRIDC%NPLWV))
#ifdef MPI
         IF (SP_GRIDC%COMM%NODE_ME==SP_GRIDC%COMM%IONODE) THEN
#endif
!               print *,'integral, ncav', sum(n_r)*(SP_LATT_CUR%OMEGA/SP_GRIDC%NPLWV)
!               print *,'integral, rho', sum(rho_r)*(SP_LATT_CUR%OMEGA/SP_GRIDC%NPLWV)
!               print *, shape(n_r),shape(rho_r),shape(A_n_r),shape(A_rho_r)
               print *, 'Calling minimizefluid'
               call minimizefluid(Adiel, n_r, rho_r, A_n_r, A_rho_r)
               print *, "Adiel=", Adiel
!               do i = 1,SP_GRIDC%NPLWV
!                  A_n_r_tmp(i) = A_n_r(i)
!                  A_rho_r_tmp(i) = A_rho_r(i)
!               end do

#ifdef MPI
         ENDIF
#endif

          CALLMPI( M_bcast_d(SP_GRIDC%COMM, Adiel, 1)) 

!derivatives wrt free energy 
!          A_n_c = CMPLX(A_n_r)
!          A_rho_c = CMPLX(A_rho_r)

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
!        DEALLOCATE(A_n_r_tmp,A_rho_r_tmp)
        DEALLOCATE(A_n_c,A_rho_c)
        DEALLOCATE(n_r,rho_r)
#endif
     ENDIF
!-------------------------------------------------------------
!END JDFTx
!-------------------------------------------------------------     

!-------------------------------------------------------------
!
!
!              END Solvation minimization
!
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
!    CALL IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,5,Vcorr,Vdiel)

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
!     print *,'Adiel, Ecorr3, Ecorr4', Adiel, Ecorr3, Ecorr4

     Ecorr = Adiel - Ecorr3 - Ecorr4

#ifdef debugsol
     IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,'(/a, (es15.6))')"E%Ediel = Ecorr", Ecorr
#endif
!     WRITE(*,'(/a, (es15.6))')"E%Ediel = Ecorr", Ecorr

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


!>=====================SUBROUTINE IONIC_CHARGE==========================
!!
!! evaluates the gaussian ionic charge(flag=1) and pseudo-core charge
!!(flag=3)in k space
!!also calculates force corrections(flag=5) in real space
!!======================================================================

SUBROUTINE IONIC_CHARGE (SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_N,flag,Vcorr,Vdiel)
! USE mpimy
! USE prec
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

 COMPLEX(q),INTENT(OUT) :: SP_N(SP_GRIDC%MPLWV)!ionic charge in k space 
! local variables   
 INTEGER,INTENT(IN) :: flag
! REAL(q),INTENT(IN) :: sigma_nc_k,sigma_rc_k,CORE_C
 INTEGER :: SP_NC,SP_N1,SP_N2,SP_N3,SP_NT,SP_NIS,SP_NI,SP_NG
 REAL(q) :: dx, dy,dz,dv, SP_G1, SP_G2, SP_G3,SP_GX, SP_GY, SP_GZ, SP_GSQU
 COMPLEX(q) :: SP_CEXPF
 COMPLEX(q) :: SP_FOR_gion, SP_FOR1_gion, SP_FOR2_gion, SP_FOR3_gion
 COMPLEX(q) :: SP_FOR_cb, SP_FOR1_cb, SP_FOR2_cb, SP_FOR3_cb
 COMPLEX(q) :: imag_i
 COMPLEX(q), OPTIONAL, INTENT(IN) ::  Vcorr(SP_GRIDC%MPLWV),Vdiel(SP_GRIDC%MPLWV)
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
!-------------------------------------------------------------
!phase factor
!-------------------------------------------------------------
           SP_CEXPF=EXP(-CITPI*(SP_G3+SP_G2+SP_G1))*SP_T_INFO%VCA(SP_NT)

          SP_GX=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(1,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
       &_CUR%B(1,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(1,3))
          SP_GY=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(2,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
       &_CUR%B(2,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(2,3))
         SP_GZ=(SP_GRIDC%LPCTX(SP_N1)*SP_LATT_CUR%B(3,1)+SP_GRIDC%LPCTY(SP_N2)*SP_LATT&
       &_CUR%B(3,2)+SP_GRIDC%LPCTZ(SP_N3)*SP_LATT_CUR%B(3,3))
!-------------------------------------------------------------
!G**2
!-------------------------------------------------------------
          SP_GSQU=(SP_GX*SP_GX+SP_GY*SP_GY+SP_GZ*SP_GZ)*(TPI**2)
!-------------------------------------------------------------
!flag=1, gaussian ionic charge
!-------------------------------------------------------------
             IF ( flag == 1) THEN

               SP_N(SP_NG) = SP_N(SP_NG)-(SP_P(SP_NT)%ZVALF)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_nc_k**2))
!-------------------------------------------------------------
!flag=3, pseudo core cHarge density to repel fluid from core
!-------------------------------------------------------------
             ELSE IF ( flag == 3 ) THEN

               SP_N(SP_NG) = SP_N(SP_NG)+CORE_C*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_rc_k(SP_NT)**2))

             ENDIF
!-------------------------------------------------------------------------
!flag=5(not used), force correction due the gaussian ions and chrge ball
!-------------------------------------------------------------------------
             IF (PRESENT(Vcorr) .AND. PRESENT(Vdiel) ) THEN

               SP_CEXPF=EXP(CITPI*(SP_G3+SP_G2+SP_G1))*SP_T_INFO%VCA(SP_NT)
               SP_FOR_gion = (Vcorr(SP_NG))*imag_i*((SP_P(SP_NT)%ZVALF)*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_nc_k**2)))
               SP_FOR_cb = -(Vdiel(SP_NG))*imag_i*(CORE_C*SP_CEXPF*EXP(-0.5_q*SP_GSQU*(sigma_rc_k(SP_NT)**2)))

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

!-------------------------------------------------------------------------
!force corrections to real space
!-------------------------------------------------------------------------
   IF (PRESENT(Vcorr) .AND. PRESENT(Vdiel) ) THEN
!      CALLMPI( M_sum_d(SP_GRIDC%COMM, EIFOR_SOL(1,1),SP_T_INFO%NIONS*3))
      CALL  DIRKAR(SP_T_INFO%NIONS,EIFOR_SOL,SP_LATT_CUR%B)
#ifdef debugsol
      IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*)'FORCE CORRECTIONS'
      IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,'(1X,3F12.8)')((EIFOR_SOL(I,J),I=1,3),J=1,SP_T_INFO%NIONS)
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
!!======================================================================

 SUBROUTINE SOL_INIT(INFO,LATT_CUR,GRIDC,T_INFO,P)

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
   
   INTEGER :: IDUM, N, IERR, K, I
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
#ifdef debugsol
   IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Setting up the paramters.Unless explicity specified in the INCAR file these values will be used as the defaults'
#endif

!   delta_grd = ((LATT_CUR%OMEGA)/(GRIDC%NPLWV))**(1._q/3._q)
   delta_grd = MAXVAL(LATT_CUR%ANORM/GRIDC%NGPTAR)
#ifdef debugsol      
   IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'delta_grd = ', delta_grd      
#endif

!   sigma_nc_k = 1.5_q * delta_grd
!   sigma_rc_k = 1.1_q * delta_grd

!   sigma_rc_k = 3.0_q * delta_grd
!   sigma_nc_k = 3.0_q * delta_grd !sigma for point charges

   sigma_nc_k = 1.6_q * delta_grd !sigma for point charges

   CALL GETRCORE(T_INFO%NTYP,RCS)

!   RCS = (/(P(K)%PSRMAX,K=1,T_INFO%NTYP)/) !PSDMAX, maximal r for nonlocal contribution (in fact rmax=PSRMAX/NPSNL*(NPSNL-1))

!   sigma_rc_k = MINVAL(RCS)*2.0/4.0_q   ! sigma for the pseudo core
#ifdef debugsol
  IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'RCOREs = ', RCS      
#endif
!STOP
   sigma_rc_k = RCS/2.5_q   ! sigma for the pseudo core

!   print *,RCS,sigma_rc_k
!   stop

   CORE_C = 1.0_q

   LOPEN = .FALSE.
   
   OPEN(UNIT=SP_IO%IU5,FILE='INCAR',STATUS='OLD')      

!the following three paramters are essential for all pcm models
!---------------------------------------------------------------
!sigma_nc_k
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'SIGMA_NC_K','=','#',';','F',IDUM,sigma_nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol      
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''SIGMA_NC_K'',the width of ionic gaussians, not specified, using the default'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'sigma_nc_k = ',sigma_nc_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : sigma_nc_k(width of ionic gaussians) = ',sigma_nc_k 

      ENDIF
#endif
      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''SIGMA_NC_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!sigma_rc_k
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'SIGMA_RC_K','=','#',';','F',IDUM,sigma_rc_k0,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''SIGMA_RC_K'',the width of pseudo-core gaussians, not specified, using the default'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'sigma_rc_k = ',sigma_rc_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : sigma_rc_k(width of pseudo-core gaussians) = ',sigma_rc_k0
         DO I=1,T_INFO%NTYP
            sigma_rc_k(I) = sigma_rc_k0
         ENDDO
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''SIGMA_RC_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!CORE_C
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'CORE_C','=','#',';','F',IDUM,CORE_C,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''CORE_C'', magnitude of pseudo-core charge, not specified, using the default'


         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'CORE_C = ', CORE_C          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : CORE_C (magnitude of pseudo-core charge) = ', CORE_C 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''CORE_C''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
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

 END SUBROUTINE SOL_INIT

!>==========================SUBROUTINE SET_SOL_IO============================
!!
!!gets the IO object and ediff from main, sets LSOL to .TRUE. 
!!and initializes the solvennt parameters(eb_k,nc_k,sigma_k,tau, by callingsol_reader)
!!called from main.F(called only if INFO%LSOL is true)
!!
!!===========================================================================

SUBROUTINE SET_SOL_IO(IO,EDIFF)

  IMPLICIT NONE

  TYPE (in_struct),INTENT(IN) :: IO
  REAL(q) :: EDIFF

  SP_IO = IO
  LSOL = .TRUE. 
  EDIFFSOL = EDIFF/100._q

!default values
  sigma_k = 0.6_q
!   nc_k = 0.00473_q
  nc_k = 0.0025_q
  eb_k = 80.0_q
  tau = 5.25E-4_q

!default values for nonlinear sol
  p0_k = 1.0 ! not used
  einfty_k = 1.33 * 1.33 ! square of water refractive index 
  Nsol_k = 55.5 !mol/L, for water
  T_k = 300 ! temperature i kelvin

#ifdef debugsol
  IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,'(/A/)')'Reading in solvation model parameters from INCAR. If not present, default values will be used'
#endif

  CALL SOL_READER()

END SUBROUTINE SET_SOL_IO

!>==========================SUBROUTINE SOL_READER============================
!!
!!Reads in the linear pcm model parameters(eb_k,nc_k,sigma_k,tau,ediffsol))
!!for the solvent from INCAR
!!
!!===========================================================================

 SUBROUTINE SOL_READER()

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
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'LNLSOL','=','#',';','L',IDUM,RDUM,CDUM,LNLSOL,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''LNLSOL'', switch for nonlinear solvation'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'LNLSOL = ', LNLSOL          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : LNLSOL(Non-linear solvation) = ', LNLSOL 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''LNLSOL''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!P0_K
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'P0_K','=','#',';','F',IDUM,P0_K,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''P0_K'', dipole moment of the molecule'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'P0_K = ', P0_K          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : P0_K = ', P0_K
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''P0_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!EINFTY_K
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'EINFTY_K','=','#',';','F',IDUM,EINFTY_K,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''EINFTY_K'', \epsilon_{\infty}'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'EINFTY_K = ', EINFTY_K         
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : EINFTY_K = ', EINFTY_K
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''EINFTY_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!NSOL_K
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'NSOL_K','=','#',';','F',IDUM,NSOL_K,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''NSOL_K'', Number density of the solvent molecules'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'NSOL_K = ', NSOL_K         
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : NSOL_K = ', NSOL_K
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''NSOL_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!T_K
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'T_K','=','#',';','F',IDUM,T_K,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''T_K'', Temperature in Kelvin'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'T_K = ', T_K         
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : T_K = ', T_K
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''T_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!sigma_k
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'SIGMA_K','=','#',';','F',IDUM,sigma_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''SIGMA_K'',the width of dielectric cavity, not specified, using the default(fit parameter for solvent = water)'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'sigma_k = ',sigma_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : sigma_k(width of dielectric cavity) = ',sigma_k 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''SIGMA_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!nc_k
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'NC_K','=','#',';','F',IDUM,nc_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''NC_K'', cavity turn-on charge density, not specified, using the default(fit parameter for solvent = water)'

         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'nc_k = ', nc_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : nc_k (cavity turn-on charge density) = ', nc_k 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''NC_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!eb_k
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'EB_K','=','#',';','F',IDUM,eb_k,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''EB_K'', relative permittivity of the bulk solvent, not specified, using the default(for solvent = water)'

         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'eb_k = ', eb_k          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : eb_k (relative permittivity of the bulk solvent) = ', eb_k 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''EB_K''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!tau
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'TAU','=','#',';','F',IDUM,tau,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''TAU'',cavity surface tension, not specified, using the default'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'tau = ',tau
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : tau(cavity surface tension) = ',tau
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''TAU''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF

!---------------------------------------------------------------
!EDIFFSOL
!---------------------------------------------------------------
      CALL RDATAB(LOPEN,'INCAR',SP_IO%IU5,'EDIFFSOL','=','#',';','F',IDUM,EDIFFSOL,CDUM,LDUM,CHARAC,N,1,IERR)
#ifdef debugsol
      IF (IERR == 3 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'''EDIFFSOL'', solvation convergence tolerance, not specified, using the default(EDIFF/100)'

         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'EDIFFSOL = ', EDIFFSOL          
         
      ELSE IF (IERR == 0 ) THEN
      
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Found : EDIFFSOL = ', EDIFFSOL 
         
      ENDIF
#endif

      IF ((IERR/=0).AND.(IERR/=3)) THEN
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)'Error reading ''EDIFFSOL''from INCAR.'
         IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) ' Error code = ', IERR
         STOP
      ENDIF


END SUBROUTINE SOL_READER

!>==========================SUBROUTINE SOL_WRITER============================
!!
!!Writes the linear pcm model parameters(eb_k,nc_k,sigma_k,tau))
!!for the solvent to outcar
!!
!!===========================================================================

SUBROUTINE SOL_WRITER()

 IMPLICIT NONE
 
 IF (SP_IO%IU6>=0) THEN
   WRITE(SP_IO%IU6,100) eb_k, sigma_k, nc_k, tau, LNLSOL
 ENDIF

100   FORMAT( &
             ' Solvation parameters'/ &
             '   EB_K    =',F10.6,  '     relative permittivity of the bulk solvent' /&
             '   SIGMA_K =',F10.6,  '     width of the dielectric cavity' /&
             '   NC_K    =',F10.6,  '     cutoff charge density' /&
             '   TAU     =',F10.6,  '     cavity surface tension' /&
             '   LNLSOL     =',L10,  '    Nonlinear solvation' /)

END SUBROUTINE SOL_WRITER 

!>==========================SUBROUTINE Fcorrection============================
!!
!!Force Correction
!!
!!===========================================================================

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
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
!=======================================================================
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
        PCWORK (N) = (CORE_C)*EXP(-0.5_q*(G**2)*(sigma_rc_k(NT)**2))


      ENDDO

      ion: DO NI=NIS,NIADD+NIS-1
!=======================================================================
! initialise the force on the ion to zero
!=======================================================================
         FOR1=0
         FOR2=0
         FOR3=0
!=======================================================================
! CGXDX,Y,Z = I* the changes in the phase factor g.r on moving one
! reciprocal lattice vector in the x,y,z directions, respectively
!=======================================================================
!=======================================================================
! calculate the total force on the ions by summing over reciprocal
! lattice vectors
! first calculate phase factor:
! there are two version for calculating the phase factor
! on vector machines you might try the first version
! (see stufak.F)
!=======================================================================
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
!=======================================================================
! add the contribution to the force from the present reciprocal lattice
! vector  and multiply by i (ie take imaginary part)
!=======================================================================
           FOR = GWORK(NG)*MULFACT AIMAG(CONJG(Vcorr(NG))*CEXPF)
           FOR = FOR + PCWORK(NG)*MULFACT AIMAG(CONJG(Vdiel(NG))*CEXPF)
           FOR1=FOR1-GRIDC%LPCTX_(N1)*FOR
           FOR2=FOR2-GRIDC%LPCTY_(N2)*FOR
           FOR3=FOR3-GRIDC%LPCTZ_(N3)*FOR
         ENDDO
#ifndef vector
         ENDDO
#endif
!=======================================================================
! multiply forces by 2*Pi
!=======================================================================
         EIFOR_SOL(1,NI)=FOR1*TPI
         EIFOR_SOL(2,NI)=FOR2*TPI
         EIFOR_SOL(3,NI)=FOR3*TPI

      ENDDO ion
      NIS=NIS+NIADD
      ENDDO typ

!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
      CALLMPI( M_sum_d(GRIDC%COMM, EIFOR_SOL(1,1),T_INFO%NIONS*3))

      CALL  DIRKAR(T_INFO%NIONS,EIFOR_SOL,LATT_CUR%B)

      DEALLOCATE(GWORK,PCWORK)
#ifdef debugsol
      IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*)'FORCE CORRECTIONS'
      IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,'(1X,3F12.8)')((EIFOR_SOL(I,J),I=1,3),J=1,T_INFO%NIONS)
#endif
      RETURN
      END SUBROUTINE Fcorrection

!>==========================SUBROUTINE GETRCORE============================
!!
!!get the core radii from potcar
!!
!!===========================================================================

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


END MODULE POT_K

