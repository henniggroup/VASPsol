!*********************************************************************
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
!*********************************************************************

#include "symbol.inc"

MODULE POT_LPCM_K

 USE prec
 USE constant 
 USE base
 USE mpimy 
 USE poscar, ONLY : type_info 
 USE lattice, ONLY : latt, DIRKAR 
 USE mgrid, ONLY : grid_3d
 USE pseudo, ONLY : potcar

 PUBLIC :: Vcorrection_lpcm, MY_D_PROD, GET_VAC_LEVEL
 PUBLIC :: PACK_R2C, UNPACK_C2R, WRITE_TO_FILE, SimplePoisson
 PRIVATE :: CREATECAVITY, MINIMIZE, Leps, EPS_POT
 PRIVATE :: SHAPE_GRADN, NABLAN, LAPLN
 PRIVATE :: COMPUTE_RHOB, COMPUTE_RHOION

 INTEGER, PRIVATE, SAVE :: PK_counter = 1
 REAL(q), PRIVATE, SAVE :: CHGMIN = 1E-10_q
 COMPLEX(q), PRIVATE, ALLOCATABLE, SAVE :: Sol_phi(:)
 REAL(q), ALLOCATABLE, PRIVATE, SAVE :: sigma_rc_k(:), CORE_C(:)
 REAL(q), PRIVATE, SAVE :: sigma_nc_k, sigma_k, nc_k, eb_k
 REAL(q), PRIVATE, SAVE :: kappa2b_k, tau, EDIFFSOL
 LOGICAL, PRIVATE, SAVE :: LRHOB = .FALSE., LRHOION = .FALSE.
 TYPE (in_struct), PUBLIC, SAVE :: SP_IO
 REAL(q), PUBLIC, SAVE :: VACPOT_POIS = 0._q

 CONTAINS

!>==========================SUBROUTINE Vcorrection_lpcm===============
!!
!!Solves the poisson's eqn, returns Vcorr, Vdiel and Atot
!! Vcorr = phi_tilda - phi, only the electrostatic contrib.
!! Vdiel = epsilon and kappa derivative terms + pot. due to cavitation
!! Atot  = Ael + Acav = electrostatic contribution to the solv.energy 
!!+ cavitation 
!!
!!====================================================================
SUBROUTINE Vcorrection_lpcm(SP_GRIDC, SP_LATT_CUR, SP_P, &
     & SP_T_INFO,SP_CHTOT, TOTN, PS_TOTN, Vcorr, Vdiel, Atot, &
     & sigma_nc_k1, sigma_k1,nc_k1, eb_k1, kappa2b_k1, CORE_C1, &
     & sigma_rc_k1, tau1, ediffs, LRHOB1, LRHOION1 )

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: SP_T_INFO
 TYPE(potcar),INTENT(IN) :: SP_P (SP_T_INFO%NTYP)
 TYPE (latt),INTENT(IN) :: SP_LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: SP_GRIDC
 !TYPE (in_struct),INTENT(IN) :: IO
 LOGICAL :: LRHOB1, LRHOION1

 COMPLEX(q), INTENT(OUT) ::  Vcorr(SP_GRIDC%MPLWV)
 COMPLEX(q),ALLOCATABLE :: phi_tilda(:),phi(:),Vcav(:)
 COMPLEX(q),INTENT(OUT) :: Vdiel(SP_GRIDC%MPLWV)
 COMPLEX(q),INTENT(IN) :: SP_CHTOT(SP_GRIDC%MPLWV)
 COMPLEX(q),INTENT(IN) :: TOTN(SP_GRIDC%MPLWV),PS_TOTN(SP_GRIDC%MPLWV) 
 REAL(q),ALLOCATABLE :: SP_d_eps(:), SP_kappa2(:)
 RGRID,ALLOCATABLE :: RWORK(:)
 INTEGER :: I,J, NCGSOL
 REAL(q)  :: Atot, Ael, Acav, Ecorr1, Ecorr2, Ecorr_cav
 REAL(q) :: sigma_nc_k1,sigma_k1,nc_k1,eb_k1,kappa2b_k1
 REAL(q) :: tau1,ediffs
 REAL(q), INTENT(IN) :: sigma_rc_k1(SP_T_INFO%NTYP), CORE_C1(SP_T_INFO%NTYP)
 !vacuum level
#ifdef debugsol_pbz
 RGRID,ALLOCATABLE :: RWORK_PHI(:), RWORK_CHG(:)
 COMPLEX(q),ALLOCATABLE :: CWORK_K(:)
#endif
 !vacuum level
 REAL(q),ALLOCATABLE :: RHOB(:), RHOION(:)

 ALLOCATE(SP_d_eps(DIMREAL(SP_GRIDC%MPLWV)) )
 ALLOCATE(SP_kappa2(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(phi_tilda(SP_GRIDC%MPLWV), phi(SP_GRIDC%MPLWV))
 ALLOCATE(Vcav(SP_GRIDC%MPLWV), RWORK(DIMREAL(SP_GRIDC%MPLWV)))
 !vacuum level
#ifdef debugsol_pbz
 ALLOCATE(RWORK_PHI(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(RWORK_CHG(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(CWORK_K(SP_GRIDC%MPLWV))
#endif
 !vacuum level
 ALLOCATE(RHOB(DIMREAL(SP_GRIDC%MPLWV)) )
 ALLOCATE(RHOION(DIMREAL(SP_GRIDC%MPLWV)) )
      
 IF (PK_COUNTER == 1)  THEN
   ALLOCATE(Sol_phi(SP_GRIDC%MPLWV))
   Sol_phi = 0._q
 ENDIF

!---------------------------------------------------------------
!set values to the parameters
!---------------------------------------------------------------
 sigma_nc_k = sigma_nc_k1 
 sigma_k =  sigma_k1
 nc_k =  nc_k1
 eb_k =  eb_k1
 kappa2b_k =  kappa2b_k1
 CORE_C =  CORE_C1
 sigma_rc_k =  sigma_rc_k1
 tau = tau1
 EDIFFSOL = ediffs
 LRHOB =  LRHOB1
 LRHOION = LRHOION1

!---------------------------------------------------------------
! initialize
!---------------------------------------------------------------
 Ecorr1 = 0.0_q
 Ecorr2 = 0.0_q
 phi_tilda = 0._q
 phi = 0._q
 Vdiel = 0._q
 Vcorr = 0._q
 Vcav = 0._q
! SP_IO = IO

!-------------------------------------------------------------
!create the positon-dependent dielectric constant in real space
!ps_totn is unscaled and in real space after this call
!returns Vcav(cavitaion pot in rec space) and Acav
!-------------------------------------------------------------
 CALL CREATECAVITY(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, PS_TOTN,&
      & SP_d_eps, SP_kappa2, Vcav, Acav)

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, 
!returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
 CALL MINIMIZE(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, SP_d_eps,&
      & SP_kappa2, TOTN, phi_tilda, NCGSOL)

!-------------------------------------------------------------
!ps_totn is unscaled and in real space
!unpack ps_totn from complex array to real array
!-------------------------------------------------------------
 CALL UNPACK_C2R(PS_TOTN, 1.0_q, PS_TOTN, 0.0_q, RWORK, SP_GRIDC)

!-------------------------------------------------------------
!compute the interface contribution to the local pot, Vdiel
!basically the terms involving the derivatives wrt the charge density
!Vdiel includes the terms containing the epsilon derivative and
!the kappa derivative
!ps_totn is unscaled and in real space
!-------------------------------------------------------------
 CALL EPS_POT(SP_GRIDC, SP_LATT_CUR, SP_T_INFO, SP_P,RWORK,phi_tilda,&
      & SP_d_eps, SP_kappa2, Vdiel)

!-------------------------------------------------------------
! multiply by 0.5*e/permittivity of free space
!-------------------------------------------------------------
 Vdiel = 0.5 * Vdiel / EDEPS

!-------------------------------------------------------------
!compute phi, pot due to (n+N) with epsilon=1
!returns properly scaled phi
!-------------------------------------------------------------
 CALL SimplePoisson(SP_GRIDC, SP_LATT_CUR, TOTN, phi, 1)

!---------------------------------------------------------------------
!compute V_correction,Vcorr
!the hartree potential difference between the vacuum case and 
! the solvent case
!---------------------------------------------------------------------
 Vcorr = phi_tilda - phi

!---------------------------------------------------------------------
!add the cavitation correction
!! Vcorr = Vcorr + Vcav
!both Vdiel and Vcav are the derivatives of the free enrgy wrt 
!n_cavity
!where n_cavity=n+pseudo core
!the reason for seperating Vcorr(coulomb) and Vdiel is becasue for the
!force corrections due Vcorr we need to multiply by the derivative of
!total charge density with respect to the positions wheras for 
!Vdiel we need to multiply by the derivativer of the n_cavity wrt 
!the positions
!---------------------------------------------------------------------
 Vdiel = Vdiel + Vcav

!-------------------------------------------------------------
!compute the components of E_correction, Ecorr
!all energies are computed with chargedensities scaled by the 
!lattice volume inorder to be consistent with vasp
!-------------------------------------------------------------
!  Ecorr1 = 0.5*phi_tilda*(n+N), where phi_tilda is the CG soln
!-------------------------------------------------------------
 CALL MY_D_PROD(Ecorr1, TOTN, phi_tilda, SP_GRIDC)
 Ecorr1 = Ecorr1 * 0.5_q
 CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr1,1))

!---------------------------------------------------------------
!Ecorr2 = 0.5*phi*(n+N), where phi is the soln with epsilon = 1
!---------------------------------------------------------------
 CALL MY_D_PROD(Ecorr2, TOTN, phi, SP_GRIDC)
 Ecorr2 = Ecorr2 * 0.5_q
 CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr2,1))

!---------------------------------------------------------------
!Electrostatic contribution to the solvation energy
!---------------------------------------------------------------
 Ael = Ecorr1 - Ecorr2

!---------------------------------------------------------------
!add the cavitaion energy
!Atot = total solavtion energy
!---------------------------------------------------------------
 Atot = Ael + Acav
 
 333 FORMAT('SOL:',' ',I3,'  ',E12.5,'   ',E12.5,'   ',E12.5,'   ',I3)
!---------------------------------------------------------------
! write to OUTCAR and stdout
!---------------------------------------------------------------
 IF (SP_IO%IU0>=0) THEN
   WRITE(17,333)  PK_COUNTER, Ael, Acav, Atot, NCGSOL
   WRITE(SP_IO%IU0,333)  PK_COUNTER, Ael, Acav, Atot, NCGSOL
 ENDIF

 !rhob = div((epsilon-1)grad phi)
 !phi in rec space
 IF (LRHOB) CALL COMPUTE_RHOB(SP_GRIDC, SP_LATT_CUR, SP_T_INFO, & 
      & SP_D_EPS, PHI_TILDA, RHOB)

 !rhion = -kappa^2 * phi, in real space
 !phi in rec space
 !mind: after this the phi will be in real space
 IF (LRHOION) CALL COMPUTE_RHOION(SP_GRIDC, SP_LATT_CUR, SP_T_INFO, &
      & SP_KAPPA2, PHI_TILDA, RHOION)
 
 !vacuum level
#ifdef debugsol_pbz

 !phi to real space
 CALL FFT3D(phi, SP_GRIDC, 1) 
 CALL UNPACK_C2R(phi, 1.0_q, phi, 0.0_q, RWORK_PHI, SP_GRIDC)
#ifdef debugsolwf
 IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) 'Writing PHI_POIS'
 CALL WRITE_TO_FILE(SP_GRIDC, SP_LATT_CUR, SP_T_INFO, 'PHI_POIS', &
      & RWORK_PHI)
#endif
 !chtot to real space
 CALL RC_ADD(SP_CHTOT, 1.0_q, SP_CHTOT, 0.0_q, CWORK_K, SP_GRIDC)
 CALL FFT3D(CWORK_K, SP_GRIDC, 1) 
 CALL UNPACK_C2R(CWORK_K, 1.0_q, CWORK_K, 0.0_q, RWORK_CHG, SP_GRIDC)
 CALL GET_VAC_LEVEL(SP_GRIDC, RWORK_CHG, RWORK_PHI, VACPOT_POIS)
#ifdef debugsol
 IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) 'VACPOT_POIS = ', VACPOT_POIS
#endif

#endif
 !vacuum level

 PK_COUNTER = PK_COUNTER + 1
      
 DEALLOCATE(SP_d_eps, SP_kappa2, phi_tilda, phi, Vcav, RWORK)
 !vacuum level
#ifdef debugsol_pbz
 DEALLOCATE(RWORK_PHI)
 DEALLOCATE(RWORK_CHG)
 DEALLOCATE(CWORK_K)
#endif
 !vacuum level
 DEALLOCATE(RHOB, RHOION)

 RETURN

END SUBROUTINE Vcorrection_lpcm


!>=======================SUBROUTINE CREATECAVITY======================
!!
!!initializes epsilon in the grid
!!calculates the cavitation energy and cavitation potential
!!Acav and Vcav
!!
!!====================================================================
SUBROUTINE CREATECAVITY(CE_GRIDC, CE_LATT_CUR, CE_P, CE_T_INFO, &
     & CE_PS_TOTN,CE_d_eps,CE_kappa2,Vcav,Acav)  

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: CE_T_INFO
 TYPE(potcar),INTENT(IN)    :: CE_P (CE_T_INFO%NTYP)
 TYPE(grid_3d),INTENT(IN)   :: CE_GRIDC
 TYPE(latt),INTENT(IN)      :: CE_LATT_CUR

 COMPLEX(q)  :: CE_PS_TOTN(CE_GRIDC%MPLWV), Vcav(CE_GRIDC%MPLWV)
 REAL(q),INTENT(OUT) :: CE_d_eps(DIMREAL(CE_GRIDC%MPLWV))
 REAL(q),INTENT(OUT) :: CE_kappa2(DIMREAL(CE_GRIDC%MPLWV))

 REAL(q),ALLOCATABLE :: CE_S(:), VWORK(:) ! shape func & Vcav
 !\nabla n and |\naba n|^2
 REAL(q),ALLOCATABLE :: CE_gn(:,:), CE_gn2_tmp(:) 
 COMPLEX(q),ALLOCATABLE :: CE_gn2(:) ! |\naba n|^2
 REAL(q),ALLOCATABLE :: CE_ggn(:,:) !\nabla |\nabla n|
 REAL(q),ALLOCATABLE :: CE_ggn2(:) !| \nabla |\nabla n| |^2
 REAL(q),ALLOCATABLE :: CE_divn(:) ! laplacian of n
 REAL(q),ALLOCATABLE :: CE_term1(:) ! grad shape func wrt n
 ! REAL(q), ALLOCATABLE       :: RWORK(:)
 RGRID, ALLOCATABLE :: RWORK(:)
 !compimentary error function implemented in errf.F in vasp.5.lib
 REAL(q), EXTERNAL :: ERRFC 
 INTEGER(q)             :: I, J, NPC, NPR, RLNP, RCNP
 REAL(q)                :: Acav,CE_qs,dv,tmp,tmp1,tmp2
 !! test 7/8/14
 INTEGER(q)             :: NG, NC, N1, N2, N3, NX, NY, NZ
 REAL(q)                :: R(3)
 !! test 7/8/14

 NPC =  CE_GRIDC%MPLWV
 RLNP = CE_GRIDC%RL%NP
 RCNP = CE_GRIDC%RC%NP
 NPR =  DIMREAL(CE_GRIDC%MPLWV) !NPR = NPC, in complex mode

#ifdef debugsol
 IF (SP_IO%IU0>=0) print *, 'NPLWV, MPLWV, RL%NP, RC%NP',&
      & CE_GRIDC%NPLWV, CE_GRIDC%MPLWV, CE_GRIDC%RL%NP, CE_GRIDC%RC%NP
#endif

 ALLOCATE(CE_gn2(NPC))
 ALLOCATE(CE_S(RLNP))
 ALLOCATE(CE_gn2_tmp(NPR), CE_gn(3, RLNP) )
 ALLOCATE(CE_ggn2(NPR), CE_ggn(3, RLNP) )
 ALLOCATE(CE_divn(NPR), CE_term1(RLNP))
 ALLOCATE(RWORK(NPR), VWORK(NPR))

 !-------------------------------------------------------------
 !volume element
 !-------------------------------------------------------------
 dv = CE_LATT_CUR%OMEGA / CE_GRIDC%NPLWV

 !-------------------------------------------------------------
 !remove the scaling in the reciprocal charge density
 !-------------------------------------------------------------
 CE_PS_TOTN = CE_PS_TOTN / CE_LATT_CUR%OMEGA

 !-------------------------------------------------------------
 !charge density in rec space, returns \nabla n and |\nabla n|^2 
 !in real space
 !-------------------------------------------------------------
 CALL NABLAN(CE_PS_TOTN, CE_GRIDC, CE_LATT_CUR, CE_gn2_tmp, CE_gn)

 !-------------------------------------------------------------
 !charge density in rec space, returns laplacian of n, Lap(n)
 !Lap(n) = (\nabla \dot \nabla n) in real space
 !-------------------------------------------------------------
 CALL LAPLN(CE_PS_TOTN, CE_GRIDC, CE_LATT_CUR, CE_divn)

 !-------------------------------------------------------------
 !cavitaion pot in real space
 !-------------------------------------------------------------
 Vcav = 0._q
 VWORK = 0._q

 !-------------------------------------------------------------
 !add -Lap(n)/|\nabla n| to Vcav and copy \sqrt(|\nabla n|^2) 
 !to CE_gn2_tmp
 !-------------------------------------------------------------
 DO I=1,DIMREAL(CE_GRIDC%MPLWV)
   tmp = SQRT(MAX(REAL(CE_gn2_tmp(I),KIND=q),CHGMIN))
   VWORK(I) = VWORK(I) - (CE_divn(I))/tmp
   CE_gn2_tmp(I) = tmp
 ENDDO

 !-------------------------------------------------------------
 !CE_PS_TOTN to the real space
 !-------------------------------------------------------------
 CALL FFT3D(CE_PS_TOTN,CE_GRIDC,1)

 !-------------------------------------------------------------
 !unpack complex array to a real array
 !in complex mode, just copies CE_PS_TOTN to RWORK
 !-------------------------------------------------------------
 CALL UNPACK_C2R(CE_PS_TOTN, 1.0_q, CE_PS_TOTN, 0.0_q, RWORK,CE_GRIDC)

 !-------------------------------------------------------------
 !charge density must be in real space(unscaled), 
 !returns derivative of shape function wrt the charge density 
 !in the real space,
 ! CE_term1 = gamma*A / n, where
 ! gamma * A = exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )
 !-------------------------------------------------------------
 CALL SHAPE_GRADN(CE_GRIDC, CE_T_INFO, CE_P, CE_LATT_CUR, RWORK, &
      & CE_term1)

 !-------------------------------------------------------------
 ! CE_S = shapefunction = N_lq/N_b, density of the solvent as a 
 ! function of n(r) over the bulk solvent density.
 ! shapefunction value varies from 0 in the solute to 1
 ! in the solvent
 !
 ! CE_d_eps = epsilon  = 1.0 + (eb_k - 1) * shape function
 !
 ! CE_kappa2 = kappa2  = kappa2b_k * shape function
 ! kappa2b_k: inverse of the square of debye length in angstroms
 !-------------------------------------------------------------
 DO I=1, CE_GRIDC%RL%NP

   CE_S(I) = ERRFC( log(MAX(REAL(RWORK(I),KIND=q),CHGMIN)/nc_k) / SQRT(2._q)/sigma_k ) / 2._q

   CE_d_eps(I) = 1._q + (eb_k - 1._q) * CE_S(I)

   CE_kappa2(I) = kappa2b_k * CE_S(I)

 END DO

 !-------------------------------------------------------------
 !quantum surface area
 !-------------------------------------------------------------
 CE_qs = 0._q

 !-------------------------------------------------------------
 !quantum surface area, integral of (gamma*A / n) * |\nabla n| 
 != CE_term1 * \sqrt(CE_gn2_tmp)
 !-------------------------------------------------------------
 DO I=1, CE_GRIDC%RL%NP
   CE_qs = CE_qs + (CE_term1(I)) * (CE_gn2_tmp(I))
 ENDDO

 CE_qs = CE_qs * dv

 CALLMPI( M_sum_d(CE_GRIDC%COMM,CE_qs,1))

 !-------------------------------------------------------------
 !cavitation energy
 !-------------------------------------------------------------
 Acav = tau * CE_qs

#ifdef debugsol
 IF (SP_IO%IU0>=0) THEN
   print *,'tau ',tau
   print *,'surface area ',CE_qs
   print *,'Acav ',Acav
 ENDIF
#endif

 !-------------------------------------------------------------
 !1/ |nabla n|
 !-------------------------------------------------------------
 DO I=1, DIMREAL(CE_GRIDC%MPLWV)
   CE_gn2_tmp(I) = 1._q / ( MAX(CE_gn2_tmp(I), CHGMIN) ) 
 ENDDO

 !-------------------------------------------------------------
 !Compute cavitation potential, Vcav
 !-------------------------------------------------------------
 !packs the real array into a complex one
 !in complex mode, same as  CE_gn2 = CMPLX(CE_gn2_tmp)
 CALL PACK_R2C(CE_gn2_tmp,CE_gn2,CE_GRIDC)

 !-------------------------------------------------------------
 !1/ |nabla n| to reciprocal space to take the gradient
 !-------------------------------------------------------------
 CALL FFT3D(CE_gn2,CE_GRIDC,-1)
 CALL SETUNB(CE_gn2,CE_GRIDC)
 CE_gn2 = CE_gn2/CE_GRIDC%NPLWV

 !-------------------------------------------------------------
 !\nabla(1 / |\nabla n|), CE_ggn in real space
 !-------------------------------------------------------------
 CALL NABLAN(CE_gn2, CE_GRIDC, CE_LATT_CUR, CE_ggn2, CE_ggn)

 !-------------------------------------------------------------
 !add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
 !and multiply by CE_term1 = gamma*A/n in real space
 !-------------------------------------------------------------
 DO I=1, CE_GRIDC%RL%NP
  tmp = ( CE_gn(1,I)*CE_ggn(1,I) + CE_gn(2,I)*CE_ggn(2,I) +&
        & CE_gn(3,I)*CE_ggn(3,I) ) * (CE_term1(I))
  VWORK(I) = VWORK(I) - tmp
 ENDDO

 CALL PACK_R2C(VWORK, Vcav, CE_GRIDC)

 !-------------------------------------------------------------
 !Vcav to reciprocal space
 !-------------------------------------------------------------
 CALL FFT3D(Vcav,CE_GRIDC,-1)
 CALL SETUNB(Vcav,CE_GRIDC)
 Vcav = Vcav/CE_GRIDC%NPLWV

 !-------------------------------------------------------------
 !multiply by the surface tension
 !-------------------------------------------------------------
 Vcav = Vcav * tau

 DEALLOCATE(CE_S, RWORK, VWORK)
 DEALLOCATE(CE_gn2, CE_gn2_tmp)
 DEALLOCATE(CE_gn, CE_ggn2)
 DEALLOCATE(CE_ggn, CE_divn, CE_term1)

 RETURN

END SUBROUTINE CREATECAVITY


!>=====================SUBROUTINE MINIMIZE============================
!!
!!routine for solving poissons eqn using conjugate gradient or 
!!steepest descent
!!
!!====================================================================
SUBROUTINE MINIMIZE(MN_GRIDC, MN_LATT_CUR, MN_P, MN_T_INFO, MN_d_eps,&
     & MN_kappa2, MN_TOTN, MN_phi, NCGSOL)

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: MN_GRIDC
 TYPE (latt),INTENT(IN) :: MN_LATT_CUR
 TYPE (type_info),INTENT(IN) :: MN_T_INFO
 TYPE (potcar),INTENT(IN) :: MN_P (MN_T_INFO%NTYP)

 COMPLEX(q),INTENT(IN) :: MN_TOTN(MN_GRIDC%MPLWV)
 COMPLEX(q),INTENT(OUT) :: MN_phi(MN_GRIDC%MPLWV)
!local variables
 COMPLEX(q),ALLOCATABLE :: MN_RESID(:),MN_Z(:),MN_d(:)
 COMPLEX(q),ALLOCATABLE :: MN_Lp(:)
!parameters in CG
 REAL(q) :: MN_beta,MN_alpha,GX,GY,GZ,rinvLr,eps_bar,r2
 REAL(q),INTENT(IN) ::  MN_d_eps(DIMREAL(MN_GRIDC%MPLWV))
 REAL(q),INTENT(IN) ::  MN_kappa2(DIMREAL(MN_GRIDC%MPLWV))
 COMPLEX(q),ALLOCATABLE :: MN_GSQU(:)

 COMPLEX(q),ALLOCATABLE :: gradphi_x_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_y_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_z_in(:)

 COMPLEX(q),ALLOCATABLE :: phi_work(:)

 INTEGER :: I,NC,N1,N2,N3
 INTEGER :: flag,count
 INTEGER,INTENT(OUT) :: NCGSOL
 REAL(q) :: sum_c,sum_n,sum_t,step,phi_diff,r20,maximum,maximum1,rand
 REAL(q) :: max_lp_r,max_lp_i,max_n_r,max_n_i

 ! ALLOCATING the local arrays
 ALLOCATE(MN_RESID(MN_GRIDC%MPLWV),MN_Z(MN_GRIDC%MPLWV))
 ! ALLOCATE(MN_PH(MN_GRIDC%MPLWV))
 ALLOCATE(MN_Lp(MN_GRIDC%MPLWV),MN_GSQU(MN_GRIDC%MPLWV))
 ALLOCATE(MN_d(MN_GRIDC%MPLWV))

 ALLOCATE(gradphi_x_in(MN_GRIDC%MPLWV),gradphi_y_in(MN_GRIDC%MPLWV))
 ALLOCATE(gradphi_z_in(MN_GRIDC%MPLWV), phi_work(MN_GRIDC%MPLWV))

 count=0
 MN_beta=0._q
 MN_alpha=0._q
! MN_PH=0._q
 MN_phi=0._q
 rinvLr=0._q
 r2 = 0._q
 
 I = 0

 !-------------------------------------------------------------
 !calculating 1/G^2 for preconditioning
 !-------------------------------------------------------------
 col: DO NC=1,MN_GRIDC%RC%NCOL
    N2= MN_GRIDC%RC%I2(NC)
    N3= MN_GRIDC%RC%I3(NC)
    row: DO N1=1,MN_GRIDC%RC%NROW
       I = I +1
       !GX,GY,GZ
       GX = (MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(1,1)+MN_GRIDC%LPCTY(N2)*&
            & MN_LATT_CUR%B(1,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(1,3))

       GY = (MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(2,1)+MN_GRIDC%LPCTY(N2)*&
            & MN_LATT_CUR%B(2,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(2,3))

       GZ = (MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(3,1)+MN_GRIDC%LPCTY(N2)*&
            & MN_LATT_CUR%B(3,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(3,3))

       !set G=0 component to 1/kappa2b_k if kappa2b_k is not zero
       IF ((MN_GRIDC%LPCTX(N1)==0).AND.(MN_GRIDC%LPCTY(N2)==0).AND.&
            &(MN_GRIDC%LPCTZ(N3)==0)) THEN
          
          IF (kappa2b_k<1E-10) THEN
             MN_GSQU(I)= (0.0_q,0.0_q)
          ELSE
             MN_GSQU(I)= 1._q/kappa2b_k 
          ENDIF

       ELSE
          ! 1/(G^2 + kappa2b_k)
          MN_GSQU(I) = 1._q/((GX*GX+GY*GY+GZ*GZ) * (TPI**2) + kappa2b_k)
       ENDIF
       
    ENDDO row
 ENDDO col

 !*********************************************************
 ! start of Pre-conditioned Conjugate Gradient
 !**********************************************************

 !-------------------------------------------------------------
 !intial guess for phi
 !phi0 = (totn*EDEPS/omega) / (G^2 + kappa2b_k)
 !unscaled charge density used
 !-------------------------------------------------------------
 IF (PK_counter ==1) THEN
    DO I=1,MN_GRIDC%RC%NP
       MN_phi(I) = REAL(MN_GSQU(I)) * (MN_TOTN(I)) !* EDEPS / MN_LATT_CUR%OMEGA)
    ENDDO
 ELSE
    MN_phi = Sol_phi !phi from previous vasp iteration, unscaled
 ENDIF

 !-------------------------------------------------------------
 !Leps returns  MN_Lp
 !MN_Lp = div ( epsilon * grad ) phi - kappa^2 * phi
 !in k space
 !Both kappa^2 and epsilon are in real space
 !Lp = - A*phi
 !-------------------------------------------------------------
 CALL Leps(MN_GRIDC, MN_LATT_CUR, MN_T_INFO, MN_phi, MN_d_eps, & 
      & MN_kappa2, gradphi_x_in, gradphi_y_in, gradphi_z_in, &
      & phi_work, MN_Lp)

 !-------------------------------------------------------------
 ! the residue
 !r = A*phi + (CHTOT + N)
 !-------------------------------------------------------------
 DO I = 1, MN_GRIDC%RC%NP
   MN_RESID(I) = MN_Lp(I) + MN_TOTN(I) !* EDEPS / MN_LATT_CUR%OMEGA
 ENDDO

 !-------------------------------------------------------------
 !preconditioning of the residue, MN_Z = invL(r)
 !z = precond * r
 !-------------------------------------------------------------
 MN_Z = 0._q

 DO I = 1, MN_GRIDC%RC%NP
   MN_Z(I) = REAL(MN_GSQU(I)) * MN_RESID(I)
 ENDDO

 !-------------------------------------------------------------
 !rinvLr = transpose(r) . z
 !-------------------------------------------------------------
 CALL MY_D_PROD(rinvLr, MN_RESID, MN_Z, MN_GRIDC)
 CALLMPI( M_sum_d(MN_GRIDC%COMM, rinvLr, 1))

 !-------------------------------------------------------------
 !square of the residue
 !r2 = r^2
 !-------------------------------------------------------------
 CALL MY_D_PROD(r2, MN_RESID, MN_RESID, MN_GRIDC)
 CALLMPI( M_sum_d(MN_GRIDC%COMM, r2, 1))

 r20 = r2

 MN_d = MN_Z

 !-------------------------------------------------------------
 !CG loop
 !-------------------------------------------------------------

 DO WHILE ( (count < 20000) .AND. (sqrt(r2) > EDIFFSOL).AND. (sqrt(rinvLr) > 1e-10) )

    CALL Leps(MN_GRIDC, MN_LATT_CUR, MN_T_INFO, MN_d, MN_d_eps, &
         & MN_kappa2, gradphi_x_in, gradphi_y_in, gradphi_z_in, &
         & phi_work, MN_Lp)

    ! MN_alpha = DOT_PRODUCT(MN_d,MN_Lp)
    CALL MY_D_PROD(MN_alpha, MN_d, MN_Lp, MN_GRIDC)

    CALLMPI( M_sum_d(MN_GRIDC%COMM, MN_alpha, 1))

    MN_alpha = 1._q/MN_alpha

    MN_alpha = -rinvLr*MN_alpha

    MN_phi = MN_phi + MN_alpha * MN_d

    MN_RESID = MN_RESID +  (MN_alpha*MN_Lp)

    DO I = 1, MN_GRIDC%RC%NP
       MN_Z(I) = (REAL(MN_GSQU(I))) * MN_RESID(I)
    ENDDO

    MN_beta = 1._q/rinvLr

    CALL MY_D_PROD(rinvLr, MN_RESID, MN_Z, MN_GRIDC)

    CALLMPI( M_sum_d(MN_GRIDC%COMM, rinvLr, 1))

    MN_beta = MN_beta * rinvLr

    MN_d = MN_Z + MN_beta * MN_d

    r2 = 0._q

    CALL MY_D_PROD(r2, MN_RESID, MN_RESID, MN_GRIDC)

    CALLMPI( M_sum_d(MN_GRIDC%COMM, r2, 1))

#ifdef debugsol
    IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) SQRT(rinvLr), SQRT(r2)
#endif
    count = count+1

 END DO
 !-------------------------------------------------------------
 !end of CG loop
 !-------------------------------------------------------------
 NCGSOL = count

#ifdef debugsol
 IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*)"CG loop count", count
#endif

 !-------------------------------------------------------------
 !save phi for the intial guess for the next iteration
 !-------------------------------------------------------------
 Sol_phi = MN_phi 

 !-------------------------------------------------------------
 !multiply by e/epsilon_0 and divide by volume for CHTOT
 !-------------------------------------------------------------
 MN_phi = MN_phi * EDEPS / MN_LATT_CUR%OMEGA

 !-------------------------------------------------------------
 !deallocate all the work arrays
 !-------------------------------------------------------------
 DEALLOCATE(MN_RESID,MN_Z,MN_Lp,MN_GSQU)
 DEALLOCATE(MN_d)
 DEALLOCATE(gradphi_x_in,gradphi_y_in,gradphi_z_in,phi_work)

 RETURN

END SUBROUTINE MINIMIZE


!>======================SUBROUTINE Leps===============================
!!
!!subroutine that evaluates the fourier transform of 
!!div(epsilon*grad phi)
!!Ax in Ax=b
!!
!!====================================================================
SUBROUTINE Leps(LPS_GRIDC, LPS_LATT_CUR, LPS_T_INFO, LPS_phi, &
     & LPS_d_eps, LPS_kappa2, gradphi_x_in, gradphi_y_in, &
     & gradphi_z_in, phi_work, LPS_Lp)

 IMPLICIT NONE

 TYPE (grid_3d) LPS_GRIDC
 TYPE (type_info) LPS_T_INFO
 TYPE (latt) LPS_LATT_CUR

 REAL(q),INTENT(IN) ::  LPS_d_eps(DIMREAL(LPS_GRIDC%MPLWV))
 REAL(q),INTENT(IN) ::  LPS_kappa2(DIMREAL(LPS_GRIDC%MPLWV))
 COMPLEX(q),INTENT(INOUT) ::  LPS_phi(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(OUT) :: LPS_Lp(LPS_GRIDC%MPLWV)
 REAL(q) :: RINPL,GX,GY,GZ

 COMPLEX(q),INTENT(INOUT) :: gradphi_x_in(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(INOUT) :: gradphi_y_in(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(INOUT) :: gradphi_z_in(LPS_GRIDC%MPLWV)

 COMPLEX(q),INTENT(INOUT) :: phi_work(LPS_GRIDC%MPLWV)

 INTEGER(q) :: N,I,J,k,N1,N2,N3,NC
      
 N = LPS_GRIDC%NPLWV

 !---------------------------------------------------------------------
 !Fourier transform of grad phi
 !---------------------------------------------------------------------
 DO I=1,LPS_GRIDC%RC%NP

    N1= MOD((I-1),LPS_GRIDC%RC%NROW) +1
    NC= (I-1)/LPS_GRIDC%RC%NROW+1
    N2= LPS_GRIDC%RC%I2(NC)
    N3= LPS_GRIDC%RC%I3(NC)

    GX = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(1,1)+LPS_GRIDC%LPCTY(N2)*&
         &LPS_LATT_CUR%B(1,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(1,3))
    GY = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(2,1)+LPS_GRIDC%LPCTY(N2)*&
         &LPS_LATT_CUR%B(2,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(2,3))
    GZ = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(3,1)+LPS_GRIDC%LPCTY(N2)*&
         &LPS_LATT_CUR%B(3,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(3,3))

    gradphi_x_in(I) = CITPI*GX*LPS_phi(I)
    gradphi_y_in(I) = CITPI*GY*LPS_phi(I)
    gradphi_z_in(I) = CITPI*GZ*LPS_phi(I)

 ENDDO 

 !---------------------------------------------------------------------
 !gradphi to real space
 !---------------------------------------------------------------------
 CALL FFT3D(gradphi_x_in,LPS_GRIDC,1) 
 CALL FFT3D(gradphi_y_in,LPS_GRIDC,1) 
 CALL FFT3D(gradphi_z_in,LPS_GRIDC,1) 

 !---------------------------------------------------------------------
 !phi to real space
 !---------------------------------------------------------------------
 CALL FFT3D(LPS_phi, LPS_GRIDC, 1) 

 !---------------------------------------------------------------------
 !multiply grad phi with epsilon in real space
 ! Also compute:
 !!phi_work = phi * eb_k * kappa^2, both phi and kappa2 in real space
 !--------------------------------------------------------------------
#ifdef realmode

 k = 1

 DO j = 1, DIMREAL(LPS_GRIDC%MPLWV), 2
    gradphi_x_in(k) = CMPLX( REAL(gradphi_x_in(k),KIND=q)*LPS_d_eps(j),&
         & AIMAG(gradphi_x_in(k))*LPS_d_eps(j+1),KIND=q)
    gradphi_y_in(k) = CMPLX( REAL(gradphi_y_in(k),KIND=q)*LPS_d_eps(j),&
         & AIMAG(gradphi_y_in(k))*LPS_d_eps(j+1),KIND=q)
    gradphi_z_in(k) = CMPLX( REAL(gradphi_z_in(k),KIND=q)*LPS_d_eps(j),&
         & AIMAG(gradphi_z_in(k))*LPS_d_eps(j+1),KIND=q)

    phi_work(k) = CMPLX( REAL(LPS_phi(k),KIND=q) * eb_k * LPS_kappa2(j),&
         & AIMAG(LPS_phi(k)) * eb_k * LPS_kappa2(j+1), KIND=q )

    k = k + 1

 ENDDO

#else

 DO j = 1,LPS_GRIDC%RL%NP
    gradphi_x_in(j) = gradphi_x_in(j)*LPS_d_eps(j)
    gradphi_y_in(j) = gradphi_y_in(j)*LPS_d_eps(j)
    gradphi_z_in(j) = gradphi_z_in(j)*LPS_d_eps(j)
    !phi_work = phi * eb_k * kappa^2, both phi and kappa2 in real space
    phi_work(j) = LPS_phi(j) * eb_k * LPS_kappa2(j)

 ENDDO
#endif
 !---------------------------------------------------------------------
 !epsilon*grad phi in kspace
 !---------------------------------------------------------------------
 CALL FFT3D(gradphi_x_in,LPS_GRIDC,-1)   
 CALL FFT3D(gradphi_y_in,LPS_GRIDC,-1)   
 CALL FFT3D(gradphi_z_in,LPS_GRIDC,-1)   

 CALL SETUNB(gradphi_x_in,LPS_GRIDC)   
 CALL SETUNB(gradphi_y_in,LPS_GRIDC)   
 CALL SETUNB(gradphi_z_in,LPS_GRIDC)   

 !---------------------------------------------------------------------
 !phi_work in kspace
 !---------------------------------------------------------------------
 CALL FFT3D(phi_work, LPS_GRIDC, -1)   

 CALL SETUNB(phi_work, LPS_GRIDC)   

 !---------------------------------------------------------------------
 !phi to k space
 !---------------------------------------------------------------------
 CALL FFT3D(LPS_phi, LPS_GRIDC, -1) 
 CALL SETUNB(LPS_phi, LPS_GRIDC)   

 !---------------------------------------------------------------------
 !scaling
 !---------------------------------------------------------------------
 DO I=1,LPS_GRIDC%RC%NP
   gradphi_x_in(I) = gradphi_x_in(I)/N
   gradphi_y_in(I) = gradphi_y_in(I)/N
   gradphi_z_in(I) = gradphi_z_in(I)/N

   phi_work(I) = phi_work(I)/N

   LPS_phi(I) = LPS_phi(I)/N

 ENDDO

 !---------------------------------------------------------------------
 !div(epsilon*grad phi) in kspace
 !add the kappa^2 contrib, -phi_work
 !---------------------------------------------------------------------
 DO I=1,LPS_GRIDC%RC%NP
   N1= MOD((I-1),LPS_GRIDC%RC%NROW) +1
   NC= (I-1)/LPS_GRIDC%RC%NROW+1
   N2= LPS_GRIDC%RC%I2(NC)
   N3= LPS_GRIDC%RC%I3(NC)

   GX = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(1,1)+LPS_GRIDC%LPCTY(N2)*&
        &LPS_LATT_CUR%B(1,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(1,3))
   GY = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(2,1)+LPS_GRIDC%LPCTY(N2)*&
        &LPS_LATT_CUR%B(2,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(2,3))
   GZ = (LPS_GRIDC%LPCTX(N1)*LPS_LATT_CUR%B(3,1)+LPS_GRIDC%LPCTY(N2)*&
        &LPS_LATT_CUR%B(3,2)+LPS_GRIDC%LPCTZ(N3)*LPS_LATT_CUR%B(3,3))

   LPS_Lp(I)= CITPI*(GX*gradphi_x_in(I) + GY*gradphi_y_in(I) + &
              &GZ*gradphi_z_in(I)) - phi_work(I)

 ENDDO
 
 END SUBROUTINE Leps


!>===========================SUBROUTINE  EPS_POT======================
!!
!!subroutine that evaluates the V
!!returns -(depsilon(n)/dn)*|grad phi|^2 
!! - eb_k*(dkappa^2(n)/dn)*(phi)^2
!! in fourier space
!! Both terms are computed in real space first and then FFT to k space
!!
!!====================================================================
SUBROUTINE  EPS_POT(EPP_GRIDC, EPP_LATT_CUR, EPP_T_INFO, EPP_P, &
     & EPP_PS_TOTN, EPP_phi, EPP_deps, EPP_kappa2, EPP_CWORK)

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN)   :: EPP_GRIDC
 TYPE (latt),INTENT(IN)      :: EPP_LATT_CUR
 TYPE (type_info),INTENT(IN) :: EPP_T_INFO
 TYPE(potcar),INTENT(IN)     :: EPP_P (EPP_T_INFO%NTYP)

 COMPLEX(q),INTENT(OUT) :: EPP_CWORK(EPP_GRIDC%MPLWV)
 RGRID,INTENT(IN)       :: EPP_PS_TOTN(DIMREAL(EPP_GRIDC%MPLWV))
 COMPLEX(q),INTENT(INOUT)  :: EPP_phi(EPP_GRIDC%MPLWV)
 REAL(q),INTENT(IN)     :: EPP_deps(EPP_GRIDC%MPLWV)
 REAL(q),INTENT(IN)     :: EPP_kappa2(EPP_GRIDC%MPLWV)

 REAL(q), ALLOCATABLE   :: eprime(:), kprime(:)
 REAL(q),ALLOCATABLE    :: phisq(:), RWORK(:)
 RGRID, ALLOCATABLE     :: RWORK_1(:)
 REAL(q)                :: EPP_CWORK_temp
 INTEGER(q)             :: I,J,K  

 ALLOCATE(eprime(EPP_GRIDC%RL%NP), kprime(EPP_GRIDC%RL%NP) )
 ALLOCATE(phisq(DIMREAL(EPP_GRIDC%MPLWV)))
 ALLOCATE(RWORK(DIMREAL(EPP_GRIDC%MPLWV)))
 ALLOCATE(RWORK_1(DIMREAL(EPP_GRIDC%MPLWV)))

 !---------------------------------------------------------------------
 !SHAPE_GRADN returns eprime, the derivative of shape function wrt 
 !the charge density in real space
 ! shapefunction = erfc = 1 - erf
 ! ==> -d erfc/dx = d erf/dx
 !eprime = d erf/dx, 
 !so the minus sign in the potenials are taken care of
 !---------------------------------------------------------------------
 CALL SHAPE_GRADN(EPP_GRIDC, EPP_T_INFO, EPP_P, EPP_LATT_CUR, &
       & EPP_PS_TOTN, eprime)
 !soltest
 !CALL WRITE_TO_FILE(EPP_GRIDC, EPP_LATT_CUR, EPP_T_INFO, 'EPRIME', eprime)
 !soltest
 kprime = eprime * kappa2b_k
 eprime = eprime * (eb_k - 1._q)

 !---------------------------------------------------------------------
 !EPP_phi,the solution to the generalized poisson equation,
 ! is in reciprocal space, 
 !NABLAN returns |\nabla phi| ^2 in real space
 !---------------------------------------------------------------------
 CALL NABLAN(EPP_phi,EPP_GRIDC,EPP_LATT_CUR,phisq)

 !---------------------------------------------------------------------
 !phi to real space
 !---------------------------------------------------------------------
 CALL FFT3D(EPP_phi, EPP_GRIDC, 1) 

 CALL UNPACK_C2R(EPP_phi, 1.0_q, EPP_phi, 0.0_q, RWORK_1, EPP_GRIDC)

#ifdef debugsolwf
 IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing PHI_PBZ'
 CALL WRITE_TO_FILE(EPP_GRIDC, EPP_LATT_CUR, EPP_T_INFO, 'PHI_PBZ', &
      & RWORK_1)
#endif

 !---------------------------------------------------------------------
 ! eps potential in real space
 !---------------------------------------------------------------------
 EPP_CWORK_temp = 0._q

 !---------------------------------------------------------------------
 !multiply eprime and phisq in real space
 !add kappa^2 term, eb_k * kprime * phi^2 in real space
 !---------------------------------------------------------------------
 RWORK = 0._q

 DO I = 1, EPP_GRIDC%RL%NP
    EPP_CWORK_temp = eprime(I) * phisq(I) + &
         & eb_k * kprime(I) * RWORK_1(I) * RWORK_1(I)
    RWORK(I) = EPP_CWORK_temp
 ENDDO

 !---------------------------------------------------------------------
 !phi to k space
 !---------------------------------------------------------------------
 CALL FFT3D(EPP_phi, EPP_GRIDC, -1) 

 !---------------------------------------------------------------------
 !scale phi
 !---------------------------------------------------------------------
 DO I=1,EPP_GRIDC%RC%NP
   EPP_phi(I) = EPP_phi(I)/EPP_GRIDC%NPLWV
 ENDDO

 CALL PACK_R2C(RWORK, EPP_CWORK, EPP_GRIDC)
   
 !---------------------------------------------------------------------
 !take the eps potential to reciprocal space
 !---------------------------------------------------------------------
 CALL FFT3D(EPP_CWORK,EPP_GRIDC,-1)
 CALL SETUNB(EPP_CWORK,EPP_GRIDC)
 EPP_CWORK = EPP_CWORK/EPP_GRIDC%NPLWV

 DEALLOCATE(eprime, kprime, phisq, RWORK, RWORK_1)
 
 RETURN

END SUBROUTINE EPS_POT


!>======================SUBROUTINE SHAPE_GRADN========================
!!
!!calculates first derivative of the shape function  wrt CHTOT 
!!in realspace
!!
!!====================================================================
SUBROUTINE SHAPE_GRADN(EPR_GRIDC, EPR_T_INFO, EPR_P, EPR_LATT_CUR, &
     & EPR_PS_TOTN, eprime)

 IMPLICIT NONE

 !dummy variables
 TYPE(grid_3d),INTENT(IN)   :: EPR_GRIDC
 TYPE(type_info),INTENT(IN) :: EPR_T_INFO
 TYPE(potcar),INTENT(IN)    :: EPR_P(EPR_T_INFO%NTYP)
 TYPE(latt),INTENT(IN)      :: EPR_LATT_CUR
 ! COMPLEX(q),INTENT(IN) ::  EPR_PS_TOTN(EPR_GRIDC%MPLWV)
 RGRID,INTENT(IN)           ::  EPR_PS_TOTN(DIMREAL(EPR_GRIDC%MPLWV))
 REAL(q),INTENT(OUT) ::  eprime(EPR_GRIDC%RL%NP)

 !local variables
 REAL(q) ::EPR_C,eprime_temp,EPR_Z
 INTEGER I 

 EPR_C = 1._q / sqrt(2*PI) / sigma_k

 DO I=1, EPR_GRIDC%RL%NP
    EPR_Z = log(MAX(REAL(EPR_PS_TOTN(I),KIND=q),CHGMIN)/nc_k) / sqrt(2._q) / sigma_k
    eprime(I) = EPR_C * (EXP(-(EPR_Z)**2)) / (MAX(REAL(EPR_PS_TOTN(I),KIND=q),CHGMIN))
 END DO

 RETURN

END SUBROUTINE SHAPE_GRADN


!>==========================SUBROUTINE NABLAN=========================
!!
!! input : D_PHI , in reciprocal space, unscaled
!! calculate |\nabla D_PHI |^2 = D_PHI2 in real space
!! optional arg, D_grad = \nabla D_PHI, real space vector
!!
!!====================================================================
SUBROUTINE NABLAN(D_PHI, D_GRIDC, D_LATT_CUR, D_PHI2, D_grad)

 IMPLICIT NONE

!dummy variables
 TYPE(grid_3d),INTENT(IN)      :: D_GRIDC
 TYPE(latt),INTENT(IN)         :: D_LATT_CUR
 COMPLEX(q),INTENT(IN)         :: D_PHI(D_GRIDC%MPLWV)
 REAL(q),INTENT(OUT)           :: D_PHI2(DIMREAL(D_GRIDC%MPLWV))
 REAL(q), OPTIONAL,INTENT(OUT) :: D_grad(3, D_GRIDC%RL%NP)
!local variables
 COMPLEX(q),ALLOCATABLE        :: CR(:)  
 RGRID,ALLOCATABLE             :: RWORK(:)

 INTEGER(q)                    :: I,N1,N2,N3,NC
 REAL(q)                       :: GX,GY,GZ,tmp
 
 ALLOCATE( CR(D_GRIDC%MPLWV), RWORK(DIMREAL(D_GRIDC%MPLWV)) )

 D_PHI2 = 0._q

 !---------------------------------------------------------------------
 ! x-component:
 !---------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GX = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(1,1)+D_GRIDC%LPCTY(N2)*&
        & D_LATT_CUR%B(1,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(1,3))

   CR(I) = D_PHI(I) * GX * CITPI
 ENDDO
 !---------------------------------------------------------------------
 ! gradphi_x to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
   tmp = REAL(RWORK(I),KIND=q)
   D_PHI2(I) = D_PHI2(I) + tmp * tmp
   IF (PRESENT(D_grad)) D_grad(1,I) = tmp
 ENDDO
 !--------------------------------------------------------------------
 ! y-component:
 !--------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GY = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(2,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(2,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(2,3))

   CR(I) = D_PHI(I) * GY * CITPI
 ENDDO

 !---------------------------------------------------------------------
 ! grad_y to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
   tmp = REAL(RWORK(I),KIND=q)
   D_PHI2(I) = D_PHI2(I) + tmp * tmp
   IF (PRESENT(D_grad)) D_grad(2,I) = tmp
 ENDDO
 !--------------------------------------------------------------------
 ! z-component:
 !--------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GZ = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(3,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(3,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(3,3))

   CR(I)=D_PHI(I) * GZ * CITPI
 ENDDO

 !---------------------------------------------------------------------
 ! grad_z to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
   tmp = REAL(RWORK(I),KIND=q)
   D_PHI2(I) = D_PHI2(I) + tmp * tmp
   IF (PRESENT(D_grad)) D_grad(3,I) = tmp
  ENDDO

 DEALLOCATE(CR,RWORK)

 RETURN

END SUBROUTINE NABLAN


!>==========================SUBROUTINE LAPLN==========================
!!
!! calculate the Laplacian of D_PHI in real space
!! D_PHI, the input in reciprocal space
!! D_PHI2, the output in real space
!!
!!====================================================================
SUBROUTINE LAPLN(D_PHI,D_GRIDC,D_LATT_CUR,D_PHI2)

 IMPLICIT NONE

 !dummy variables
 TYPE(grid_3d),INTENT(IN) ::  D_GRIDC
 TYPE(latt),INTENT(IN)    :: D_LATT_CUR
 COMPLEX(q),INTENT(IN)    :: D_PHI(D_GRIDC%MPLWV)
 REAL(q),INTENT(OUT)      :: D_PHI2(DIMREAL(D_GRIDC%MPLWV))
 !local variables
 COMPLEX(q),ALLOCATABLE   :: CR(:)
 RGRID, ALLOCATABLE       :: RWORK(:)
 INTEGER(q)               :: I,N1,N2,N3,NC
 REAL(q)                  :: GX,GY,GZ
 
 ALLOCATE(CR(D_GRIDC%MPLWV), RWORK(DIMREAL(D_GRIDC%MPLWV)))

 D_PHI2 = 0._q

 !---------------------------------------------------------------------
 ! x-component:
 !---------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GX = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(1,1)+D_GRIDC%LPCTY(N2)*&
       &D_LATT_CUR%B(1,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(1,3))

   CR(I) = D_PHI(I) * (GX*CITPI)**2.0
 ENDDO

 !---------------------------------------------------------------------
 ! gradphi to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
    D_PHI2(I) = D_PHI2(I) + REAL(RWORK(I),KIND=q)
 ENDDO

 !---------------------------------------------------------------------
 ! y-component:
 !---------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GY = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(2,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(2,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(2,3))

   CR(I) = D_PHI(I) * (GY*CITPI)**2
 ENDDO

 !---------------------------------------------------------------------
 ! grad_y to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
    D_PHI2(I) = D_PHI2(I) + REAL(RWORK(I),KIND=q)
 ENDDO
 !---------------------------------------------------------------------
 ! z-component:
 !---------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GZ = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(3,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(3,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(3,3))

   CR(I)=D_PHI(I) * (GZ*CITPI)**2.0
 ENDDO

 !---------------------------------------------------------------------
 ! grad_z to real space:
 !---------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
  D_PHI2(I) = D_PHI2(I) + REAL(RWORK(I),KIND=q)
 ENDDO

 DEALLOCATE(CR,RWORK)

 RETURN

END SUBROUTINE LAPLN


!>======================SUBROUTINE SimplePoisson=================
!!
!!returns n+N/G^2
!!
!!===============================================================
SUBROUTINE SimplePoisson(GC_GRIDC,GC_LATT_CUR,GC_IN,GC_OUT,GFLAG)

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: GC_GRIDC
 TYPE (latt),INTENT(IN) :: GC_LATT_CUR
 COMPLEX(q),INTENT(OUT) ::  GC_OUT(GC_GRIDC%MPLWV)
 COMPLEX(q),INTENT(IN) :: GC_IN(GC_GRIDC%MPLWV)  
 REAL(q) GX,GY,GZ,GSQU,SCALE
 INTEGER I,NC,N1,N2,N3,GFLAG

 SCALE = EDEPS/GC_LATT_CUR%OMEGA/TPI**2
 I = 0

 col: DO NC=1,GC_GRIDC%RC%NCOL

    N2= GC_GRIDC%RC%I2(NC)
    N3= GC_GRIDC%RC%I3(NC)

    row: DO N1=1,GC_GRIDC%RC%NROW
       I=I+1

       GX = (GC_GRIDC%LPCTX(N1)*GC_LATT_CUR%B(1,1)+&
            &GC_GRIDC%LPCTY(N2)*GC_LATT_CUR%B(1,2)+&
            &GC_GRIDC%LPCTZ(N3)*GC_LATT_CUR%B(1,3))
       GY = (GC_GRIDC%LPCTX(N1)*GC_LATT_CUR%B(2,1)+&
            &GC_GRIDC%LPCTY(N2)*GC_LATT_CUR%B(2,2)+&
            &GC_GRIDC%LPCTZ(N3)*GC_LATT_CUR%B(2,3))
       GZ = (GC_GRIDC%LPCTX(N1)*GC_LATT_CUR%B(3,1)+&
            &GC_GRIDC%LPCTY(N2)*GC_LATT_CUR%B(3,2)+&
            &GC_GRIDC%LPCTZ(N3)*GC_LATT_CUR%B(3,3))
       
       GSQU=GX**2+GY**2+GZ**2
       
       IF (GFLAG == 1) THEN

          IF ((GC_GRIDC%LPCTX(N1)==0).AND.&
               &(GC_GRIDC%LPCTY(N2)==0).AND.&
               &(GC_GRIDC%LPCTZ(N3)==0)) THEN
             
             GC_OUT(I)=(0.0_q,0.0_q)
          ELSE
             GC_OUT(I)=(GC_IN(I))/GSQU*SCALE
          ENDIF

       ELSEIF (GFLAG == 0) THEN
          GC_OUT(I)=(GC_IN(I))*GSQU/SCALE         
       ENDIF

    ENDDO row
 ENDDO col

 RETURN

END SUBROUTINE SimplePoisson


!>==========================SUBROUTINE MY_D_PROD====================
!!
!! helper routine for dot product
!! also used in pot_k.F
!!
!!==================================================================
SUBROUTINE MY_D_PROD(D_P_val,ARRAY1,ARRAY2,GRIDC)

 IMPLICIT NONE

 TYPE(grid_3d),INTENT(IN) :: GRIDC

 COMPLEX(q),INTENT(IN) ::  ARRAY1(GRIDC%RC%NP)
 COMPLEX(q),INTENT(IN) :: ARRAY2(GRIDC%RC%NP)
 REAL(q),INTENT(OUT) :: D_P_val
 REAL(q) :: FACTM
 INTEGER :: I, NC, N1, N2, N3

 D_P_val = 0._q
 I = 1
 col: DO NC=1,GRIDC%RC%NCOL
    N2= GRIDC%RC%I2(NC)
    N3= GRIDC%RC%I3(NC)
    row: DO N1=1,GRIDC%RC%NROW
       SETFACT1
       SETFACT
       D_P_val = D_P_val + MULFACT (ARRAY1(I)) * CONJG(ARRAY2(I))
       I=I+1
    ENDDO row
 ENDDO col

 RETURN

END SUBROUTINE MY_D_PROD


!>==========================SUBROUTINE PACK_R2C======================
!!
!!  helper routine to pack real array(RA) to a complex array(CB)
!!
!!===================================================================
SUBROUTINE PACK_R2C(RA,CB,GRIDC)

 IMPLICIT NONE

 TYPE(grid_3d),INTENT(IN) :: GRIDC

 REAL(q)    :: RA(DIMREAL(GRIDC%MPLWV))
 COMPLEX(q) :: CB(GRIDC%MPLWV)
 INTEGER(q) :: I,J

#ifdef realmode
 J = 1
 DO I= 1,DIMREAL(GRIDC%MPLWV),2
   CB(J) = CMPLX(RA(I),RA(I+1))
   J = J + 1
 ENDDO
#else
 DO I= 1,DIMREAL(GRIDC%MPLWV)
  CB(I) = CMPLX(RA(I))
 ENDDO
#endif

 RETURN

END SUBROUTINE PACK_R2C


!>==========================SUBROUTINE UNPACK_C2R====================
!!
!!  helper routine to unpack and maybe scale
!!  complex array(A & B) to real array(C)
!!
!!===================================================================
SUBROUTINE UNPACK_C2R(A,SCALE1,B,SCALE2,C,GRID)

 IMPLICIT NONE

 TYPE (grid_3d) :: GRID
 COMPLEX(q) :: A(GRID%MPLWV),B(GRID%MPLWV)
 RGRID  :: C(DIMREAL(GRID%MPLWV))
 REAL(q) :: SCALE1, SCALE2
 INTEGER(q) :: NP, I ,J

 NP = GRID%MPLWV
 
 J = 1
 DO I=1,NP
#ifdef realmode
    C(J) = REAL(A(I),KIND=q)* SCALE1 + REAL(B(I),KIND=q)* SCALE2
    C(J+1) = AIMAG(A(I))* SCALE1 + AIMAG(B(I))* SCALE2 
    J = J + 2
#else
    C(I) = A(I) * SCALE1 + B(I) * SCALE2
#endif
 ENDDO

 RETURN

END SUBROUTINE UNPACK_C2R


!>==========================SUBROUTINE WRITE_TO_FILE==================
!!
!! write real space quantities to file
!! writes poscar first and then CVTOT
!! could be any quantity defined on RL grid
!!
!!====================================================================
SUBROUTINE WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, FNAME, CVTOT)

 USE prec
 USE mgrid
 USE poscar

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: GRIDC
 TYPE(type_info),INTENT(IN) :: T_INFO
 TYPE (latt),INTENT(IN) :: LATT_CUR
 RGRID CVTOT(GRIDC%RL%NP)

 LOGICAL LLONG
 CHARACTER (40) FORM, SZNAM
 CHARACTER (*) FNAME
 REAL(q),ALLOCATABLE ::  WORK(:)
 INTEGER :: ISTAT, NODE_ME, IONODE, IU
 INTEGER :: NALLOC, NZ, NWRITE, NWRITTEN, N
 !REAL(q),ALLOCATABLE ::  AVG_XY(:)

 NODE_ME = 0
 IONODE = 0
 IU = 468
 SZNAM = 'Test'
 FORM = '(1(1X,E17.11))'
 NWRITE = 5
 NALLOC = GRIDC%NGX*GRIDC%NGY

#ifdef MPI
 NODE_ME = GRIDC%COMM%NODE_ME
 IONODE = GRIDC%COMM%IONODE
#endif

 ALLOCATE(WORK(NALLOC), STAT=ISTAT)
 !ALLOCATE(AVG_XY(GRIDC%NGZ))

 io_begin
     OPEN(IU, FILE = FNAME, STATUS = 'UNKNOWN')
     REWIND IU
     CALL OUTPOS(IU, .FALSE., SZNAM, T_INFO, LATT_CUR%SCALE, &
          & LATT_CUR%A, .FALSE., T_INFO%POSION)
 io_end

 IF (ISTAT>0) RETURN ! can not write the potential immediate exit

 IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
    WRITE(*,*)'internal ERROR: OUTPOT GRIDC%NPLWV is not compatible',&
         &   ' with  GRIDC%NGX,GRIDC%NGY,NGZC'
    WRITE(*,*)'   ',GRIDC%NPLWV,GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ
    STOP
 ENDIF

 do_io WRITE(IU,'(3I5)') GRIDC%NGX, GRIDC%NGY, GRIDC%NGZ

 NWRITTEN = 0
 !AVG_XY = 0

 DO NZ = 1, GRIDC%NGZ
    CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CVTOT, NZ)
    io_begin
    DO N = 1, NALLOC
       !AVG_XY(NZ) = AVG_XY(NZ) + WORK(N)/NALLOC
       NWRITTEN=NWRITTEN+1
       IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
          IF (ABS(WORK(N)) < 1E-10 ) THEN
             WRITE(IU, FORM) 0._q
          ELSE
             WRITE(IU, FORM) WORK(N)
          ENDIF
       ELSE
          IF (ABS(WORK(N)) < 1E-10 ) THEN
             WRITE(IU, FORM, ADVANCE='NO') 0._q
          ELSE
             WRITE(IU, FORM, ADVANCE='NO') WORK(N)
          ENDIF
       ENDIF
    ENDDO
    io_end
 ENDDO

 IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

 !IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) 'MAX AVG_XY = ', MAXVAL(AVG_XY)
 !IF (SP_IO%IU0>=0) WRITE(SP_IO%IU0,*) 'AVG_XY averaged along vacuum z = ', SUM(AVG_XY(GRIDC%NGZ-9:GRIDC%NGZ)) / 10.

 !DEALLOCATE(AVG_XY)

 DEALLOCATE(WORK)

 do_io CLOSE(IU)

 RETURN

END SUBROUTINE WRITE_TO_FILE


!>========================== SUBROUTINE GET_VAC_LEVEL ================
!!
!! Get the vacuum level
!! uses the AVERAG subroutine from dipol.F
!! both CHTOT and CVTOT must in real space
!!
!!====================================================================
SUBROUTINE GET_VAC_LEVEL(GRIDC, CHTOT, CVTOT, VACPOT)

 USE prec
 USE mgrid
 USE poscar

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: GRIDC

 RGRID CVTOT(DIMREAL(GRIDC%MPLWV))
 RGRID CHTOT(DIMREAL(GRIDC%MPLWV))

 REAL(q),ALLOCATABLE :: DENLIN(:)

 REAL(q) :: DENMIN , DEN, VACPOT
 INTEGER(q) :: DIRECTION, IND, NOUT, II, INDMINDEN, II_IND
 INTEGER,PARAMETER :: NWIN = 3
 REAL(q) :: WINDOW(-NWIN:NWIN) = (/0.1,0.2,0.3,0.4,0.3,0.2,0.1/)

 ALLOCATE(DENLIN(MAX(GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ)) )
 
 DENMIN = 1E9
 !x-->1, y-->2, z-->3
 DIRECTION = 3
 NOUT = GRIDC%NGPTAR(DIRECTION) 

 DENLIN = 0._q
 !average chargedensity along direction, DIRECTION
 CALL AVERAG(GRIDC, DIRECTION, CHTOT, DENLIN, .TRUE.)

 !set vacuum to be the point in space where the electronic charge 
 !density is the minimum 
 !get the index corresponding to the minimum charge density
 DO IND = 1, NOUT
    !use a sliding average to smoothen in charge density
    DEN = 0
    DO II = -NWIN, NWIN
       II_IND = MOD(IND+II-1+NOUT, NOUT) + 1
       DEN = DEN + DENLIN(II_IND) * WINDOW(II)
    ENDDO
    IF (DEN < DENMIN) THEN
       DENMIN = DEN
       INDMINDEN = IND
    ENDIF
 ENDDO

 DENLIN = 0._q
 !average potential along direction, DIRECTION
 CALL AVERAG(GRIDC, DIRECTION, CVTOT, DENLIN, .FALSE.)

 VACPOT = DENLIN(MOD(INDMINDEN + NOUT, NOUT))
 !VACPOT=DENLIN(MOD(INDMIN+INT(WIDTH)+NOUT,NOUT))

#ifdef debugsol
 IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) 'DENMIN = ', DENMIN
 IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) 'INDMINDEN = ', &
      & INDMINDEN, MOD(INDMINDEN + NOUT, NOUT)
#endif

 DEALLOCATE(DENLIN)

 RETURN

END SUBROUTINE GET_VAC_LEVEL


!>========================== SUBROUTINE COMPUTE_RHOB ================
!!
!! bound charge in real space
!! phi in rec space
!!
!!====================================================================
SUBROUTINE COMPUTE_RHOB(GRIDC, LATT_CUR, T_INFO, EPSILON, PHI, RHOB)

  IMPLICIT NONE

  TYPE(grid_3d),INTENT(IN) :: GRIDC
  TYPE (latt),INTENT(IN) :: LATT_CUR
  TYPE(type_info),INTENT(IN) :: T_INFO

  COMPLEX(q) :: PHI(GRIDC%MPLWV)
  REAL(q) :: EPSILON(DIMREAL(GRIDC%MPLWV)), RHOB(DIMREAL(GRIDC%MPLWV))

  COMPLEX(q),ALLOCATABLE :: gradphi_x_in(:)
  COMPLEX(q),ALLOCATABLE :: gradphi_y_in(:)
  COMPLEX(q),ALLOCATABLE :: gradphi_z_in(:)
  COMPLEX(q),ALLOCATABLE :: phi_work(:), C_RHOB(:)

  REAL(q),ALLOCATABLE :: DUMMY_KAPPA2(:)

  ALLOCATE(C_RHOB(GRIDC%MPLWV))
  ALLOCATE(DUMMY_KAPPA2(DIMREAL(GRIDC%MPLWV)))
  ALLOCATE(gradphi_x_in(GRIDC%MPLWV),gradphi_y_in(GRIDC%MPLWV))
  ALLOCATE(gradphi_z_in(GRIDC%MPLWV), PHI_WORK(GRIDC%MPLWV))

  DUMMY_KAPPA2 = 0._q
  !epsilon -1 
  EPSILON = EPSILON - 1._q
  
  !div(epsilon grap phi)
  CALL Leps(GRIDC, LATT_CUR, T_INFO, PHI, EPSILON, & 
      & DUMMY_KAPPA2, gradphi_x_in, gradphi_y_in, gradphi_z_in, &
      & phi_work, C_RHOB)

  !rhob in real space
  CALL FFT3D(C_RHOB, GRIDC, 1) 
  CALL UNPACK_C2R(C_RHOB, 1.0_q, C_RHOB, 0.0_q, RHOB, GRIDC)

  !write to file
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing RHOB'
  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'RHOB', RHOB)

  DEALLOCATE(C_RHOB, DUMMY_KAPPA2)
  DEALLOCATE(gradphi_x_in, gradphi_y_in)
  DEALLOCATE(gradphi_z_in, PHI_WORK)

  RETURN

END SUBROUTINE COMPUTE_RHOB


!>========================== SUBROUTINE COMPUTE_RHOION ===============
!!
!! ionic charge in real space
!! phi in rec space
!!
!!====================================================================
SUBROUTINE COMPUTE_RHOION(GRIDC, LATT_CUR, T_INFO, KAPPA2, PHI,RHOION)

  IMPLICIT NONE

  TYPE(grid_3d),INTENT(IN) :: GRIDC
  TYPE (latt),INTENT(IN) :: LATT_CUR
  TYPE(type_info),INTENT(IN) :: T_INFO

  COMPLEX(q) :: PHI(GRIDC%MPLWV)
  REAL(q) :: KAPPA2(DIMREAL(GRIDC%MPLWV))
  REAL(q) :: RHOION(DIMREAL(GRIDC%MPLWV))

  REAL(q),ALLOCATABLE :: R_PHI(:)

  INTEGER(q) :: I

  ALLOCATE(R_PHI(DIMREAL(GRIDC%MPLWV)))

  !phi to real space
  CALL FFT3D(PHI, GRIDC, 1) 
  CALL UNPACK_C2R(PHI, 1.0_q, PHI, 0.0_q, R_PHI, GRIDC)
  
  ! rhoion = - kappa^2 * phi
  ! in real space
  DO I = 1, GRIDC%RL%NP
     RHOION(I) = -KAPPA2(I) * R_PHI(I)
  ENDDO

  !write to file
  IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing RHOION'
  CALL WRITE_TO_FILE(GRIDC, LATT_CUR, T_INFO, 'RHOION', RHOION)

  DEALLOCATE(R_PHI)

  RETURN 

END SUBROUTINE COMPUTE_RHOION


END MODULE POT_LPCM_K
