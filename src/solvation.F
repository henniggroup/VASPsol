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


MODULE POT_LPCM_K

 USE prec
 USE constant 
 USE base
 USE mpimy 
 USE poscar, ONLY : type_info 
 USE lattice, ONLY : latt, DIRKAR 
 USE mgrid, ONLY : grid_3d
 USE pseudo, ONLY : potcar
#ifdef sol_compat
 USE c2f_interface, ONLY : ERRFC
#endif

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
 RGRID,ALLOCATABLE :: SP_d_eps(:), SP_kappa2(:)
 RGRID,ALLOCATABLE :: RWORK(:)
 INTEGER :: I,J, NCGSOL
 REAL(q)  :: Atot, Ael, Acav, Ecorr1, Ecorr2, Ecorr_cav
 REAL(q) :: sigma_nc_k1,sigma_k1,nc_k1,eb_k1,kappa2b_k1
 REAL(q) :: tau1,ediffs
 REAL(q), INTENT(IN) :: sigma_rc_k1(SP_T_INFO%NTYP), CORE_C1(SP_T_INFO%NTYP)
 !vacuum level
 RGRID,ALLOCATABLE :: RWORK_PHI(:), RWORK_CHG(:)
 COMPLEX(q),ALLOCATABLE :: CWORK_K(:)
 !vacuum level
 RGRID,ALLOCATABLE :: RHOB(:), RHOION(:)

 ALLOCATE(SP_d_eps(DIMREAL(SP_GRIDC%MPLWV)) )
 ALLOCATE(SP_kappa2(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(phi_tilda(SP_GRIDC%MPLWV), phi(SP_GRIDC%MPLWV))
 ALLOCATE(Vcav(SP_GRIDC%MPLWV), RWORK(DIMREAL(SP_GRIDC%MPLWV)))
 !vacuum level
 ALLOCATE(RWORK_PHI(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(RWORK_CHG(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(CWORK_K(SP_GRIDC%MPLWV))
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
 !vacuum level

 PK_COUNTER = PK_COUNTER + 1
      
 DEALLOCATE(SP_d_eps, SP_kappa2, phi_tilda, phi, Vcav, RWORK)
 !vacuum level
 DEALLOCATE(RWORK_PHI)
 DEALLOCATE(RWORK_CHG)
 DEALLOCATE(CWORK_K)
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
 RGRID,INTENT(OUT) :: CE_d_eps(DIMREAL(CE_GRIDC%MPLWV))
 RGRID,INTENT(OUT) :: CE_kappa2(DIMREAL(CE_GRIDC%MPLWV))

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
#ifndef sol_compat
 ! compimentary error function implemented in errf.F in vasp.5.lib      
 REAL(q), EXTERNAL :: ERRFC
#endif
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

#ifdef debugsol
 IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing EPS'
 CALL WRITE_TO_FILE(CE_GRIDC, CE_LATT_CUR, CE_T_INFO, 'EPS', CE_d_eps)
 IF (SP_IO%IU0>=0)  WRITE(SP_IO%IU0,*) 'Writing KAPPA2'
 CALL WRITE_TO_FILE(CE_GRIDC, CE_LATT_CUR, CE_T_INFO, 'KAPPA2', CE_kappa2)
#endif 

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
 RGRID,INTENT(IN) ::  MN_d_eps(DIMREAL(MN_GRIDC%MPLWV))
 RGRID,INTENT(IN) ::  MN_kappa2(DIMREAL(MN_GRIDC%MPLWV))
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

 RGRID,INTENT(IN) ::  LPS_d_eps(DIMREAL(LPS_GRIDC%MPLWV))
 RGRID,INTENT(IN) ::  LPS_kappa2(DIMREAL(LPS_GRIDC%MPLWV))
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
 RGRID,INTENT(IN)     :: EPP_deps(EPP_GRIDC%MPLWV)
 RGRID,INTENT(IN)     :: EPP_kappa2(EPP_GRIDC%MPLWV)

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
  RGRID :: EPSILON(DIMREAL(GRIDC%MPLWV)), RHOB(DIMREAL(GRIDC%MPLWV))

  COMPLEX(q),ALLOCATABLE :: gradphi_x_in(:)
  COMPLEX(q),ALLOCATABLE :: gradphi_y_in(:)
  COMPLEX(q),ALLOCATABLE :: gradphi_z_in(:)
  COMPLEX(q),ALLOCATABLE :: phi_work(:), C_RHOB(:)

  RGRID,ALLOCATABLE :: DUMMY_KAPPA2(:)

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
  RGRID :: KAPPA2(DIMREAL(GRIDC%MPLWV)), RHOION(DIMREAL(GRIDC%MPLWV))

  RGRID,ALLOCATABLE :: R_PHI(:)

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
 PUBLIC :: GET_FERMISHIFT
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
! vaclevel

END MODULE POT_K

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
  PUBLIC :: GET_FERMISHIFT, LSOL
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

!******************** SUBROUTINE SOL_READER ***************************
!
!
! Reads in the solvation model parameters
!
!      
!**********************************************************************
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

       !---------------------------------------------------------------
       !sigma_k
       !--------------------------------------------------------------
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

       !---------------------------------------------------------------
       !nc_k
       !---------------------------------------------------------------
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

       !---------------------------------------------------------------
       !eb_k
       !---------------------------------------------------------------
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

       !---------------------------------------------------------------
       !tau
       !---------------------------------------------------------------
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

       !---------------------------------------------------------------
       !lambda_d_k = debye screening length in angstrom : the user 
       !input from INCAR.
       !kappa2b_k = (1/debye screening length)^2 
       ! if lambd_d_k not set in INCAR, kappa2b_k will be set to 
       ! zero(default value)
       !--------------------------------------------------------------
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

    ENDIF ! close lsol
      
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
             '   LRHOB     =',L6,  '     write boundcharge  density' /&
             '   LRHOION     =',L6,  '     write ioniccharge density' /)


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
