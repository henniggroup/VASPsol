!***********************************************************************
! Copyright 2014
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

MODULE POT_NLN_K

 USE pot_lpcm_k

 PUBLIC :: VCORRECTION_NLN_PCM
 PRIVATE :: SET_PARAMS_NLN, COMPUTE_EFIELD_NLN, COMPUTE_LANG_NLN 
 PRIVATE :: COMPUTE_CHI_NLN, COMPUTE_EPSILON_NLN, COMPUTE_RHOB_NLN, GET_V_E_NLN
 PRIVATE :: MINIMIZE_NLN

 REAL(q), PRIVATE, SAVE                 :: p0_k, eps_infty_k, Nsol_k, kb_k, T_k
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: NDIM_EFIELD(:) !non dimensional electri field in real space
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: R_PS_TOTN(:) !ps_totn in real space

 CONTAINS

!>========================== SUBROUTINE VCORRECTION_NLN_PCM ======================
!!
!!
!!
!!===============================================================================

  SUBROUTINE Vcorrection_nln_pcm(SP_IO, SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, &
                                 & TOTN, PS_TOTN, Vcorr, Vdiel, Adiel )

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: SP_T_INFO
 TYPE(potcar),INTENT(IN) :: SP_P (SP_T_INFO%NTYP)
 TYPE (latt),INTENT(IN) :: SP_LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: SP_GRIDC
 TYPE (in_struct),INTENT(IN) :: IO

 COMPLEX(q), INTENT(OUT) ::  Vcorr(SP_GRIDC%MPLWV)
 COMPLEX(q),ALLOCATABLE :: phi_tilda(:),phi(:),Vcav(:)
 COMPLEX(q),INTENT(OUT) :: Vel(SP_GRIDC%MPLWV)
 COMPLEX(q),INTENT(IN) ::TOTN(SP_GRIDC%MPLWV),PS_TOTN(SP_GRIDC%MPLWV) 
 REAL(q),ALLOCATABLE :: SP_d_eps(:)
 RGRID,ALLOCATABLE :: RWORK(:)
 INTEGER :: I,J, NCGSOL
 REAL(q)  :: Atot, Ael, Acav, Ecorr1, Ecorr2, Ecorr_cav
 REAL(q) :: sigma_nc_k1,sigma_k1,nc_k1,eb_k1,CORE_C1,tau1,ediffs
 REAL(q), INTENT(IN) :: sigma_rc_k1(SP_T_INFO%NTYP)

 ALLOCATE(SP_d_eps(DIMREAL(SP_GRIDC%MPLWV)),RWORK(DIMREAL(SP_GRIDC%MPLWV)))
 ALLOCATE(phi_tilda(SP_GRIDC%MPLWV),phi(SP_GRIDC%MPLWV))
 ALLOCATE(Vcav(SP_GRIDC%MPLWV))
      
 IF (PK_COUNTER == 1)  THEN
   ALLOCATE(Sol_phi(SP_GRIDC%MPLWV))
   Sol_phi = 0._q
 ENDIF

 Ecorr1 = 0.0_q
 Ecorr2 = 0.0_q
 phi_tilda = 0._q
 phi = 0._q
 Vel = 0._q
 Vcorr = 0._q
 Vcav = 0._q
 SP_IO = IO


 CALL SET_PARAMS_NLN(p0_k1, eps_infty_k1, Nsol_k1, T_k1)

! calls compute_epsilon_nln and minimize(from pot_lpcm_k module)
 CALL MINIMIZE_NLN(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, PS_TOTN, TOTN, SP_d_eps, phi_tilda, NCGSOL) 

! returns the potential and the energy contributoins, sets Vcorr and Atot
 CALL GET_V_E_NLN(phi_tilda, phi, Vcorr, Atot) 


 RETURN

END SUBROUTINE Vcorrection_nln_pcm




!>========================== SUBROUTINE SET_PARAMS_NLN =================
!!
!! Sets the params specifi to the nonlinear model
!!
!!=====================================================================
SUBROUTINE  SET_PARAMS_NLN(p0_k1, eps_infty_k1, Nsol_k1, T_k1)

 IMPLICIT NONE

 REAL(q) :: p0_k1, eps_infty_k1, Nsol_k1, kb_k1, T_k1

 p0_k = p0_k1 
 eps_infty = eps_infty_k1
 Nsol_k = Nsol_k1
 T_k = T_k1

END SUBROUTINE  SET_PARAMS_NLN




!>========================== FUNCTION SHAPE_FUNC =================
!!
!! shape function
!!
!!=======================================================================
FUNCTION SHAPE_FUNC(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: SHAPE_FUNC, PS_TOTTN_I
 REAL(q), EXTERNAL :: ERRFC !compimentary error function implemented in errf.F in vasp.5.lib

 PS_TOTN_I = R_PS_TOTN(I)

 SHAPE_FUNC = ERRFC( log(MAX(REAL(PS_TOTN_I, KIND=q), CHGMIN)/nc_k) / SQRT(2._q)/sigma_k ) / 2._q

 RETURN

END FUNCTION SHAPE_FUNC




!>========================== FUNCTION LANG_NLN =================
!!
!! Compute the langevin function
!! input: nondimensional electric field, NDIM_E = \frac{p0 * |E|}{k * T}
!!
!!=======================================================================
FUNCTION LANG_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, LANG_NLN

 NDIM_E = NDIM_EFIELD(I)

!the calling function/subroutine must ensure that NDIM_E is not zer
 LANG_NLN = COTH(NDIM_E)/NDIM_E  - 1._q / NDIM_E

 RETURN

END FUNCTION LANG_NLN




!>========================== FUNCTION CHI_NLN =================
!!
!! Compute the electric susceptibiltiy
!! input: PS_TOTN = valence + pseudo core charge density in real space
!! input: nondimensional electric field, NDIM_E = \frac{p0 * |E|}{k * T}
!! output: CHI_NLN
!!
!!=======================================================================
FUNCTION CHI_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, CHI_NLN

 NDIM_E = NDIM_EFIELD(I)
 
!if the electric fiedl is non zero
 IF (ABS(NDIM_E) > CHGMIN) THEN

      CHI_NLN = SHAPE_FUNC(I) * ( 3 * (EB_K - EINFTY_K) * LANG_NLN(I)/NDIM_E + (EINFTY - 1._q) )

!electric field --> 0 , L(x)/x --> 1/3
 ELSE

      CHI_NLN = SHAPE_FUNC(I) * ( EB_K - 1._q )

 ENDIF

 RETURN

END FUNCTION CHI_NLN




!>========================== FUNCTION RHOB_NLN_G =================
!!
!! Compute the bound charge in fourier space
!! rho_b = -\grad \codt \grad \phi - \rho_s, where \phi is net potential
!! input: G2, square of G reciprocal lattice vector
!! input: PHI_G, Gth fourier component of PHI i.e PHI(G)
!! input: TOTN_G, Gth fourier component of TOTN i.e TOTN(G)
!! output: RHOB_NLN_G
!! used to compute the 0.5 * \int \phi * \rho_b d^3r energy term
!!
!!=======================================================================
FUNCTION RHOB_NLN_G(G, PHI_G, TOTN_G)

 IMPLICIT NONE

 REAL(q) :: G(3)
 REAL(q) :: G2, PHI_G, TOTN_G, RHOB_NLN_G

 G2 = G(1)**2 + G(2)**2 + G(3)**2 

!dont forgert the epsilon_0, permittivity of free space factor 
 RHOB_NLN_G = G2 * PHI_G - TOTN_G

 RETURN

END SUBROUTINE COMPUTE_RHOB_NLN




!>==========================SUBROUTINE UPDATE_EFIELD =======================
!!
!! Sets the module private variable NDIM_EFIELD(in real sapce), given the potential
!! input: PHI, potential in reciprocal space
!! NDIM_EFIELD = |\grad \phi| *p0 / k/ T; non dimensional electricfield in real space
!!
!!===========================================================================
SUBROUTINE UPDATE_EFIELD(GRIDC, LATT_CUR, PHI)


 IMPLICIT NONE

 TYPE(grid_3d),INTENT(IN)      :: GRIDC
 TYPE(latt),INTENT(IN)         :: LATT_CUR
 COMPLEX(q),INTENT(IN)         :: PHI(D_GRIDC%MPLWV)
 REAL(q), ALLOCATABLE           :: EFIELD(:)

 INTEGER(q)                    :: I,N1,N2,N3,NC
 REAL(q)                       :: GX,GY,GZ,tmp
 
!reurns EFIELD = |\grad \phi|^2
 CALL NABLAN(PHI, GRIDC, LATT_CUR, EFIELD)

 DO I = 1, GRIDC%RL%NP
   NDIM_EFIELD(I) = SQRT(EFIELD(I)) * p0_k / kb_k /T_k
 ENDDO

 RETURN

END SUBROUTINE UPDATE_EFIELD




!>==========================SUBROUTINE COMPUTE_EPSILON_NLN=================
!!
!!  
!!
!!==========================================================================

SUBROUTINE COMPUTE_EPSILON_NLN(EPSILON_NLN)

 IMPLICIT NONE

 REAL(q), INTENT(OUT) ::  EPSILON_NLN(CE_GRIDC%RL%NP)

 INTEGER(q) :: I
 REAL(q) :: CE_a, CE_b, CE_x, CE_n, CE_LG, DCE_LG


 EPSILON_NLN(I) = 1._q + CHI_NLN(I)

 RETURN

END SUBROUTINE COMPUTE_EPSILON_NLN




!>==========================SUBROUTINE MINIMIZE_NLN =================
!!
!!  nonlinear minimizer
!!
!!===================================================================
SUBROUTINE MINIMIZE_NLN(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, PS_TOTN, TOTN, SP_d_eps, phi_tilda, NCGSOL)

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: SP_T_INFO
 TYPE(potcar),INTENT(IN)    :: SP_P (SP_T_INFO%NTYP)
 TYPE (latt),INTENT(IN)     :: SP_LATT_CUR
 TYPE(grid_3d),INTENT(IN)   :: SP_GRIDC

 COMPLEX(q),INTENT(IN)     ::TOTN(SP_GRIDC%MPLWV)
 COMPLEX(q)                ::PS_TOTN(SP_GRIDC%MPLWV)
 REAL(q), INTENT(INOUT)    ::  SP_d_eps(SP_GRIDC%MPLWV)
 COMPLEX(q), INTENT(INOUT) ::  Vcav(SP_GRIDC%MPLWV), phi_tilda(SP_GRIDC%MPLWV)

 INTEGER(q) :: NCGSOL
 REAL(q)    :: Acav, diffphi

 diffphi = 1000._q

 MIN_counter = 1


 CALL UPDATE_EFIELD()
 
 CALL COMPUTE_EPSILON_NLN()

 DO WHILE ( (sqrt(diffphi) > 1e-6) )

#ifdef debugsol
     IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) "PK_counter, MIN_counter", PK_counter, MIN_counter
#endif

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
     CALL MINIMIZE(SP_GRIDC, SP_LATT_CUR, SP_P, SP_T_INFO, SP_d_eps, TOTN, phi_tilda, NCGSOL, diffphi)

     CALL UPDATE_EFIELD()
 
     CALL COMPUTE_EPSILON_NLN()

 ENDDO

 RETURN

END SUBROUTINE MINIMIZE_NLN


!>==========================SUBROUTINE VACORR_ELECTROSTATIC =================
!!
!!  computes the electroststic potential and energy corrections in reciprocal space
!!
!!===============================================================================
SUBROUTINE VACORR_ELECTROSTATIC(Vcorr, Acorr)

 IMPLICIT NONE

 REAL(q) :: Ael, Arhob, Acav, Adm, Acorr

 CALL EPS_POT(SP_GRIDC,SP_LATT_CUR,SP_T_INFO,SP_P,RWORK,phi_tilda,SP_d_eps,Vrhob)

 Vcorr = Vel + Vrhob
 Acorr = Ael + Arhob

Arhob = sum_g rhob_g * phi_g

END SUBROUTINE VACORR_ELECTROSTATIC




!>==========================SUBROUTINE VACORR_NON_ELECTROSTATIC =================
!!
!!  computes the non electroststic potential and energy correctionsin real space
!!
!!===============================================================================
SUBROUTINE VACORR_NON_ELECTROSTATIC(Vcorr, Acorr)

 IMPLICIT NONE

!Adm , energy of the dielectric medium
 REAL(q) :: Ael, Arhob, Acav, Adm, Acorr


 CALL CAVITY_POT_A(Vcav, Acav)
 CALL DIELMEDIUM_POT_A(Vdm, Adm)

 Vcorr = Vcav + Vdm
 Acorr = Acav + Adm - Vcorr*n


END SUBROUTINE VACORR_NON_ELECTROSTATIC




END MODULE POT_NLN_K
