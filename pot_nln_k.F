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

 PUBLIC  :: CONSTRUCTOR_NLN, VCORRECTION_NLN_PCM, DESTRUCTOR_NLN
 PRIVATE :: SHAPE_FUNC, SHAPE_FUNC_DER, LANG_NLN, LANG_X_NLN, CHI_NLN, RHOB_NLN_G
 PRIVATE :: SET_PARAMS_NLN, UPDATE_EFIELD, MINIMIZE_NLN, COMPUTE_EPSILON_NLN, SET_GVEC_NLN
 PRIVATE :: VACORR_ELECTROSTATIC, VACORR_NON_ELECTROSTATIC, VA_RHOB, VA_CAV, VA_DMED
 PRIVATE :: INIT_EFIELD, TEST_K

 TYPE(type_info), PRIVATE, SAVE :: T_INFO_NLN
 TYPE(potcar), PRIVATE, ALLOCATABLE, SAVE :: P_NLN(:)
 TYPE (latt), PRIVATE, SAVE :: LATT_CUR_NLN
 TYPE(grid_3d), PRIVATE, SAVE :: GRIDC_NLN
 TYPE (in_struct), PRIVATE, SAVE :: IO_NLN

 REAL(q), PRIVATE, SAVE    :: p0_k, einfty_k, Nsol_k, kb_k, T_k, KT, dv
 REAL(q), PRIVATE, SAVE    :: pkT !p0/k/T
 REAL(q), PRIVATE, PARAMETER    :: HART_TO_EV_K = 2 * RYTOEV
 REAL(q), PRIVATE, SAVE    :: EDIFFSOL_INNER, EDIFFSOL_OUTER, DIFFPHI=1._q
 COMPLEX(q), PRIVATE, ALLOCATABLE, SAVE  :: PHI_TILDA(:) !soln to the poisson eqn
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: EFIELD(:)!, EFIELD_TMP(:) !modulus of electri field vector in real space
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: R_PS_TOTN(:) !ps_totn in real space
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: GVEC_NLN(:,:) !array og G vectors

 INTEGER(q), PRIVATE, SAVE :: SIZE_REAL, SIZE_REC !allocation size
 INTEGER(q), PRIVATE, SAVE :: SIZE_RLNP, SIZE_RCNP !actual size



 CONTAINS



!>========================== SUBROUTINE CONSTRUCTOR_NLN ======================
!!
!!
!!
!!===============================================================================
SUBROUTINE CONSTRUCTOR_NLN(IO, GRIDC, LATT_CUR, P, T_INFO, sigma_nc_k1, sigma_k1, nc_k1, &
                           & eb_k1, CORE_C1, sigma_rc_k1, tau1, ediffsol1, &
                           & p0_k1, einfty_k1, Nsol_k1, T_k1)

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: T_INFO
 TYPE(potcar),INTENT(IN) :: P(T_INFO%NTYP)
 TYPE (latt),INTENT(IN) :: LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: GRIDC
 TYPE (in_struct),INTENT(IN) :: IO

 REAL(q) :: sigma_nc_k1, sigma_k1, nc_k1
 REAL(q) :: eb_k1, CORE_C1, sigma_rc_k1(T_INFO%NTYP), tau1, ediffsol1
 REAL(q) :: p0_k1, einfty_k1, Nsol_k1, T_k1


 IO_NLN = IO
 SP_IO = IO
 T_INFO_NLN = T_INFO

 LATT_CUR_NLN = LATT_CUR
 GRIDC_NLN = GRIDC

 EDIFFSOL_INNER = 1E-6
 EDIFFSOL_OUTER = EDIFFSOL_INNER * 10

 SIZE_REC = GRIDC_NLN%MPLWV

 SIZE_REAL = DIMREAL(GRIDC_NLN%MPLWV)

 SIZE_RCNP = GRIDC_NLN%RC%NCOL * GRIDC_NLN%RC%NROW

 SIZE_RLNP = GRIDC_NLN%RL%NCOL * GRIDC_NLN%RL%NROW

 ALLOCATE(P_NLN(T_INFO_NLN%NTYP))
! ALLOCATE(EFIELD(SIZE_REAL), EFIELD_TMP(SIZE_REAL), R_PS_TOTN(SIZE_REAL))
 ALLOCATE(EFIELD(SIZE_REAL), R_PS_TOTN(SIZE_REAL))
 ALLOCATE(PHI_TILDA(SIZE_REC), GVEC_NLN(3, SIZE_RCNP))

 IF (PK_COUNTER == 1)  THEN
   ALLOCATE(SOL_PHI(SIZE_REC))
   SOL_PHI = 0._q
 ENDIF

 P_NLN = P
 PHI_TILDA = 0._q
! EFIELD = 0._q
! EFIELD_TMP = 0._q
 R_PS_TOTN = 0._q
 GVEC_NLN = 0._q

 CALL SET_PARAMS_LPCM(sigma_nc_k1, sigma_k1, nc_k1, eb_k1, CORE_C1, sigma_rc_k1, tau1, &
                      & ediffsol1, T_INFO_NLN%NTYP)
 CALL SET_PARAMS_NLN(p0_k1, einfty_k1, Nsol_k1, T_k1)

 CALL SET_GVEC_NLN()

 RETURN

END SUBROUTINE CONSTRUCTOR_NLN




!>========================== SUBROUTINE DESTRUCTOR_NLN ======================
!!
!!
!!
!!===============================================================================
SUBROUTINE DESTRUCTOR_NLN()

 IMPLICIT NONE

! DEALLOCATE(P_NLN, EFIELD, EFIELD_TMP, R_PS_TOTN, PHI_TILDA, GVEC_NLN)
 DEALLOCATE(P_NLN, EFIELD, R_PS_TOTN, PHI_TILDA, GVEC_NLN)

 PK_COUNTER = PK_COUNTER + 1

 RETURN

END SUBROUTINE DESTRUCTOR_NLN




!>========================== SUBROUTINE SET_PARAMS_NLN =================
!!
!! Sets the params specifi to the nonlinear model
!!
!!=====================================================================
SUBROUTINE  SET_PARAMS_NLN(p0_k1, einfty_k1, Nsol_k1, T_k1)

 IMPLICIT NONE

 REAL(q) :: p0_k1, einfty_k1, Nsol_k1, kb_k1, T_k1

!-----------------------------------------------------------------
!dimensionless quatity
!einfty_k = (refractive index)**2
!-----------------------------------------------------------------

 einfty_k = einfty_k1

!-----------------------------------------------------------------
!number of solvent molecules per unit volume
!INACR input in mol/liter
! to convert to molecules per Angstrom^3 multiply by 6.022141 * 1E-4
! eg:- molar concentration of water = 55.5 mol/L
! ==> 0.033 molecules/Angstrom^3
! also multiply by AUTOA**3 to convert to atomic units
!for the above example the value is approximately 0.005 in atomic units
!-----------------------------------------------------------------

 Nsol_k   = Nsol_k1 * 6.022141 * 1E-4_q * (AUTOA**3)

!-----------------------------------------------------------------
!temperature in kelvin
!-----------------------------------------------------------------

 T_k      = T_k1

!-----------------------------------------------------------------
! boltzmann const in eV/K from constant.inc
!-----------------------------------------------------------------

 kb_k     = BOLKEV 

!-----------------------------------------------------------------
!in atomic units, hartree_to_ev = 2 * rydberg_to_ev
!-----------------------------------------------------------------

 KT = kb_k * T_k/(2 * RYTOEV)


!Note:
!for the above example, N*k*T = 0.005 * 8.6173857E-5 * 300/ 2 / 13.61 
! = 4.75e-06 hartree/a0^3, where a0 = bohr radius 
! ~= 0.9 meV/Angstrom^3


!-----------------------------------------------------------------
!1 Debye = 0.393430 e*a0
! where the e is the charge of the electron and a0 the bohr radius.
!in atomic units e=1 and a0=1 so 1 Debye = 0.393430 au
! p0_k is not independent, it can be computed from eb_k and einfty_k
! the value of p0_k obtained from the INCAR is not used
! p0_k     = p0_k1
!
!in atomic units \epsilon_0 = 1/4/\pi
!
!p0_k in atomic units
!-----------------------------------------------------------------

 p0_k =  SQRT( (3._q * (eb_k - einfty_k) * KT)/ (4._q * PI * Nsol_k )  )

!-----------------------------------------------------------------
!multiply by pkT to make the electric field nondimensional
!in atomic units
!note: electric field must be in atomic units too
!-----------------------------------------------------------------

 pkT      = p0_k/KT

 RETURN

END SUBROUTINE  SET_PARAMS_NLN



!>========================== FUNCTION SHAPE_FUNC =================
!!
!! shape function in real space
!!
!!=======================================================================
FUNCTION SHAPE_FUNC(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: SHAPE_FUNC, PS_TOTN_I
 REAL(q), EXTERNAL :: ERRFC 

!R_PS_TOTN has the units of 1/Angstrom^3, so nc_k must be of the same unit
 PS_TOTN_I = R_PS_TOTN(I)

 SHAPE_FUNC = ERRFC( log(MAX(REAL(PS_TOTN_I, KIND=q), CHGMIN)/nc_k) / SQRT(2._q)/sigma_k ) / 2._q

 RETURN

END FUNCTION SHAPE_FUNC



!>========================== FUNCTION SHAPE_FUNC_DER =================
!!
!! shape function derivative in real space
!!
!!=======================================================================
FUNCTION SHAPE_FUNC_DER(I)

 IMPLICIT NONE

 INTEGER(q) :: I

 REAL(q) :: EPR_C, EPR_Z, PS_TOTN_I, SHAPE_FUNC_DER

 PS_TOTN_I = R_PS_TOTN(I)

 EPR_C = 1._q / sqrt(2*PI) / sigma_k

 EPR_Z = log(MAX(REAL(PS_TOTN_I,KIND=q),CHGMIN)/nc_k) / sqrt(2._q) / sigma_k

 SHAPE_FUNC_DER = EPR_C * (EXP(-(EPR_Z)**2)) / (MAX(REAL(PS_TOTN_I,KIND=q),CHGMIN))

 RETURN

END FUNCTION SHAPE_FUNC_DER



!>========================== FUNCTION LANG_NLN ==========================
!!
!! Compute the langevin function in real space, coth(x) - 1/x
!! x = NDIM_E = \frac{p0 * |E|}{k * T}
!!
!!=======================================================================
FUNCTION LANG_NLN(NDIM_E)

 IMPLICIT NONE

 REAL(q) :: NDIM_E, LANG_NLN

 LANG_NLN = ( COSH(NDIM_E)/ SINH(NDIM_E) )  - ( 1._q / NDIM_E )

!langevin function cannot be greater than 1
 IF (LANG_NLN > 1._q) THEN 
    WRITE(*,*) "LANG_NLN > 1._q", LANG_NLN;
    STOP
 ENDIF

 RETURN

END FUNCTION LANG_NLN



!>========================== FUNCTION LANG_X_NLN ==========================
!!
!! Compute ( coth(x) - 1/x) /x, in real space
!! x = NDIM_E = \frac{p0 * |E|}{k * T}
!!
!!=======================================================================
FUNCTION LANG_X_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, LANG_X_NLN, SINHX, COSHX, EXPX, EXP_X

 NDIM_E = EFIELD(I) * pkT

!only interested in the magnitude of the electric field
 IF (NDIM_E < 0) THEN
    WRITE(*,*) "NDIM_E < 0 ", NDIM_E
    STOP
 ENDIF

!coth(x) and 1/x, individually diverges as x-->0
! but coth(x) - 1/x --> 0 as x-->0
! when x is very samll, numerical errors dominate. 
! L(x)/x --> 0 as x --> \infty since L(x) --> 1 as x --> \infty
! L(x)/x --> 1/3 as x--> 0, 
!note: L(x)/x ~= 1/3 - 1/45 x^2  ==> for x = 1e-4 , L(x)/x - (1/3) ~= 2e-10
!on the other hand ( (cosh(x)/sinh(x) - 1/x)/x - (1/3) ) ~= 5e-9


!case 1: x -->0, L/x --> 1/3
 IF (NDIM_E < 1E-4_q) THEN

    LANG_X_NLN = 1._q/3._q - 1./45 * (NDIM_E**2)

!case 2: x --> \infty, L/x --> 1/x since L --> 1
!but exp(x) = 1e30 corresponds to x ~= 70, if the L/x is approximated as 1/x then the error
! introduced is pretty significant. But, with L/x ~= 1/x - 1/x^2, the error is negligible for x >= 70
 ELSE IF (EXP(NDIM_E) > 1E30_q) THEN

        LANG_X_NLN = 1._q/NDIM_E - 1._q/(NDIM_E**2)

 ELSE

     LANG_X_NLN = LANG_NLN(NDIM_E) / NDIM_E

 ENDIF

! WRITE(*,*) 'NDIM_E, LANG_X_NLN  ', NDIM_E, LANG_X_NLN

!better safe than sorry
!0 < L(x)/x < 1/3
 IF ( ( (LANG_X_NLN - 1._q/3._q) > 1E-6) .OR. (LANG_X_NLN < 0) ) THEN 
    WRITE(*,*) "ERROR: 0 < L(x)/x < 1/3 ", LANG_X_NLN
    STOP
 ENDIF

 RETURN

END FUNCTION LANG_X_NLN




!>========================== FUNCTION CHI_NLN ===========================
!!
!! Compute the electric susceptibiltiy in real space
!!
!!=======================================================================
FUNCTION CHI_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, CHI_NLN

 NDIM_E = EFIELD(I) * pkT

 CHI_NLN =  ( 3 * (EB_K - EINFTY_K) * LANG_X_NLN(I) + (EINFTY_K - 1._q) )

! WRITE(*,*) 'shape_func, chi', SHAPE_FUNC(I), CHI_NLN

 CHI_NLN = CHI_NLN * SHAPE_FUNC(I) 

 RETURN

END FUNCTION CHI_NLN



!>========================== FUNCTION RHOB_NLN_G =================
!!
!! Gth the bound charge component
!! bound charge, rho_b = \epsilon_0 \grad \cdot E - \rho_s
!!    = -\epsilon_0 \grad \codt \grad \phi - \rho_s, where \phi is net potential
!!
!!=======================================================================
FUNCTION RHOB_NLN_G(I, TOTN_G)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: G(3)
 REAL(q) :: G2
 COMPLEX(q) :: TOTN_G, TOTN_G_SCALED, PHI_G, RHOB_NLN_G

!in atomic units and remove the units of energy
 PHI_G = PHI_TILDA(I) * 4 * PI * AUTOA /  EDEPS

 G(:) = GVEC_NLN(:,I)

!WRITE(*,*) 'G(1), G(2), G(3) ', G(1), G(2), G(3)

 G2 = ( G(1)**2 + G(2)**2 + G(3)**2 )

!to atomic units
 G2 = G2 * AUTOA * AUTOA

!remove the volume and convert to atomic units
 TOTN_G_SCALED = TOTN_G * (AUTOA**3) / LATT_CUR_NLN%OMEGA

 RHOB_NLN_G = G2 * PHI_G/4/PI - TOTN_G_SCALED

 RETURN

END FUNCTION RHOB_NLN_G




!>==========================SUBROUTINE INIT_EFIELD =======================
!!
!! initialozes EFILED = |E| = |\grad \phi|, using Sol_phi from pot_lpcm_cav.F
!! Sol_phi doesnt have edeps/omega scaling 
!! EFIELD in real space
!!
!!===========================================================================
SUBROUTINE INIT_EFIELD()

 IMPLICIT NONE


 INTEGER(q)   :: I,N1,N2,N3,NC
 REAL(q)      :: GX,GY,GZ,tmp
 
!reurns EFIELD = |\grad \phi|^2 in real space
 CALL NABLAN(SOL_PHI, GRIDC_NLN, LATT_CUR_NLN, EFIELD)

! |E| = |\grad \phi| in real space and in atomic units
 DO I = 1, SIZE_RLNP

!we only need the magnitude of the electric field
!in atomic units
    EFIELD(I) = SQRT(ABS(EFIELD(I))) * 4 * PI * AUTOA * AUTOA / LATT_CUR_NLN%OMEGA
!    EFIELD_TMP(I) = EFIELD(I)
    EFIELD(I) = 0._q

 ENDDO

 RETURN

END SUBROUTINE INIT_EFIELD



!>==========================SUBROUTINE UPDATE_EFIELD =======================
!!
!! updates EFILED = |E| = |\grad \phi|
!! EFIELD in real space
!! input: PHI, potential in reciprocal space
!!
!!===========================================================================
SUBROUTINE UPDATE_EFIELD()

 IMPLICIT NONE

 INTEGER(q)   :: I,N1,N2,N3,NC
 REAL(q)      :: GX,GY,GZ,tmp
 
!reurns EFIELD = |\grad \phi|^2 in real space
 CALL NABLAN(PHI_TILDA, GRIDC_NLN, LATT_CUR_NLN, EFIELD)

! |E| = |\grad \phi| in real space and in atomic units
 DO I = 1, SIZE_RLNP

!we only need the magnitude of the electric field
!in atomic units
    EFIELD(I) = SQRT(ABS(EFIELD(I))) * 4 * PI * AUTOA * AUTOA / EDEPS
!    EFIELD_TMP(I) = EFIELD(I)

 ENDDO

 RETURN

END SUBROUTINE UPDATE_EFIELD



!>==========================SUBROUTINE COMPUTE_EPSILON_NLN=================
!!
!!  \espilon = 1 + \chi, in real space
!!
!!==========================================================================

SUBROUTINE COMPUTE_EPSILON_NLN(EPSILON_NLN)

 IMPLICIT NONE

 REAL(q), INTENT(OUT) ::  EPSILON_NLN(SIZE_REAL)

 INTEGER(q) :: I

 EPSILON_NLN = 0._q

 DO I = 1, SIZE_RLNP

    EPSILON_NLN(I) = 1._q + CHI_NLN(I)

!    WRITE(*,*) 'R_PS_TOTN(I), EFIELD(I), EPSILON_NLN(I) ', R_PS_TOTN(I), EFIELD(I), EPSILON_NLN(I)

 ENDDO

 RETURN

END SUBROUTINE COMPUTE_EPSILON_NLN



!>==========================SUBROUTINE MINIMIZE_NLN =================
!!
!!  nonlinear minimizer
!!
!!===================================================================
SUBROUTINE MINIMIZE_NLN(TOTN, NCGSOL)

 IMPLICIT NONE

 COMPLEX(q),INTENT(IN)   :: TOTN(SIZE_REC)
 REAL(q), ALLOCATABLE    :: D_EPS(:)

 INTEGER :: NCGSOL, MIN_COUNTER

 ALLOCATE(D_EPS(SIZE_REAL))

 MIN_COUNTER = 1

!uncomment the statemt if efields need to initialized from 
!the previously computed phi, sol_phi
 CALL INIT_EFIELD()
 
 CALL COMPUTE_EPSILON_NLN(D_EPS)

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
     CALL MINIMIZE(GRIDC_NLN, LATT_CUR_NLN, P_NLN, T_INFO_NLN, D_EPS, TOTN, PHI_TILDA, NCGSOL, DIFFPHI)

 DO WHILE ( (SQRT(DIFFPHI) > 1e-6) )

     CALL UPDATE_EFIELD()
 
     CALL COMPUTE_EPSILON_NLN(D_EPS)

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
     CALL MINIMIZE(GRIDC_NLN, LATT_CUR_NLN, P_NLN, T_INFO_NLN, D_EPS, TOTN, PHI_TILDA, NCGSOL, DIFFPHI)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "PK_COUNTER, MIN_COUNTER, DIFFPHI", PK_COUNTER, MIN_COUNTER, DIFFPHI
#endif

     MIN_COUNTER = MIN_COUNTER + 1

 ENDDO

 DEALLOCATE(D_EPS)

 RETURN

END SUBROUTINE MINIMIZE_NLN



!>==========================SUBROUTINE VACORR_ELECTROSTATIC =================
!!
!!  computes the electroststic potential and energy corrections in reciprocal space
!!
!!===============================================================================
SUBROUTINE VACORR_ELECTROSTATIC(TOTN, VEL, VCORR, ATOT)

 IMPLICIT NONE

 REAL(q) :: AEL, ARHOB, ACAV, ADM, ATOT, Ecorr1, Ecorr2
 COMPLEX(q) :: TOTN(SIZE_REC), VEL(SIZE_REC), VCORR(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: PHI(:), VRHOB(:)

 ALLOCATE(PHI(SIZE_REC), VRHOB(SIZE_REC))

 CALL VA_RHOB(TOTN, VRHOB, ARHOB)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ARHOB(in a.u) = ", ARHOB
#endif

 CALL SimplePoisson(GRIDC_NLN, LATT_CUR_NLN, TOTN, PHI, 1)

!!!in atomic units
 VEL = (PHI_TILDA - PHI) * 4 * PI * AUTOA / EDEPS

 VCORR = VCORR + VRHOB

!!!TOTN has the volume factor in it
 CALL MY_D_PROD(Ecorr1, TOTN, PHI_TILDA, GRIDC_NLN)

!!!in atomic units
 Ecorr1 = Ecorr1 * 0.5_q * 4 * PI * AUTOA / EDEPS

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, Ecorr1, 1))

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "Ecorr1(in a.u) = ", Ecorr1
#endif

 CALL MY_D_PROD(Ecorr2, TOTN, PHI, GRIDC_NLN)

!!!in atomic units
 Ecorr2 = Ecorr2 * 0.5_q * 4 * PI * AUTOA / EDEPS

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, Ecorr2, 1))

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "Ecorr2(in a.u) = ", Ecorr2
#endif

!!!in atomic units
 AEL = Ecorr1 - Ecorr2

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "AEL(in a.u) = ", AEL
#endif

!net correction in atomic units
 ATOT = ATOT + AEL !!!!!+ ARHOB

 DEALLOCATE(PHI, VRHOB)

 RETURN

END SUBROUTINE VACORR_ELECTROSTATIC



!>==========================SUBROUTINE VACORR_NON_ELECTROSTATIC =================
!!
!!  computes the non electroststic potential and energy correctionsin recip space
!!
!!===============================================================================
SUBROUTINE VACORR_NON_ELECTROSTATIC(VCORR, ATOT)

 IMPLICIT NONE

 REAL(q) :: ACAV, ADM, ECORR1, ATOT, AEPS
 COMPLEX(q) :: CHTOT
 COMPLEX(q) :: VCORR(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: VCAV(:), VDM(:)

 ALLOCATE(VCAV(SIZE_REC), VDM(SIZE_REC))

 VCAV = 0._q
 VDM = 0._q

 ACAV = 0._q
 ADM = 0._q
 AEPS = 0._q

!!!returns pot and energy in atomic units
 CALL VA_CAV(VCAV, ACAV)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ACAV(in a.u) = ", ACAV
#endif

!!!returns pot and energy in atomic units
 CALL VA_DMED(VDM, ADM, AEPS)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ADM(in a.u) = ", ADM
#endif

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "AEPS(in a.u) = ", AEPS
#endif

 VCORR = VCORR + VCAV + VDM

 ATOT = ATOT + ACAV + ADM + AEPS

 DEALLOCATE(VCAV, VDM)

 RETURN

END SUBROUTINE VACORR_NON_ELECTROSTATIC



!>========================== SUBROUTINE VA_RHOB =================
!!
!!  computes the potential and energy due to the bound charge
!! VRHOB = -\frac{\partial \epsilon}{\partial n} |E|^2
!! ARHOB = 0.5 \int \phi \rho_b
!! returns the potential in reciprocal space
!!
!!===============================================================
SUBROUTINE VA_RHOB(TOTN, VRHOB, ARHOB)

 IMPLICIT NONE

 INTEGER(q) :: I, FACT_K
 COMPLEX(q) :: TOTN(SIZE_REC)
 REAL(q) :: EPSILON_DER_I, ARHOB, NDIM_E
 RGRID, ALLOCATABLE :: VRHOB_WORK(:)
 COMPLEX(q) :: VRHOB(SIZE_REC)

 ALLOCATE(VRHOB_WORK(SIZE_REAL))

 ARHOB = 0._q
 VRHOB_WORK = 0._q

#ifdef realmode
 FACT_K = 2
#elif
 FACT_K = 1
#endif

 DO I = 1, SIZE_RLNP
 
!      NDIM_E = EFIELD_TMP(I) * pkT
      NDIM_E = EFIELD(I) * pkT

     EPSILON_DER_I = SHAPE_FUNC_DER(I) * ( 3 * (EB_K - EINFTY_K) * LANG_X_NLN(I) + (EINFTY_K - 1._q) )

!\epsilon_0 d\epsilon/dn in atomic units
     EPSILON_DER_I = EPSILON_DER_I / (4 * PI * (AUTOA**3) )

!     VRHOB_WORK(I) = - EPSILON_DER_I * EFIELD_TMP(I) * EFIELD_TMP(I)
     VRHOB_WORK(I) = - EPSILON_DER_I * EFIELD(I) * EFIELD(I)

 ENDDO 

 CALL PACK_R2C(VRHOB_WORK, VRHOB, GRIDC_NLN)

!to fourier space
 CALL FFT3D(VRHOB, GRIDC_NLN, -1)
 CALL SETUNB(VRHOB, GRIDC_NLN)
 VRHOB = VRHOB/ GRIDC_NLN%NPLWV

!iterate over the G vectors
 DO I = 1, SIZE_RCNP

    ARHOB = ARHOB + 0.5_q * FACT_K *CONJG( RHOB_NLN_G(I, TOTN(I)) ) &
                    & * ( PHI_TILDA(I) * 4 * PI * AUTOA / EDEPS ) &
                    & * LATT_CUR_NLN%OMEGA

 ENDDO

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, ARHOB, 1))

 DEALLOCATE(VRHOB_WORK)

 RETURN

END SUBROUTINE VA_RHOB



!>========================== SUBROUTINE VA_CAV =================
!!
!!  computes the potential and energy due to the cavity formation
!!
!!===============================================================
SUBROUTINE VA_CAV(VCAV, ACAV)

 IMPLICIT NONE

 REAL(q) :: ACAV
 COMPLEX(q) :: VCAV(SIZE_REC)

 VCAV = 0._q
 ACAV = 0._q

 RETURN

END SUBROUTINE VA_CAV



!>========================== SUBROUTINE VA_DMED =================
!!
!!  computes the potential and energy due to the dielectric medium
!!
!!===============================================================
SUBROUTINE VA_DMED(VDM, ADM, AEPS)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, dv, ADM, ADM_1
 REAL(q) :: AEPS, AEPS_1
 REAL(q) :: EXPX, EXP_X, LOG_TERM
 COMPLEX(q) :: VDM(SIZE_REC)
 RGRID, ALLOCATABLE :: VDM_WORK(:)

 ALLOCATE(VDM_WORK(SIZE_REAL))

!volume element
 dv =  LATT_CUR_NLN%OMEGA / GRIDC_NLN%NPLWV

!to atomic units
 dv = dv / (AUTOA**3)

 VDM_WORK = 0._q

 ADM = 0._q
 AEPS = 0._q

 DO I = 1, SIZE_RLNP

!    NDIM_E = EFIELD_TMP(I) * pkT
    NDIM_E = EFIELD(I) * pkT

    EXPX = EXP(NDIM_E)
    EXP_X = EXP(-NDIM_E)

!case 1: x --> 0
!log(sinh(x)/x) --> 0 as x --> 0
    IF ( (ABS(EXPX - 1._q) < 1E-6_q) .AND. (ABS(EXP_X - 1._q) < 1E-6_q) )  THEN
       LOG_TERM = 0._q

!case 2: x --> \infty, sinh(x) --> exp(x)/2 ==> log(sinh(x)/x) --> ( x - log(2x) )
!exp(70) ~= 1E30
    ELSEIF ( (EXPX > 1E30_q) .AND. (ABS(EXP_X) < 1E-6_q) )  THEN
       LOG_TERM = NDIM_E - LOG(2 * NDIM_E)

    ELSE

       LOG_TERM = LOG(SINH(NDIM_E) / NDIM_E)
    ENDIF

!    WRITE(*,*) 'NDIM_E, LOG_TERM, L/x', NDIM_E, LOG_TERM, LANG_X_NLN(I)

    ADM_1 = ( (einfty_k - 1._q) / 6._q / (eb_k - einfty_k) + LANG_X_NLN(I) ) * (NDIM_E**2)  - LOG_TERM

    ADM_1 = SHAPE_FUNC(I) * Nsol_k * KT *  ADM_1

    ADM = ADM + ADM_1 * dv

!-0.5 \epsilon_0 \int \chi |E|^2
    AEPS = AEPS - 0.5_q * CHI_NLN(I) * EFIELD(I) * EFIELD(I) / 4 / PI * dv

    VDM_WORK(I) = SHAPE_FUNC_DER(I) * ADM_1

 ENDDO

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, ADM, 1))
 CALLMPI( M_sum_d(GRIDC_NLN%COMM, AEPS, 1))

 WRITE(*,*) 'AEPS, ADM, ADM - AEPS', AEPS, ADM, ADM - AEPS

 CALL PACK_R2C(VDM_WORK, VDM, GRIDC_NLN)

 CALL FFT3D(VDM, GRIDC_NLN, -1)
 CALL SETUNB(VDM, GRIDC_NLN)
 VDM = VDM/ GRIDC_NLN%NPLWV
 
 DEALLOCATE(VDM_WORK)

 RETURN

END SUBROUTINE VA_DMED



!>========================== SUBROUTINE SET_GVEC_NLN ======================
!!
!! set the the reciprocal lattice vectors
!!
!!===============================================================================

SUBROUTINE SET_GVEC_NLN()

 IMPLICIT NONE

 INTEGER(q) :: I, NC, N1, N2, N3
 
 I = 0

 col: DO NC = 1, GRIDC_NLN%RC%NCOL
         N2 = GRIDC_NLN%RC%I2(NC)
         N3 = GRIDC_NLN%RC%I3(NC)
         row: DO N1 = 1, GRIDC_NLN%RC%NROW

               I = I +1
            
               GVEC_NLN(1, I) = ( GRIDC_NLN%LPCTX(N1) * LATT_CUR_NLN%B(1,1) &
                                & + GRIDC_NLN%LPCTY(N2) * LATT_CUR_NLN%B(1,2) &
                                & + GRIDC_NLN%LPCTZ(N3) * LATT_CUR_NLN%B(1,3) )

               GVEC_NLN(2, I) = ( GRIDC_NLN%LPCTX(N1) * LATT_CUR_NLN%B(2,1) &
                                 & + GRIDC_NLN%LPCTY(N2) * LATT_CUR_NLN%B(2,2) &
                                 & + GRIDC_NLN%LPCTZ(N3) * LATT_CUR_NLN%B(2,3) )

               GVEC_NLN(3, I) = ( GRIDC_NLN%LPCTX(N1) * LATT_CUR_NLN%B(3,1) &
                                 & + GRIDC_NLN%LPCTY(N2) * LATT_CUR_NLN%B(3,2) &
                                 & + GRIDC_NLN%LPCTZ(N3) * LATT_CUR_NLN%B(3,3) )
               
               GVEC_NLN(:,I) = GVEC_NLN(:,I) * TPI
 
         ENDDO row
      ENDDO col

 RETURN

END SUBROUTINE SET_GVEC_NLN



!>========================== SUBROUTINE VCORRECTION_NLN_PCM ======================
!!
!! inputs TOTN and PS_TOTN are in reciprocal space and has the volume scaling 
!!returns VEL, just the poiss diff, and VCORR, the rest of the pot, in reciprocal space
!!note: in pot_k.F VDIEL refers to VCORR and VCORR in pot_k.F refers to VEL here
!!
!!===============================================================================

SUBROUTINE VCORRECTION_NLN_PCM(TOTN, PS_TOTN, VCORR, VEL, ATOT )

 IMPLICIT NONE

 COMPLEX(q), INTENT(OUT) ::  VCORR(SIZE_REC)
 COMPLEX(q),INTENT(OUT) :: VEL(SIZE_REC)
 COMPLEX(q),INTENT(IN) :: TOTN(SIZE_REC), PS_TOTN(SIZE_REC) 
 COMPLEX(q), ALLOCATABLE :: PS_TOTN_TMP(:)

 INTEGER :: NCGSOL
 REAL(q)  :: ATOT

 ALLOCATE(PS_TOTN_TMP(SIZE_REC))

 VEL   = 0._q
 VCORR = 0._q
 ATOT  = 0._q

!!!set R_PS_TOTN, ps-totn in real space with volume scaling removed

 PS_TOTN_TMP = PS_TOTN / LATT_CUR_NLN%OMEGA

 CALL FFT3D(PS_TOTN_TMP, GRIDC_NLN, 1)

 R_PS_TOTN = 0._q

!!!r_ps_totn in the units of 1/Angstrom^3
 CALL UNPACK_C2R(PS_TOTN_TMP, 1.0_q, PS_TOTN_TMP, 0.0_q, R_PS_TOTN, GRIDC_NLN)

!!!end of setting R_PS_TOTN

!!! CALL TEST_K()
 
 CALL MINIMIZE_NLN(TOTN, NCGSOL) 

 CALL UPDATE_EFIELD()

!!!poisson diff and rhob contributions
!!!returns pot and energy in atomic units
!!!atot = atot + ael + arhob, arhob not used
!!!atot = atot + ael
 CALL VACORR_ELECTROSTATIC(TOTN, VEL, VCORR, ATOT)

!!!cav and dielectric medium
!!!atot = atot + acav + adm + aeps
 CALL VACORR_NON_ELECTROSTATIC(VCORR, ATOT)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ATOT(in a.u) = AEL + ACAV + ADM + AEPS  =  ", ATOT
#endif

!!!to eV
 VEL   = VEL * HART_TO_EV_K
 VCORR = VCORR * HART_TO_EV_K
 ATOT  = ATOT * HART_TO_EV_K

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ATOT(in eV) = ", ATOT
#endif

 DEALLOCATE(PS_TOTN_TMP)

 RETURN

END SUBROUTINE VCORRECTION_NLN_PCM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!for testing L(x)/x limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TEST_K()

 IMPLICIT NONE

 REAL(q) :: x, LX

 x = 1e6

 WRITE(*,*) 'X , X_NICE ', X, X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX

 x = 1e-3

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX - 1._q/3._q


 x = 1e-4

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX - 1._q/3._q

 x = 1e-5

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX - 1._q/3._q

 x = 1e-6

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX - 1._q/3._q

 x = 1e-7 

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x = ', LX - 1._q/3._q

 x = 1e-15

 WRITE(*,*) 'X, X_NICE ', X, X_NICE(X)
 LX = LANG_NLN(X_NICE(X)) / X_NICE(X)
! LX = LANG_NLN(X_NICE(X)) / X
 WRITE(*,*) 'L/x - 1/3 = ', LX - 1._q/3._q



 CONTAINS

   FUNCTION X_NICE(X1)
     IMPLICIT NONE
     REAL(q) :: X1,X_NICE
     X_NICE = MAX( X1, 1E-6_q )
     X_NICE = MIN( X_NICE, 700._q )
   END FUNCTION X_NICE


END SUBROUTINE TEST_K


END MODULE POT_NLN_K
