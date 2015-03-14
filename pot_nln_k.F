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
 !private subroutines and functions
 PRIVATE :: SET_PARAMS_NLN, COMPUTE_EPSILON_NLN, SET_GVEC_NLN
 PRIVATE :: VA_RHOB, VA_CAV, VA_DMED
 PRIVATE :: VACORR_ELECTROSTATIC, VACORR_NON_ELECTROSTATIC
 PRIVATE :: INIT_EFIELD, UPDATE_EFIELD, SET_EEFF_TABLE, MINIMIZE_NLN
 !for test and debug
 PRIVATE :: TEST_K, WRITE_ARRAY_TO_FILE
 !functions
 PRIVATE :: SHAPE_FUNC, SHAPE_FUNC_DER, LANG_NLN, LANG_X_NLN, CHI_NLN
 PRIVATE :: EF, F_E, G_E !e, f(e), g(e)

 TYPE(type_info), PRIVATE :: T_INFO_NLN
 TYPE(potcar), PRIVATE, ALLOCATABLE :: P_NLN(:)
 TYPE (latt), PRIVATE :: LATT_CUR_NLN
 TYPE(grid_3d), PRIVATE :: GRIDC_NLN
 TYPE (in_struct), PRIVATE :: IO_NLN

 LOGICAL, SAVE :: LBCHG = .FALSE.

 !solvation paramters
 REAL(q), PRIVATE, SAVE :: p0_k, einfty_k, Nsol_k, kb_k, T_k
 REAL(q), PRIVATE, SAVE :: KT, XX_k, dv, alpha_k
 REAL(q), PRIVATE, SAVE :: pkT !p0/k/T

 REAL(q), PRIVATE, PARAMETER :: HART_TO_EV_K = 2 * RYTOEV

 REAL(q), PRIVATE, SAVE :: EDIFFSOL_INNER, EDIFFSOL_OUTER, DIFFPHI=1._q

 !soln to the poisson eqn
 COMPLEX(q), PRIVATE, ALLOCATABLE :: PHI_TILDA(:) 
 !modulus of electri field vector in real space
 REAL(q), PRIVATE, ALLOCATABLE :: EFIELD(:)
 !ps_totn in real space
 REAL(q), PRIVATE, ALLOCATABLE :: R_PS_TOTN(:) 
 !array of G vectors
 REAL(q), PRIVATE, ALLOCATABLE :: GVEC_NLN(:,:) 

 INTEGER(q), PRIVATE :: SIZE_REAL, SIZE_REC !allocation size
 INTEGER(q), PRIVATE :: SIZE_RLNP, SIZE_RCNP !actual size

 !effective field table size
 INTEGER(q), PRIVATE, PARAMETER :: N_EF = 9000
 !effective field table
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: E_EFF_TAB(:)
 !spline coefficients for the effective field
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: E_EFF_SPL(:,:)
 !effective field corresponding to each 
 !actual field (obtained from spline fit)
 REAL(q), PRIVATE, ALLOCATABLE :: E_EFF(:)
 REAL(q), PRIVATE, SAVE :: MIN_COUNTER = 1

 CONTAINS

!>========================== SUBROUTINE VCORRECTION_NLN_PCM ===========
!!
!! inputs TOTN and PS_TOTN are in reciprocal space and has the volume scaling 
!!returns VEL, just the poiss diff, and VCORR, the rest of the pot, in reciprocal space
!!note: in pot_k.F VDIEL refers to VCORR and VCORR in pot_k.F refers to VEL here
!!
!!=====================================================================
SUBROUTINE VCORRECTION_NLN_PCM(TOTN, PS_TOTN, VCORR, VEL, ATOT )

 IMPLICIT NONE

 COMPLEX(q), INTENT(OUT) ::  VCORR(SIZE_REC)
 COMPLEX(q),INTENT(OUT) :: VEL(SIZE_REC)
 COMPLEX(q),INTENT(IN) :: TOTN(SIZE_REC), PS_TOTN(SIZE_REC) 
 COMPLEX(q), ALLOCATABLE :: PS_TOTN_TMP(:)

 INTEGER(q) :: NCGSOL
 REAL(q)  :: ADIEL, ACAV, ATOT

 ALLOCATE(PS_TOTN_TMP(SIZE_REC))

 VEL   = 0._q
 VCORR = 0._q
 ADIEL  = 0._q
 ACAV  = 0._q
 ATOT  = 0._q

!!!set R_PS_TOTN, ps-totn in real space with volume scaling removed

 PS_TOTN_TMP = PS_TOTN / LATT_CUR_NLN%OMEGA

 CALL FFT3D(PS_TOTN_TMP, GRIDC_NLN, 1)

 R_PS_TOTN = 0._q

!!!r_ps_totn in the units of 1/Angstrom^3
 CALL UNPACK_C2R(PS_TOTN_TMP, 1.0_q, PS_TOTN_TMP, 0.0_q, R_PS_TOTN, GRIDC_NLN)

!!!end of setting R_PS_TOTN

 !nonlinear solver
 CALL MINIMIZE_NLN(TOTN, NCGSOL) 

 CALL UPDATE_EFIELD()

!!!returns pot and energy in atomic units
!!! ADIEL = ADM + AEPS + AEL
!!! where AEL = 0.5 \int (\phi - \phi_H) \rho_solute
!!! and VEL is the potential corresponding to AEL
!!! VCORR = VDM + VRHOB(pot corresponding to AEPS)
 CALL VACORR_ELECTROSTATIC(TOTN, ADIEL, VEL, VCORR)

!!!returns pot and energy in atomic units
!!!cavitation energy and potential
!!! VCORR = VCORR + VCAV
 CALL VACORR_NON_ELECTROSTATIC(PS_TOTN, ACAV, VCORR)

!!! net dielectric energy = electrostatic + cavitation
 ATOT = ADIEL + ACAV

#ifdef debugsol
 IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ATOT(in a.u) = &
      & ADM + AEPS + AEL + ACAV  =  ", ATOT
#endif

!!!to eV
 VEL   = VEL * HART_TO_EV_K
 VCORR = VCORR * HART_TO_EV_K
 ATOT  = ATOT * HART_TO_EV_K

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ATOT(in eV) = ", ATOT
#endif

!!!compute and write the boundcharge to file, filename: BCHG
 IF (LBCHG) CALL COMPUTE_RHOB_NLN()

 DEALLOCATE(PS_TOTN_TMP)

 RETURN

END SUBROUTINE VCORRECTION_NLN_PCM


!>========================== SUBROUTINE CONSTRUCTOR_NLN ===============
!!
!! The constructor
!!
!!=====================================================================
SUBROUTINE CONSTRUCTOR_NLN(IO, GRIDC, LATT_CUR, P, T_INFO, sigma_nc_k1, sigma_k1, nc_k1, &
                           & eb_k1, CORE_C1, sigma_rc_k1, tau1, ediffsol1, &
                           & p0_k1, einfty_k1, Nsol_k1, T_k1, LBCHG_1)

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: T_INFO
 TYPE(potcar),INTENT(IN) :: P(T_INFO%NTYP)
 TYPE (latt),INTENT(IN) :: LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: GRIDC
 TYPE (in_struct),INTENT(IN) :: IO

 LOGICAL :: LBCHG_1
 REAL(q) :: sigma_nc_k1, sigma_k1, nc_k1
 REAL(q) :: eb_k1, CORE_C1, sigma_rc_k1(T_INFO%NTYP), tau1, ediffsol1
 REAL(q) :: p0_k1, einfty_k1, Nsol_k1, T_k1

 LBCHG = LBCHG_1
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
 ALLOCATE(EFIELD(SIZE_REAL), R_PS_TOTN(SIZE_REAL))
 ALLOCATE(PHI_TILDA(SIZE_REC), GVEC_NLN(3, SIZE_RCNP))
 ALLOCATE(E_EFF(SIZE_REAL))

 IF (PK_COUNTER == 1)  THEN
    ALLOCATE(E_EFF_TAB(N_EF))
    ALLOCATE(E_EFF_SPL(N_EF, 5))
    E_EFF_TAB = 0._q
    E_EFF_SPL = 0._q

    ALLOCATE(SOL_PHI(SIZE_REC))
    SOL_PHI = 0._q

    CALL SET_PARAMS_LPCM(sigma_nc_k1, sigma_k1, nc_k1, eb_k1, CORE_C1, sigma_rc_k1, tau1, &
                      & ediffsol1, T_INFO_NLN%NTYP)
    CALL SET_PARAMS_NLN(p0_k1, einfty_k1, Nsol_k1, T_k1)
 ENDIF

 P_NLN = P
 PHI_TILDA = 0._q
! EFIELD = 0._q
 R_PS_TOTN = 0._q
 GVEC_NLN = 0._q
 E_EFF = 0._q

 CALL SET_GVEC_NLN()

 RETURN

END SUBROUTINE CONSTRUCTOR_NLN


!>========================== SUBROUTINE DESTRUCTOR_NLN ================
!!
!! The destructor
!!
!!=====================================================================
SUBROUTINE DESTRUCTOR_NLN()

 IMPLICIT NONE

 DEALLOCATE(P_NLN, EFIELD, R_PS_TOTN, PHI_TILDA, GVEC_NLN, E_EFF)

 PK_COUNTER = PK_COUNTER + 1

 RETURN

END SUBROUTINE DESTRUCTOR_NLN


!>========================== SUBROUTINE SET_PARAMS_NLN ================
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
!INCAR input in mol/liter
! to convert to molecules per Angstrom^3 multiply by 6.022141 * 1E-4
! eg:- molar concentration of water = 55.5 mol/L
! ==> 0.033 molecules/Angstrom^3
! also multiply by AUTOA**3 to convert to atomic units
!for the above example the value is approximately 0.005 in atomic units
!-----------------------------------------------------------------
 Nsol_k   = Nsol_k1 * 6.022 * 1E-4_q * (AUTOA**3)

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
!
!p0_k in atomic units
!-----------------------------------------------------------------
! p0_k =  SQRT( (3._q * (eb_k - einfty_k) * KT)/ (4._q * PI * Nsol_k )  )
 p0_k     =  p0_k1 * 0.393430_q

!-----------------------------------------------------------------
!multiply by pkT to make the electric field nondimensional
!in atomic units
!note: electric field must be in atomic units too
!-----------------------------------------------------------------
 pkT      = p0_k/KT

!-----------------------------------------------------------------
!the alpha parameter from the paper
!-----------------------------------------------------------------
 alpha_k =  3._q - (4._q * PI * Nsol_k * p0_k**2) / KT / (eb_k - einfty_k)

!-----------------------------------------------------------------
!the X parameter from the paper
!-----------------------------------------------------------------
 XX_k =  KT * (einfty_k - 1._q) / ( 4._q * PI * Nsol_k * p0_k**2 ) 


#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "eb_k, einfty_k, Nsol_k, p0_k, T_k, KT, pkT, alpha_k, XX_k (all in a.u) ", eb_k, einfty_k, Nsol_k, p0_k, T_k, KT, pkT, alpha_k, XX_k
#endif

 RETURN

END SUBROUTINE  SET_PARAMS_NLN


!>========================== FUNCTION SHAPE_FUNC ======================
!!
!! shape function in real space
!!
!!=====================================================================
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


!>========================== FUNCTION SHAPE_FUNC_DER ==================
!!
!! shape function derivative in real space
!!
!!=====================================================================
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


!>========================== FUNCTION LANG_NLN ========================
!!
!! Compute the langevin function in real space, coth(x) - 1/x
!! NDIM_E = x
!!
!!=====================================================================
FUNCTION LANG_NLN(NDIM_E)

 IMPLICIT NONE

 REAL(q) :: NDIM_E, LANG_NLN

 LANG_NLN = ( COSH(NDIM_E)/ SINH(NDIM_E) )  - ( 1._q / NDIM_E)

!langevin function cannot be greater than 1
 IF (LANG_NLN > 1._q) THEN 
    WRITE(*,*) "LANG_NLN > 1._q", LANG_NLN;
    STOP
 ENDIF

 RETURN

END FUNCTION LANG_NLN


!>========================== FUNCTION LANG_X_NLN ==========================
!!
!! Compute f(e) = ( coth(e) - 1/e) /e, in real space
!! e = effective field s.t e(1-alpha*f(e)) = x
!! x = \frac{p0 * |E|}{k * T}
!!
!! input, I = index
!! output, LANG_X_NLN = f(e)
!!
!!=======================================================================
FUNCTION LANG_X_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, LANG_X_NLN, SINHX, COSHX, EXPX, EXP_X

 NDIM_E = E_EFF(I) !EFIELD(I) * pkT

 LANG_X_NLN = F_E(NDIM_E)

 RETURN

END FUNCTION LANG_X_NLN


!>========================== FUNCTION F_E =============================
!!
!! Compute f(e) = ( coth(e) - 1/e) /e, in real space
!! e = effective field s.t e(1-alpha*f(e)) = x
!! x = \frac{p0 * |E|}{k * T}
!!
!! input, NDIM_E = e
!! output, F_E = f(e)
!!
!!=====================================================================
FUNCTION F_E(NDIM_E)

 IMPLICIT NONE

 REAL(q) :: NDIM_E, F_E, SINHX, COSHX, EXPX, EXP_X

!only interested in the magnitude of the electric field
 IF (NDIM_E < 0) THEN
    WRITE(*,*) "NDIM_E < 0 ", NDIM_E
    STOP
 ENDIF

 !coth(x) and 1/x, individually diverges as x-->0
 ! but coth(x) - 1/x --> 0 as x-->0
 ! when x is very small, numerical errors dominate. 
 ! L(x)/x --> 0 as x --> \infty since L(x) --> 1 as x --> \infty
 ! L(x)/x --> 1/3 as x--> 0, 
 !note: L(x)/x ~= 1/3 - 1/45 x^2  ==> for x = 1e-4 , L(x)/x - (1/3) ~= 2e-10
 !on the other hand ( (cosh(x)/sinh(x) - 1/x)/x - (1/3) ) ~= 5e-9

 !case 1: x -->0, L/x --> 1/3
 IF (NDIM_E < 1E-4_q) THEN
    F_E = 1._q/3._q - 1./45._q * (NDIM_E**2)

    !case 2: x --> \infty, L/x --> 1/x since L --> 1
    !but exp(x) = 1e30 corresponds to x ~= 70, if the L/x is approximated as 1/x then the error
    ! introduced is pretty significant. But, with L/x ~= 1/x - 1/x^2, the error is negligible for x >= 70
 ELSE IF (EXP(NDIM_E) > 1E30_q) THEN
        F_E = 1._q/NDIM_E - 1._q/(NDIM_E**2)

 ELSE
     F_E = LANG_NLN(NDIM_E) / NDIM_E

 ENDIF

 !better safe than sorry
 !0 < L(x)/x < 1/3
 IF ( ( (F_E - 1._q/3._q) > 1E-6) .OR. (F_E < 0) ) THEN 
    WRITE(*,*) "ERROR: 0 < L(x)/x < 1/3 ", F_E
    STOP
 ENDIF

 RETURN

END FUNCTION F_E


!>========================== FUNCTION CHI_NLN =========================
!!
!! Compute the electric susceptibiltiy in real space
!! \chi(e) = \epsilon(e) - 1 (dimensionless)
!! where e is the efeective field
!! mind: the value returned is not modulated by the shape function
!!
!!=====================================================================
FUNCTION CHI_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, CHI_NLN

 NDIM_E = E_EFF(I) !EFIELD(I) * pkT

 ! NDIM_E = e
 ! e --> 0 , \chi --> eb_k -1 , the linear pcm \chi
 ! E_EFF is set to 0 at intialization
 !==> \chi = eb_k -1 for the first iteration
 CHI_NLN =  (EB_K - EINFTY_K) * ( (3._q - alpha_k)/alpha_k ) * ( G_E(NDIM_E) - 1._q )  + (EINFTY_K - 1._q) 
 
 RETURN

END FUNCTION CHI_NLN


!>========================== SUBROUTINE RHOB_NLN ======================
!!
!! compute the bound charge in real space
!! bound charge, rho_b = -\nabla \cdot P 
!! = - \epsilon_0 \nabla \chi \vec{E}
!! =  \epsilon_0 \nabla \cdot (\chi (\nabla \phi) )
!!=====================================================================
SUBROUTINE COMPUTE_RHOB_NLN()

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q), ALLOCATABLE :: R_CHI_WORK(:), R_RHOB_WORK(:)
 COMPLEX(q), ALLOCATABLE :: C_RHOB_WORK(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_x_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_y_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_z_in(:)

 ALLOCATE(R_CHI_WORK(SIZE_REAL), R_RHOB_WORK(SIZE_RLNP), C_RHOB_WORK(SIZE_REC))
 ALLOCATE(gradphi_x_in(SIZE_REC), gradphi_y_in(SIZE_REC))
 ALLOCATE(gradphi_z_in(SIZE_REC))

!!!initialize
 R_CHI_WORK = 0._q
 R_RHOB_WORK = 0._q

!susceptibility in real space
 DO I = 1, SIZE_RLNP
    R_CHI_WORK(I) = CHI_NLN(I) * SHAPE_FUNC(I) 
 ENDDO

!!!c_rhob_work = div(chi*grad phi) in reciprocal space
 CALL LEPS(GRIDC_NLN, LATT_CUR_NLN, T_INFO_NLN, PHI_TILDA, R_CHI_WORK, &
                & gradphi_x_in, gradphi_y_in, gradphi_z_in, C_RHOB_WORK)

!!!note: in the above subroutine call, phi_tilda used was scaled by edeps/volume
!!!in a.u, without the volume scaling
!!!C_RHOB_WORK = C_RHOB_WORK * (4 * PI * AUTOA / EDEPS) * (AUTOA**2) * 4 * PI
!!!==> in units of Angstrom^-3
C_RHOB_WORK = C_RHOB_WORK * (4 * PI / EDEPS) * 4 * PI

!to real space
 CALL FFT3D(C_RHOB_WORK, GRIDC_NLN, 1)

!!!R_CHI_WORK overwritten
 CALL UNPACK_C2R(C_RHOB_WORK, 1.0_q, C_RHOB_WORK, 0.0_q, R_CHI_WORK, GRIDC_NLN)

!!!multiply by the volume(in ang^3) to be consistent with CHG(or CHGCAR) file format
 DO I = 1, SIZE_RLNP
    R_RHOB_WORK(I) = R_CHI_WORK(I) * LATT_CUR_NLN%OMEGA
 ENDDO

  CALL WRITE_TO_FILE(GRIDC_NLN, LATT_CUR_NLN, T_INFO_NLN, 'BCHG', R_RHOB_WORK)

 DEALLOCATE(R_CHI_WORK, R_RHOB_WORK, C_RHOB_WORK)
 DEALLOCATE(gradphi_x_in,gradphi_y_in,gradphi_z_in)

 RETURN

END SUBROUTINE COMPUTE_RHOB_NLN


!>==========================SUBROUTINE INIT_EFIELD ====================
!!
!! initializes EFIELD = |E| = |\grad \phi|, using Sol_phi from pot_lpcm_cav.F
!! Sol_phi doesnt have edeps/omega scaling 
!! EFIELD in real space and in a.u
!! also sets the lookup table for the effective field as well as 
!! the spline coefficents table
!!
!!=====================================================================
SUBROUTINE INIT_EFIELD()

 IMPLICIT NONE

 INTEGER(q)   :: I
 REAL(q)      :: GX,GY,GZ,tmp, E_MAX, E_MIN
 REAL(q)      :: X, Y, YDER
 REAL(q), ALLOCATABLE :: E_X(:)

 !reurns EFIELD = |\grad \phi|^2 in real space
 CALL NABLAN(SOL_PHI, GRIDC_NLN, LATT_CUR_NLN, EFIELD)

 ! |E| = |\grad \phi| in real space and in atomic units
 DO I = 1, SIZE_RLNP
    !we only need the magnitude of the electric field
    !in atomic units
    EFIELD(I) = SQRT(ABS(EFIELD(I))) * 4 * PI * AUTOA * AUTOA / LATT_CUR_NLN%OMEGA
    !    EFIELD_TMP(I) = EFIELD(I)
 ENDDO

 !max and min of electricfield in atomic units
 E_MIN = 1E-6_q
 E_MAX = 100._q

 !set the lookup table and the spline coefficients
 !sets E_X, E_EFF_TAB, E_EFF_SPL
 IF (PK_COUNTER == 1) THEN 
    ALLOCATE(E_X(N_EF))
#ifdef debugsol
     IF (IO_NLN%IU0>=0) WRITE(IO_NLN%IU0,*) "Setting up the lookup table and the spline coeffs"
#endif
    CALL SET_EEFF_TABLE(E_MIN, E_MAX, E_X)
#ifdef debugsol
    CALL WRITE_ARRAY_TO_FILE(E_X, N_EF, GRIDC_NLN, 'E_X')
    CALL WRITE_ARRAY_TO_FILE(E_EFF_TAB, N_EF, GRIDC_NLN, 'EF_TABLE')
#endif
    CALL SPLCPY(E_X, E_EFF_TAB, E_EFF_SPL, N_EF, N_EF, 10E30_q) 
    DEALLOCATE(E_X)
 ENDIF

 !set the effective field array, E_EFF, corresponding to EFIELD 
 !using spline fit
 DO I = 1, SIZE_RLNP
    X = EFIELD(I) * pkT
    CALL SPLVAL(X, Y, YDER, E_EFF_SPL, N_EF, N_EF)
    E_EFF(I) =  Y
 ENDDO

#ifdef debugsol
  CALL WRITE_TO_FILE(GRIDC_NLN, LATT_CUR_NLN, T_INFO_NLN, 'E_EFF', E_EFF)
  CALL WRITE_TO_FILE(GRIDC_NLN, LATT_CUR_NLN, T_INFO_NLN, 'EFIELD', EFIELD)
#endif

 RETURN

END SUBROUTINE INIT_EFIELD


!>==========================SUBROUTINE UPDATE_EFIELD ==================
!!
!! updates EFIELD = |E| = |\grad \phi|
!! EFIELD in real space and in a.u
!! also updates the corresponding effective field, E_EFF
!!
!!=====================================================================
SUBROUTINE UPDATE_EFIELD()

 IMPLICIT NONE

 INTEGER(q) :: I,N1,N2,N3,NC
 REAL(q) :: GX,GY,GZ,tmp
 REAL(q) :: X, Y, YDER
 
 !returns EFIELD = |\grad \phi|^2 in real space
 CALL NABLAN(PHI_TILDA, GRIDC_NLN, LATT_CUR_NLN, EFIELD)

 ! |E| = |\grad \phi| in real space and in atomic units
 DO I = 1, SIZE_RLNP
    !we only need the magnitude of the electric field
    !in atomic units
    EFIELD(I) = SQRT(ABS(EFIELD(I))) * 4 * PI * AUTOA * AUTOA / EDEPS
 ENDDO

 !spline fit at the nondimensional electric field
 !set the effective field array, E_EFF
 DO I = 1, SIZE_RLNP
    X = EFIELD(I) * pkT
    CALL SPLVAL(X, Y, YDER, E_EFF_SPL, N_EF, N_EF)
    E_EFF(I) =  Y
 ENDDO

 RETURN

END SUBROUTINE UPDATE_EFIELD


!>==========================SUBROUTINE SET_EEFF_TABLE =================
!!
!! set up the effective field table on a log grid of the electric field, E that ranges from E_MIN to E_MAX
!! input, E_MIN and E_MAX in a.u 
!! also sets E_X, dimensionless actual electric field = E * pkT 
!! E_X is used later to set the spline coefficients table
!!
!!=====================================================================
SUBROUTINE SET_EEFF_TABLE(E_MIN, E_MAX, E_X)

 IMPLICIT NONE

 INTEGER(q)   :: I
 REAL(q)      :: E, E_MIN, E_MAX, DE, E_X(N_EF), DIMLESS_E

 DO I = 1, N_EF
    E = E_MIN * (E_MAX/E_MIN)**(REAL(I,q)/REAL(N_EF,q))
    DIMLESS_E = E * pkT !(E_MIN + DE * (I-1)) * pkT
    E_X(I) = DIMLESS_E
    E_EFF_TAB(I) = EF(DIMLESS_E, 1E-8_q)
 ENDDO

 RETURN

END SUBROUTINE SET_EEFF_TABLE


!>==========================FUNCTION EF================================
!!
!! computes the effective electric field, 
!! e = (pE/kT) / (1 - alpha * f(e))
!! input, DIMLESS_E = pE/kT 
!! ouput, EF = e
!!
!!=====================================================================
FUNCTION EF(DIMLESS_E, TOL)

  IMPLICIT NONE

  REAL(q) :: EF, EF0, D_EF, DIMLESS_E, TOL
  
  EF0 = DIMLESS_E * 3._q/(3._q - alpha_k)
  D_EF = EF0

  DO WHILE (ABS(D_EF) > TOL)
     EF = DIMLESS_E * G_E(EF0)
     D_EF = EF - EF0
     EF0 = EF
  END DO

END FUNCTION EF

!>========================== FUNCTION G_E =============================
!!
!!  g(e) = 1 / (1 - alpha*f(e))
!!  input, E = e, the effective electric field
!!  ouput, G_E = g(e)
!!
!!=====================================================================
FUNCTION G_E(E)

  IMPLICIT NONE

  REAL(q) :: E, G_E

  G_E = 1._q /(1._q - alpha_k * F_E(E))

END FUNCTION G_E


!>==========================SUBROUTINE COMPUTE_EPSILON_NLN=================
!!
!!  \espilon = 1 + \chi * s(n(r)), in real space
!!
!!==========================================================================
SUBROUTINE COMPUTE_EPSILON_NLN(EPSILON_NLN)

 IMPLICIT NONE

 REAL(q), INTENT(OUT) ::  EPSILON_NLN(SIZE_REAL)

 INTEGER(q) :: I

 EPSILON_NLN = 0._q

 DO I = 1, SIZE_RLNP
    EPSILON_NLN(I) = 1._q + CHI_NLN(I) * SHAPE_FUNC(I) 
 ENDDO

 RETURN

END SUBROUTINE COMPUTE_EPSILON_NLN


!>==========================SUBROUTINE MINIMIZE_NLN ===================
!!
!!  nonlinear minimizer
!!
!!=====================================================================
SUBROUTINE MINIMIZE_NLN(TOTN, NCGSOL)

 IMPLICIT NONE

 COMPLEX(q),INTENT(IN)   :: TOTN(SIZE_REC)
 REAL(q), ALLOCATABLE    :: D_EPS(:)

 INTEGER(q) :: NCGSOL 
 REAL(q) :: ADIEL, ADIEL0, DIFF_ADIEL

 ALLOCATE(D_EPS(SIZE_REAL))

 ADIEL = 0._q

 MIN_COUNTER = 1

 !the efields initialized from the previously computed phi, sol_phi
 ! for the scf iteration i.e whne pk_counter = 1, sol_phi = 0
 ! ==> efield = 0
 ! also updates the module variable, E_EFF, the effective field
 CALL INIT_EFIELD()
 !sets the dielectric constant using the module variable, E_EFF, the 
 !effective field computed from the actual field EFIELD
 CALL COMPUTE_EPSILON_NLN(D_EPS)

 !-------------------------------------------------------------
 !Solves the generalized poisson eqn using CG, returns phi_tilda 
 !scaled by edeps/cell vol
 !-------------------------------------------------------------
 CALL MINIMIZE(GRIDC_NLN, LATT_CUR_NLN, P_NLN, T_INFO_NLN, D_EPS, TOTN, PHI_TILDA, NCGSOL, DIFFPHI)
 ! get the net electrostatic dielectric energy contribution i.e
 ! skip the caviation energy term and vasp original vacuum contribution
 CALL VACORR_ELECTROSTATIC(TOTN, ADIEL0)

 DIFF_ADIEL = ADIEL0

 ! DO WHILE ( (SQRT(DIFFPHI) > 1e-6) )
 DO WHILE ( (ABS(DIFF_ADIEL) > 1E-6) )
    ! use the module variable PHI_TILDA to update EFIELD
    !also update E_EFF
    CALL UPDATE_EFIELD()
    !update the dielctric constant using the updated E_EFF
    CALL COMPUTE_EPSILON_NLN(D_EPS)
    !solve the poisson eqn using CG for the given dielctric constant
    ! returns phi_tilda 
    !scaled by edeps/cell vol
    CALL MINIMIZE(GRIDC_NLN, LATT_CUR_NLN, P_NLN, T_INFO_NLN, D_EPS,&
         & TOTN, PHI_TILDA, NCGSOL, DIFFPHI)
    !get the electrostatic contrib to the dielectric energy
    CALL VACORR_ELECTROSTATIC(TOTN, ADIEL)
    DIFF_ADIEL = ADIEL - ADIEL0
    ADIEL0 = ADIEL
#ifdef debugsol
    IF (IO_NLN%IU0>=0)  WRITE(IO_NLN%IU0,*) "PK_COUNTER, MIN_COUNTER, &
         & DIFFPHI, DIFF_ADIEL", PK_COUNTER, MIN_COUNTER, DIFFPHI, DIFF_ADIEL
#endif
    MIN_COUNTER = MIN_COUNTER + 1
 ENDDO

 DEALLOCATE(D_EPS)

 RETURN

END SUBROUTINE MINIMIZE_NLN


!>==========================SUBROUTINE VACORR_ELECTROSTATIC ===========
!!
!! computes the net electroststic potential and energy contributions
!! to the dielctric energy ie.e everything except Vcav and  Acav
!!
!!=====================================================================
SUBROUTINE VACORR_ELECTROSTATIC(TOTN, ADIEL, VEL, VCORR)

 IMPLICIT NONE

 REAL(q) :: AEL, ARHOB, ACAV, ADM, AEPS, Ecorr1, Ecorr2, ADIEL
 COMPLEX(q) :: TOTN(SIZE_REC)
 COMPLEX(q), OPTIONAL :: VEL(SIZE_REC), VCORR(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: PHI(:), VRHOB(:)
 COMPLEX(q), ALLOCATABLE :: VDM(:)

 ALLOCATE(PHI(SIZE_REC) )

 IF (PRESENT(VCORR)) THEN
    ALLOCATE(VRHOB(SIZE_REC))
    CALL VA_RHOB(TOTN, VRHOB)!!!, ARHOB)
    VCORR = VCORR + VRHOB
 ENDIF

 CALL SimplePoisson(GRIDC_NLN, LATT_CUR_NLN, TOTN, PHI, 1)

 !!!in atomic units
 IF (PRESENT(VEL)) VEL = (PHI_TILDA - PHI) * 4 * PI * AUTOA / EDEPS

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

 !delta_phi correction in atomic units
 ADIEL = AEL

 !rest of the electrostatic correction
 ADM = 0._q
 AEPS = 0._q

 !!!returns pot and energy in atomic units
 IF (PRESENT(VCORR)) THEN
    ALLOCATE( VDM(SIZE_REC) ) 
    VDM = 0._q
    CALL VA_DMED(ADM, AEPS, VDM)
    VCORR = VCORR + VDM
 ELSE
    CALL VA_DMED(ADM, AEPS)
 ENDIF

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ADM(in a.u) = ", ADM
#endif

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "AEPS(in a.u) = ", AEPS
#endif

!net electrostatic correction
 ADIEL = ADIEL + ADM + AEPS

 DEALLOCATE(PHI)

 IF (PRESENT(VCORR))  DEALLOCATE(VRHOB, VDM)

 RETURN

END SUBROUTINE VACORR_ELECTROSTATIC


!>==========================SUBROUTINE VACORR_NON_ELECTROSTATIC =======
!!
!!  computes the non electroststic potential and energy correctionsin recip space
!!  only cavitation
!!
!!=====================================================================
SUBROUTINE VACORR_NON_ELECTROSTATIC(PS_TOTN, ACAV, VCORR)

 IMPLICIT NONE

 REAL(q) :: ACAV, ECORR1
 COMPLEX(q) :: CHTOT
 COMPLEX(q) :: VCORR(SIZE_REC), PS_TOTN(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: VCAV(:)

 ALLOCATE(VCAV(SIZE_REC))

 VCAV = 0._q

 ACAV = 0._q

!!!returns pot and energy in atomic units
 CALL VA_CAV(PS_TOTN, VCAV, ACAV)

#ifdef debugsol
     IF (IO_NLN%IU0>=0)   WRITE(IO_NLN%IU0,*) "ACAV(in a.u) = ", ACAV
#endif

 VCORR = VCORR + VCAV

 DEALLOCATE(VCAV)

 RETURN

END SUBROUTINE VACORR_NON_ELECTROSTATIC


!>========================== SUBROUTINE VA_RHOB =======================
!!
!!  computes the potential due to the bound charge
!! i.e the potential corresponding to 
!! ARHOB = 0.5 \int \phi \rho_b = 0.5 \int \chi |E|^2 
!! VRHOB = -\frac{\partial \epsilon}{\partial n} |E|^2
!! VRHOB = -\frac{\partial \chi}{\partial n} |E|^2
!! since \chi  = \epsilon - 1
!! returns the potential in reciprocal space
!!
!!====================================================================
SUBROUTINE VA_RHOB(TOTN, VRHOB)!!!, ARHOB)

 IMPLICIT NONE

 INTEGER(q) :: I
 COMPLEX(q) :: TOTN(SIZE_REC)
 REAL(q) :: EPSILON_DER_I, ARHOB, NDIM_E
 RGRID, ALLOCATABLE :: VRHOB_WORK(:)
 COMPLEX(q) :: VRHOB(SIZE_REC)

 ALLOCATE(VRHOB_WORK(SIZE_REAL))

 VRHOB_WORK = 0._q

 DO I = 1, SIZE_RLNP
    NDIM_E = E_EFF(I) !EFIELD(I) * pkT
    EPSILON_DER_I = SHAPE_FUNC_DER(I) * CHI_NLN(I) 
    !\epsilon_0 d\epsilon/dn in atomic units
    EPSILON_DER_I = EPSILON_DER_I / (4 * PI * (AUTOA**3) )
    VRHOB_WORK(I) = - EPSILON_DER_I * EFIELD(I) * EFIELD(I)
 ENDDO 

 CALL PACK_R2C(VRHOB_WORK, VRHOB, GRIDC_NLN)

 !to fourier space
 CALL FFT3D(VRHOB, GRIDC_NLN, -1)
 CALL SETUNB(VRHOB, GRIDC_NLN)
 VRHOB = VRHOB/ GRIDC_NLN%NPLWV

 DEALLOCATE(VRHOB_WORK)

 RETURN

END SUBROUTINE VA_RHOB


!>========================== SUBROUTINE VA_CAV ========================
!!
!!  computes the potential and energy due to the cavity formation
!!
!!=====================================================================
SUBROUTINE VA_CAV(PS_TOTN, VCAV, ACAV)

 IMPLICIT NONE

 REAL(q),ALLOCATABLE        :: VWORK(:) ! Vcav in real space
 REAL(q),ALLOCATABLE        :: GN(:,:), GN2_TMP(:) !\nabla n and |\naba n|^2
 COMPLEX(q),ALLOCATABLE     :: GN2(:) ! |\naba n|^2
 REAL(q),ALLOCATABLE        :: GGN(:,:) !\nabla |\nabla n|
 REAL(q),ALLOCATABLE        :: GGN2(:) !| \nabla |\nabla n| |^2
 REAL(q),ALLOCATABLE        :: DIVN(:) ! laplacian of n
 REAL(q),ALLOCATABLE        :: TERM1(:) ! grad shape func wrt n
 RGRID, ALLOCATABLE       :: RWORK(:)

 REAL(q), EXTERNAL      :: ERRFC 
 INTEGER(q)             :: I, J, NPC, NPR, RLNP, RCNP
 REAL(q)                :: QS,dv,tmp,tmp1,tmp2
 REAL(q) :: ACAV
 COMPLEX(q) :: VCAV(SIZE_REC), PS_TOTN(GRIDC_NLN%MPLWV)

 NPC =  SIZE_REC !!GRIDC%MPLWV
 RLNP = SIZE_RLNP !!GRIDC%RL%NP
 RCNP = SIZE_RCNP !!GRIDC%RC%NP
 NPR =  SIZE_REAL !!DIMREAL(GRIDC%MPLWV) !NPR = NPC, in complex mode

 ALLOCATE(GN2(NPC)) !complex
 ALLOCATE(GN2_TMP(NPR), GN(3, RLNP) )
 ALLOCATE(GGN2(NPR), GGN(3, RLNP) )
 ALLOCATE(DIVN(NPR), TERM1(RLNP))
 ALLOCATE(RWORK(NPR), VWORK(NPR))

 VCAV = 0._q
 ACAV = 0._q
 VWORK = 0._q

 !volume element
 dv = LATT_CUR_NLN%OMEGA / GRIDC_NLN%NPLWV

 !remove the scaling in the reciprocal charge density
 PS_TOTN = PS_TOTN / LATT_CUR_NLN%OMEGA

 !charge density in rec space, returns \nabla n and |\nabla n|^2 in real space
 CALL NABLAN(PS_TOTN, GRIDC_NLN, LATT_CUR_NLN, GN2_TMP, GN)

 !charge density in rec space, returns laplacian of n, Lap(n)
 !Lap(n) = (\nabla \dot \nabla n) in real space
 CALL LAPLN(PS_TOTN, GRIDC_NLN, LATT_CUR_NLN, DIVN)

 !add -Lap(n)/|\nabla n| to Vcav and copy \sqrt(|\nabla n|^2) to gn2_tmp
 DO I = 1, SIZE_REAL !DIMREAL(GRIDC%MPLWV)

   TMP = SQRT(MAX(REAL(GN2_TMP(I),KIND=q), CHGMIN))
   VWORK(I) = VWORK(I) - (DIVN(I)) / TMP
   GN2_TMP(I) = TMP

 ENDDO

 !PS_TOTN to the real space
 CALL FFT3D(PS_TOTN, GRIDC_NLN, 1)

 !unpack complex array to a real array
 !in complex mode, just copies PS_TOTN to RWORK
 CALL UNPACK_C2R(PS_TOTN, 1.0_q, PS_TOTN, 0.0_q, RWORK, GRIDC_NLN)

 !charge density must be in real space(unscaled), 
 !returns derivative of shape function wrt the charge density in the real space,
 ! term1 = gamma*A / n, where
 ! gamma * A = exp( -( log(n/n_c) )^2 / (2 sigma^2) ) / ( sigma * sqrt(2*pi) )
 CALL SHAPE_GRADN(GRIDC_NLN, T_INFO_NLN, P_NLN, LATT_CUR_NLN, RWORK, TERM1)

 !quantum surface area
 QS = 0._q

 !quantum surface area, integral of (gamma*A / n) * |\nabla n| 
 != term1 * \sqrt(gn2_tmp)
 DO I = 1, SIZE_RLNP !GRIDC%RL%NP
   QS = QS + (TERM1(I)) * (GN2_TMP(I)) 
 ENDDO

 QS = QS * dv 

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, QS, 1))

 !cavitation energy
 !tau is in eV/Angstrom^3 and QS has the units of Angstrom^2
 !==> acav is in eV
 ACAV = tau * QS 

#ifdef debugsol
 IF (SP_IO%IU0>=0) THEN
   print *,'tau (eV/Angstrom^2)', tau
   print *,'surface area (Angstrom^2)', QS
   print *,'Acav (eV)', ACAV
 ENDIF
#endif

 !to a.u, to be consistent with the rest of the ATOT components
 ACAV = ACAV/ HART_TO_EV_K

 !1/ |nabla n|
 DO I = 1, SIZE_RLNP !DIMREAL(GRIDC%MPLWV)
   GN2_TMP(I) = 1._q / ( MAX(GN2_TMP(I), CHGMIN) ) 
 ENDDO

 !-------------------------------------------------------------
 !Compute cavitation potential, Vcav
 !-------------------------------------------------------------
 !packs the real array into a complex one
 !in complex mode, same as  gn2 = CMPLX(gn2_tmp)
 CALL PACK_R2C(GN2_TMP, GN2, GRIDC_NLN)

 !1/ |nabla n| to reciprocal space to take the gradient
 CALL FFT3D(GN2, GRIDC_NLN, -1)
 CALL SETUNB(GN2, GRIDC_NLN)
 GN2 = GN2 / GRIDC_NLN%NPLWV

 !\nabla(1 / |\nabla n|), ggn in real space
 CALL NABLAN(GN2, GRIDC_NLN, LATT_CUR_NLN, GGN2, GGN)

 !add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
 !and multiply by term1 = gamma*A/n in real space
 DO I = 1, SIZE_RLNP !GRIDC%RL%NP
  TMP = ( GN(1,I) * GGN(1,I) &
         & + GN(2,I) * GGN(2,I) &
         & + GN(3,I) * GGN(3,I) ) * ( TERM1(I) )
  VWORK(I) = VWORK(I) - TMP
 ENDDO

 CALL PACK_R2C(VWORK, VCAV, GRIDC_NLN)

 !!!Vcav to reciprocal space
 CALL FFT3D(VCAV, GRIDC_NLN, -1)
 CALL SETUNB(VCAV, GRIDC_NLN)
 VCAV = VCAV / GRIDC_NLN%NPLWV

 !!!ps_totn to rec space
 CALL FFT3D(PS_TOTN, GRIDC_NLN, -1)
 CALL SETUNB(PS_TOTN, GRIDC_NLN)
 PS_TOTN = PS_TOTN / GRIDC_NLN%NPLWV * LATT_CUR_NLN%OMEGA

 !!!multiply by the surface tension
 !!!in units of eV/Angstrom
 VCAV = VCAV * tau

 !to a.u to be consistent
 VCAV = VCAV * AUTOA / HART_TO_EV_K

 DEALLOCATE(RWORK, VWORK)
 DEALLOCATE(GN2, GN2_TMP)
 DEALLOCATE(GN, GGN2)
 DEALLOCATE(GGN, DIVN, TERM1)

 RETURN

END SUBROUTINE VA_CAV


!>========================== SUBROUTINE VA_DMED =======================
!!
!! computes  A1 = A_dm(e) + A_eps
!! where A_eps = - 0.5 * \int \chi |E|^2 
!! and the potential corresponding to A_dm
!! the potential corresponding to A_eps is computed in VA_RHOB
!! note1: VDM is optional
!! note2: net A_diel = A1 + 0.5 \int (\phi - \phi_H) \rho_solute 
!!
!!=====================================================================
SUBROUTINE VA_DMED(ADM, AEPS, VDM)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, dv, ADM, ADM_1
 REAL(q) :: AEPS, AEPS_1, FEFF
 REAL(q) :: EXPX, EXP_X, LOG_TERM
 COMPLEX(q), OPTIONAL :: VDM(SIZE_REC)
 RGRID, ALLOCATABLE :: VDM_WORK(:)

 IF (PRESENT(VDM)) ALLOCATE(VDM_WORK(SIZE_REAL))

!volume element
 dv =  LATT_CUR_NLN%OMEGA / GRIDC_NLN%NPLWV

!to atomic units
 dv = dv / (AUTOA**3)

 IF (PRESENT(VDM)) VDM_WORK = 0._q

 ADM = 0._q
 AEPS = 0._q

 DO I = 1, SIZE_RLNP
    NDIM_E = E_EFF(I) !EFIELD(I) * pkT
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

    !old
    !ADM_1 = ( (einfty_k - 1._q) / 6._q / (eb_k - einfty_k) + LANG_X_NLN(I) ) * (NDIM_E**2)  - LOG_TERM

    !f(e), e = effective field
    FEFF = LANG_X_NLN(I)

    ADM_1 = ( FEFF - 0.5_q*alpha_k* FEFF**2 + 0.5_q * XX_k * (1._q - alpha_k* FEFF)**2 ) * (NDIM_E**2)
    ADM_1 = ADM_1 - LOG_TERM

    ADM_1 = SHAPE_FUNC(I) * Nsol_k * KT *  ADM_1

    ADM = ADM + ADM_1 * dv

    !-0.5 \epsilon_0 \int \chi |E|^2
    AEPS = AEPS - 0.5_q * SHAPE_FUNC(I) * CHI_NLN(I) * EFIELD(I) * EFIELD(I) / 4 / PI * dv

    IF (PRESENT(VDM)) VDM_WORK(I) = SHAPE_FUNC_DER(I) * ADM_1

 ENDDO

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, ADM, 1))
 CALLMPI( M_sum_d(GRIDC_NLN%COMM, AEPS, 1))

#ifdef debugsol
     IF (IO_NLN%IU0>=0)  WRITE(IO_NLN%IU0,*) 'AEPS, ADM, ADM - AEPS (all in a.u) ', AEPS, ADM, ADM - AEPS
#endif

 IF (PRESENT(VDM)) THEN
    CALL PACK_R2C(VDM_WORK, VDM, GRIDC_NLN)
    CALL FFT3D(VDM, GRIDC_NLN, -1)
    CALL SETUNB(VDM, GRIDC_NLN)
    VDM = VDM/ GRIDC_NLN%NPLWV
    DEALLOCATE(VDM_WORK)
 ENDIF

 RETURN

END SUBROUTINE VA_DMED


!>========================== SUBROUTINE SET_GVEC_NLN ==================
!!
!! set the the reciprocal lattice vectors
!!
!!=====================================================================
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


!>========================== SUBROUTINE WRITE_ARRAY_TO_FILE ===========
!!
!! write a real array, ARRAY, of size N to file with name FNAME
!! GRIDC is used only to get the ionode 
!!
!!=====================================================================
SUBROUTINE WRITE_ARRAY_TO_FILE(ARRAY, N, GRIDC, FNAME)

 USE prec
 USE mgrid

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN) :: GRIDC
 INTEGER(q) :: N, I
 REAL(q) :: ARRAY(N)
 CHARACTER (*) FNAME

 CHARACTER (40) FORM, SZNAM
 INTEGER :: ISTAT, NODE_ME, IONODE, IU
 INTEGER :: NALLOC, NZ, NWRITE, NWRITTEN

 NODE_ME = 0
 IONODE = 0
 IU = 469
 SZNAM = 'untitled'
 FORM = '(1(1X,E17.11))'
 NWRITE = 5

#ifdef MPI
 NODE_ME = GRIDC%COMM%NODE_ME
 IONODE = GRIDC%COMM%IONODE
#endif

 io_begin
 OPEN(IU, FILE = FNAME, STATUS = 'UNKNOWN')
 REWIND IU
 io_end

 IF (ISTAT>0) RETURN ! can not write immediate exit

 do_io WRITE(IU,'(I5)') N

 NWRITTEN = 0

 io_begin
 DO I = 1, N
    NWRITTEN = NWRITTEN + 1
    IF ( MOD(NWRITTEN,NWRITE) == 0 ) THEN
       WRITE(IU, FORM) ARRAY(I)
    ELSE
       WRITE(IU, FORM, ADVANCE='NO') ARRAY(I)
    ENDIF
 ENDDO
 io_end

 IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

 RETURN

END SUBROUTINE WRITE_ARRAY_TO_FILE


!>========================== SUBROUTINE TEST_K ========================
!!
!! for testing L(x)/x limits
!!
!!=====================================================================
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
