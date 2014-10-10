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
 PRIVATE :: SET_PARAMS_NLN, COMPUTE_EFIELD_NLN, COMPUTE_LANG_NLN 
 PRIVATE :: COMPUTE_CHI_NLN, COMPUTE_EPSILON_NLN, COMPUTE_RHOB_NLN, GET_V_E_NLN
 PRIVATE :: MINIMIZE_NLN

 TYPE(type_info), PRIVATE, SAVE :: T_INFO_NLN
 TYPE(potcar), PRIVATE, ALLOCATABLE, SAVE :: P_NLN(:)
 TYPE (latt), PRIVATE, SAVE :: LATT_CUR_NLN
 TYPE(grid_3d), PRIVATE, SAVE :: GRIDC_NLN
 TYPE (in_struct), PRIVATE, SAVE :: IO_NLN

 REAL(q), PRIVATE, SAVE               :: p0_k, einfty_k, Nsol_k, kb_k, T_k, dv
 REAL(q), PRIVATE, SAVE    :: pkT !p0/k/T
 COMPLEX(q), PRIVATE, ALLOCATABLE, SAVE  :: PHI_TILDA(:) !soln to the poisson eqn
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: EFIELD(:) !modulus of electri field vector in real space
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: R_PS_TOTN(:) !ps_totn in real space
 REAL(q), PRIVATE, ALLOCATABLE, SAVE  :: GVEC_NLN(:,:) !array og G vectors

 INTEGER(q), PRIVATE, SAVE :: SIZE_REAL, SIZE_REC !array size in real and reciprocal space



 CONTAINS



!>========================== SUBROUTINE CONSTRUCTOR_NLN ======================
!!
!!
!!
!!===============================================================================
SUBROUTINE CONSTRUCTOR_NLN(IO, GRIDC, LATT_CUR, P, T_INFO, sigma_nc_k, sigma_k, nc_k, &
                           & eb_k, CORE_C, sigma_rc_k, tau, ediffsol, &
                           & p0_k1, einfty_k1, Nsol_k1, T_k1)

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: T_INFO
 TYPE(potcar),INTENT(IN) :: P(T_INFO%NTYP)
 TYPE (latt),INTENT(IN) :: LATT_CUR
 TYPE(grid_3d),INTENT(IN) :: GRIDC
 TYPE (in_struct),INTENT(IN) :: IO

 REAL(q) :: sigma_nc_k, sigma_k, nc_k
 REAL(q) :: eb_k, CORE_C, sigma_rc_k(T_INFO%NTYP), tau, ediffsol
 REAL(q) :: p0_k1, einfty_k1, Nsol_k1, T_k1

 INTEGER(q) :: SIZE_REAL, SIZE_REC, SIZE_GVEC

 IO_NLN = IO
 T_INFO_NLN = T_INFO

 LATT_CUR_NLN = LATT_CUR
 GRIDC_NLN = GRIDC

 SIZE_REC = GRIDC_NLN%MPLWV

 SIZE_REAL = DIMREAL(GRIDC_NLN%MPLWV)

 SIZE_GVEC = GRIDC_NLN%RC%NCOL * GRIDC_NLN%RC%NROW

 ALLOCATE(P_NLN(T_INFO_NLN%NTYP))
 ALLOCATE(EFIELD(SIZE_REAL), R_PS_TOTN(SIZE_REAL))
 ALLOCATE(PHI_TILDA(SIZE_REC), GVEC_NLN(3, SIZE_GVEC))

 IF (PK_COUNTER == 1)  THEN
   ALLOCATE(SOL_PHI(SIZE_REC))
   SOL_PHI = 0._q
 ENDIF

 P_NLN = P

 CALL SET_PARAMS_LPCM(sigma_nc_k, sigma_k, nc_k, eb_k, CORE_C, sigma_rc_k, tau, &
                      & ediffsol, T_INFO_NLN%NTYP)
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

 DEALLOCATE(EFIELD, R_PS_TOTN, PHI_TILDA, GVEC_NLN)

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

 p0_k = p0_k1 
 einfty_k = einfty_k1
 Nsol_k = Nsol_k1
 T_k = T_k1
 kb_k = 1.0
 pkT = p0_k/(kb_k * T_k)

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
!! Compute the langevin function in real space
!! NDIM_E = \frac{p0 * |E|}{k * T}
!!
!!=======================================================================
FUNCTION LANG_NLN(I)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, LANG_NLN

 NDIM_E = EFIELD(I) * pkT

!the calling function/subroutine must ensure that NDIM_E is not zer
!coth(x) - 1/x
 LANG_NLN = (COSH(NDIM_E)/ SINH(NDIM_E))  - 1._q / NDIM_E

 RETURN

END FUNCTION LANG_NLN




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
 
!if the electric fiedl is non zero
 IF (ABS(NDIM_E) > CHGMIN) THEN

      CHI_NLN = SHAPE_FUNC(I) * ( 3 * (EB_K - EINFTY_K) * LANG_NLN(I)/NDIM_E + (EINFTY_K - 1._q) )

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
FUNCTION RHOB_NLN_G(I, TOTN_G)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: G(3)
 REAL(q) :: G2
 COMPLEX(q) :: TOTN_G,PHI_G, RHOB_NLN_G

 PHI_G = PHI_TILDA(I) * EDEPS / LATT_CUR_NLN%OMEGA
 G = GVEC_NLN(:,I)

 G2 = G(1)**2 + G(2)**2 + G(3)**2 

 RHOB_NLN_G = G2 * PHI_G - TOTN_G

 RETURN

END FUNCTION RHOB_NLN_G




!>==========================SUBROUTINE UPDATE_EFIELD =======================
!!
!! Sets the module private variable NDIM_EFIELD(in real sapce), given the potential
!! input: PHI, potential in reciprocal space
!! NDIM_EFIELD = |\grad \phi| *p0 / k/ T; non dimensional electricfield in real space
!!
!!===========================================================================
SUBROUTINE UPDATE_EFIELD()


 IMPLICIT NONE


 INTEGER(q)                    :: I,N1,N2,N3,NC
 REAL(q)                       :: GX,GY,GZ,tmp
 
!reurns EFIELD = |\grad \phi|^2 in real space
 CALL NABLAN(PHI_TILDA, GRIDC_NLN, LATT_CUR_NLN, EFIELD)

 DO I = 1, GRIDC_NLN%RL%NP
   EFIELD(I) = SQRT(EFIELD(I))
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

 REAL(q), INTENT(OUT) ::  EPSILON_NLN(SIZE_REAL)

 INTEGER(q) :: I

 DO I = 1, GRIDC_NLN%RL%NP
    EPSILON_NLN(I) = 1._q + CHI_NLN(I)
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

 COMPLEX(q),INTENT(IN)     ::TOTN(SIZE_REC)
 REAL(q), ALLOCATABLE    ::  SP_d_eps(:)

 INTEGER :: NCGSOL, MIN_COUNTER
 REAL(q)    :: diffphi

 ALLOCATE(SP_D_EPS(SIZE_REAL))

 diffphi = 1000._q

 MIN_counter = 1


 CALL UPDATE_EFIELD()
 
 CALL COMPUTE_EPSILON_NLN(SP_D_EPS)

 DO WHILE ( (sqrt(diffphi) > 1e-6) )

#ifdef debugsol
     IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) "PK_counter, MIN_counter", PK_counter, MIN_counter
#endif

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
     CALL MINIMIZE(GRIDC_NLN, LATT_CUR_NLN, P_NLN, T_INFO_NLN, SP_d_eps, TOTN, PHI_TILDA, NCGSOL, diffphi)

     CALL UPDATE_EFIELD()
 
     CALL COMPUTE_EPSILON_NLN(SP_D_EPS)

     MIN_COUNTER = MIN_COUNTER + 1

 ENDDO

 DEALLOCATE(SP_D_EPS)

 RETURN

END SUBROUTINE MINIMIZE_NLN




!>==========================SUBROUTINE VACORR_ELECTROSTATIC =================
!!
!!  computes the electroststic potential and energy corrections in reciprocal space
!!
!!===============================================================================
SUBROUTINE VACORR_ELECTROSTATIC(TOTN, VEL, VCORR, ATOT)

 IMPLICIT NONE

 REAL(q) :: Ael, Arhob, Acav, Adm, ATOT, Ecorr1, Ecorr2
 COMPLEX(q) :: TOTN(SIZE_REC), VEL(SIZE_REC), VCORR(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: PHI(:), VRHOB(:)

 ALLOCATE(PHI(SIZE_REC), VRHOB(SIZE_REC))

 CALL VA_RHOB(TOTN, VRHOB, ARHOB)

 CALL SimplePoisson(GRIDC_NLN, LATT_CUR_NLN, TOTN, PHI, 1)

 VEL = PHI_TILDA - PHI

 VCORR = VCORR + VRHOB

 CALL MY_D_PROD(Ecorr1, TOTN, PHI_TILDA, GRIDC_NLN)

 Ecorr1 = Ecorr1 * 0.5_q

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, Ecorr1, 1))

 CALL MY_D_PROD(Ecorr2, TOTN, PHI, GRIDC_NLN)

 Ecorr2 = Ecorr2 * 0.5_q

 CALLMPI( M_sum_d(GRIDC_NLN%COMM, Ecorr2, 1))

 AEL = Ecorr1 - Ecorr2

!net correction
 ATOT = ATOT + AEL + ARHOB

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

 REAL(q) :: Acav, Adm, ECORR1, ATOT
 COMPLEX(q) :: CHTOT
 COMPLEX(q) :: VCORR(SIZE_REC)
 COMPLEX(q), ALLOCATABLE :: VCAV(:), VDM(:)

 ALLOCATE(VCAV(SIZE_REC), VDM(SIZE_REC))

 VCAV = 0._q
 VDM = 0._q
 VCORR = 0._q

 ACAV = 0._q
 ADM = 0._q

 CALL VA_CAV(VCAV, ACAV)
 CALL VA_DMED(VDM, ADM)

 VCORR = VCORR + VCAV + VDM

 ATOT = ATOT + ACAV + ADM

 DEALLOCATE(VCAV, VDM)

 RETURN

END SUBROUTINE VACORR_NON_ELECTROSTATIC



!>========================== SUBROUTINE VA_RHOB =================
!!
!!  computes the potential and energy due to the bound charge
!! VRHOB = -\frac{\partial \epsilon}{\partial n} |E|^2
!! ARHOB = 
!! returns the potential in reciprocal space
!!
!!===============================================================
SUBROUTINE VA_RHOB(TOTN, VRHOB, ARHOB)

 IMPLICIT NONE

 INTEGER(q) :: I
 COMPLEX(q) :: TOTN(SIZE_REC)
 REAL(q) :: EPSILON_DER_I, ARHOB, NDIM_E
 REAL(q), ALLOCATABLE :: VRHOB_WORK(:)
 COMPLEX(q) :: VRHOB(SIZE_REC)

 ALLOCATE(VRHOB_WORK(SIZE_REAL))

 ARHOB = 0._q
 VRHOB_WORK = 0._q


 DO I = 1, SIZE_REAL
 
      NDIM_E = EFIELD(I) * pkT

     EPSILON_DER_I = SHAPE_FUNC_DER(I) * ( 3 * (EB_K - EINFTY_K) * LANG_NLN(I)/NDIM_E + (EINFTY_K - 1._q) )

     VRHOB_WORK(I) = - EPSILON_DER_I * EFIELD(I) * EFIELD(I)

 ENDDO 

 CALL PACK_R2C(VRHOB_WORK, VRHOB, GRIDC_NLN)

!to fourier space
 CALL FFT3D(VRHOB, GRIDC_NLN, -1)


!iterate over the G vectors
 DO I = 1, SIZE_REC

    ARHOB = ARHOB + 0.5 * RHOB_NLN_G(I, TOTN(I)) * CONJG(PHI_TILDA(I))

 ENDDO

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
SUBROUTINE VA_DMED(VDM, ADM)

 IMPLICIT NONE

 INTEGER(q) :: I
 REAL(q) :: NDIM_E, dv, ADM, ADM_1
 COMPLEX(q) :: VDM(SIZE_REC)
 REAL(q), ALLOCATABLE :: VDM_WORK(:)

 ALLOCATE(VDM_WORK(SIZE_REAL))

 dv =  LATT_CUR_NLN%OMEGA / GRIDC_NLN%NPLWV

 VDM_WORK = 0._q

 ADM = 0._q

 DO I = 1, SIZE_REAL

    NDIM_E = EFIELD(I) * pkT
    ADM_1 = Nsol_k * kb_k * T_k *  &
          & ( ( NDIM_E * LANG_NLN(I) - LOG(SINH(NDIM_E)/NDIM_E) ) + &
          & (kb_k * T_k / Nsol_k /p0_k/p0_k) * (einfty_k - 1._q *(NDIM_E**2/2) ) )
    ADM = ADM + SHAPE_FUNC(I) * ADM_1 *dv
    VDM_WORK(I) = SHAPE_FUNC_DER(I) * ADM_1

 ENDDO

 CALL PACK_R2C(VDM_WORK, VDM, GRIDC_NLN)

 CALL FFT3D(VDM, GRIDC_NLN, -1)
 
 DEALLOCATE(VDM_WORK)

 RETURN

END SUBROUTINE VA_DMED



!set the the reciprocal lattice vectors
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
 
         ENDDO row
      ENDDO col

 RETURN

END SUBROUTINE SET_GVEC_NLN




!>========================== SUBROUTINE VCORRECTION_NLN_PCM ======================
!!
!!
!!
!!===============================================================================

SUBROUTINE VCORRECTION_NLN_PCM(TOTN, PS_TOTN, VCORR, VEL, ATOT )

 IMPLICIT NONE

 COMPLEX(q), INTENT(OUT) ::  VCORR(SIZE_REC)
 COMPLEX(q),INTENT(OUT) :: VEL(SIZE_REC)
 COMPLEX(q),INTENT(IN) :: TOTN(SIZE_REC), PS_TOTN(SIZE_REC) 
 COMPLEX(q), ALLOCATABLE :: PS_TOTN_TMP(:)

 INTEGER :: I,J, NCGSOL
 REAL(q)  :: ATOT

 ALLOCATE(PS_TOTN_TMP(SIZE_REC))

 VEL = 0._q
 VCORR = 0._q
 ATOT = 0._q

!set R_PS_TOTN

 PS_TOTN_TMP = PS_TOTN / LATT_CUR_NLN%OMEGA

 CALL FFT3D(PS_TOTN_TMP, GRIDC_NLN, 1)

 CALL UNPACK_C2R(PS_TOTN_TMP, 1.0_q, PS_TOTN_TMP, 0.0_q, R_PS_TOTN, GRIDC_NLN)
!end of setting R_PS_TOTN

 CALL MINIMIZE_NLN(TOTN, NCGSOL) 

 CALL VACORR_ELECTROSTATIC(TOTN, VEL, VCORR, ATOT)

 CALL VACORR_NON_ELECTROSTATIC(VCORR, ATOT)

 PK_COUNTER = PK_COUNTER + 1

 DEALLOCATE(PS_TOTN_TMP)

 RETURN

END SUBROUTINE VCORRECTION_NLN_PCM



END MODULE POT_NLN_K
