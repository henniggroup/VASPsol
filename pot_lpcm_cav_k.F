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

MODULE POT_LPCM_K

 USE prec
 USE constant 
 USE base
 USE mpimy 
 USE poscar, ONLY : type_info 
 USE lattice, ONLY : latt, DIRKAR 
 USE mgrid, ONLY : grid_3d
 USE pseudo, ONLY : potcar

 PUBLIC :: Vcorrection_lpcm, MY_D_PROD
 PRIVATE :: CREATECAVITY, MINIMIZE, Leps, EPS_POT
 PRIVATE :: SHAPE_GRADN, NABLAN, LAPLN, SimplePoisson
 PRIVATE :: PACK_R2C, UNPACK_C2R

 INTEGER, PRIVATE, SAVE                 :: PK_counter = 1
 REAL(q), PRIVATE, SAVE                 :: CHGMIN = 1E-10_q
 COMPLEX(q), PRIVATE, ALLOCATABLE, SAVE :: Sol_phi(:)
 REAL(q), ALLOCATABLE, PRIVATE, SAVE    :: sigma_rc_k(:)
 REAL(q), PRIVATE, SAVE                 :: sigma_nc_k,sigma_k,nc_k,eb_k,CORE_C, tau,EDIFFSOL
 TYPE (in_struct), PRIVATE, SAVE        :: SP_IO

 CONTAINS

!>==========================SUBROUTINE Vcorrection_lpcm======================
!!
!!Solves the poisson's eqn, returns Vcorr and Ecorr
!!
!!======================================================================

SUBROUTINE Vcorrection_lpcm(IO,SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,TOTN,PS_TOTN,Vcorr,Vel,Atot,sigma_nc_k1,sigma_k1,nc_k1,eb_k1,CORE_C1,sigma_rc_k1,tau1,ediffs)

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

 sigma_nc_k = sigma_nc_k1 
 sigma_k =  sigma_k1
 nc_k =  nc_k1
 eb_k =  eb_k1
 CORE_C =  CORE_C1
 sigma_rc_k =  sigma_rc_k1
 tau = tau1
 EDIFFSOL = ediffs

 Ecorr1 = 0.0_q
 Ecorr2 = 0.0_q
 phi_tilda = 0._q
 phi = 0._q
 Vel = 0._q
 Vcorr = 0._q
 Vcav = 0._q
 SP_IO = IO

!-------------------------------------------------------------
!create the positon-dependent dielectric constant in real space
!ps_totn is unscaled and in real space after this call
!returns Vcav(cavitaion pot in rec space) and Acav
!-------------------------------------------------------------
 CALL CREATECAVITY(SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,PS_TOTN,SP_d_eps,Vcav,Acav)

!-------------------------------------------------------------
!Solves the generalized poisson eqn using CG, returns phi_tilda scaled by edeps/cell vol
!-------------------------------------------------------------
 CALL MINIMIZE(SP_GRIDC,SP_LATT_CUR,SP_P,SP_T_INFO,SP_d_eps,TOTN,phi_tilda,NCGSOL)

!-------------------------------------------------------------
!compute the interface contribution to the local pot, Vel
!ps_totn is unscaled and in real space
!-------------------------------------------------------------
 CALL UNPACK_C2R(PS_TOTN,1.0_q,PS_TOTN,0.0_q,RWORK,SP_GRIDC)

 CALL EPS_POT(SP_GRIDC,SP_LATT_CUR,SP_T_INFO,SP_P,RWORK,phi_tilda,SP_d_eps,Vel)

!-------------------------------------------------------------
! multiply by 0.5*e/permittivity of free space
!-------------------------------------------------------------
 Vel = 0.5 * Vel / EDEPS

!-------------------------------------------------------------
!compute phi, pot due to (n+N) with epsilon=1
!returns properly scaled phi
!-------------------------------------------------------------
 CALL SimplePoisson(SP_GRIDC,SP_LATT_CUR,TOTN,phi,1)

!-----------------------------------------------------------------------------
!compute V_correction,Vcorr
!the hartree potential difference between the vacuum case and the solvent case
!-----------------------------------------------------------------------------
 Vcorr = phi_tilda - phi

!add the cavitation correction
! Vcorr = Vcorr + Vcav
!both Vel and Vcav are the derivatives of the free enrgy wrt the n_cavity=n+pseudo core
!the reason for seperating Vcorr(coulomb) and Vel is becasue for the force corrections due Vcorr we need to multiply by the derivative of total charge density with respect to the positions wheras for Vel we need to multiply by the derivativer of the n_cavity wrt the positions
 Vel = Vel + Vcav

!-------------------------------------------------------------
!compute the components of E_correction, Ecorr
!all energies are computed with chargedensities scaled by the lattice volme inorder to be consistent with vasp
!-------------------------------------------------------------
!0.5*phi_tilda*(n+N)
!-------------------------------------------------------------
! Ecorr1 = 0.5*DOT_PRODUCT(TOTN,phi_tilda)
 CALL MY_D_PROD(Ecorr1, TOTN, phi_tilda, SP_GRIDC)
 Ecorr1 = Ecorr1 * 0.5_q
 CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr1,1))

!-------------------------------------------------------------
!0.5*phi*(n+N)
!-------------------------------------------------------------
!  Ecorr2 = 0.5*DOT_PRODUCT(TOTN,phi)
 CALL MY_D_PROD(Ecorr2, TOTN, phi, SP_GRIDC)
 Ecorr2 = Ecorr2 * 0.5_q
 CALLMPI( M_sum_d(SP_GRIDC%COMM,Ecorr2,1))

!electrostatic contribution
 Ael = Ecorr1 - Ecorr2

!add the cavitaion energy
 Atot = Ael + Acav
 
 333 FORMAT('SOL:',' ',I3,'  ',E12.5,'   ',E12.5,'   ',E12.5,'   ',I3)

 IF (SP_IO%IU0>=0) THEN

   WRITE(17,333)  PK_COUNTER,Ael,Acav,Atot,NCGSOL
   WRITE(SP_IO%IU0,333)  PK_COUNTER,Ael,Acav,Atot,NCGSOL

 ENDIF

 PK_COUNTER = PK_COUNTER + 1
      
 DEALLOCATE(SP_d_eps,phi_tilda,phi,Vcav,RWORK)

 RETURN

END SUBROUTINE Vcorrection_lpcm

!>=======================SUBROUTINE CREATECAVITY======================
!!
!!initializes epsilon in the grid
!!calculates the cavitation energy and cavitation potential
!!Acav and Vcav
!!
!!======================================================================

SUBROUTINE CREATECAVITY(CE_GRIDC,CE_LATT_CUR,CE_P,CE_T_INFO,CE_PS_TOTN,CE_d_eps,Vcav,Acav)  

 IMPLICIT NONE

 TYPE(type_info),INTENT(IN) :: CE_T_INFO
 TYPE(potcar),INTENT(IN)    :: CE_P (CE_T_INFO%NTYP)
 TYPE(grid_3d),INTENT(IN)   :: CE_GRIDC
 TYPE(latt),INTENT(IN)      :: CE_LATT_CUR

 COMPLEX(q)             :: CE_PS_TOTN(CE_GRIDC%MPLWV), Vcav(CE_GRIDC%MPLWV)
 REAL(q),INTENT(OUT)    :: CE_d_eps(DIMREAL(CE_GRIDC%MPLWV))

 REAL(q),ALLOCATABLE        :: CE_S(:), VWORK(:) ! shape func & Vcav
 REAL(q),ALLOCATABLE        :: CE_gn(:,:), CE_gn2_tmp(:) !\nabla n and |\naba n|^2
 COMPLEX(q),ALLOCATABLE     :: CE_gn2(:) ! |\naba n|^2
 REAL(q),ALLOCATABLE        :: CE_ggn(:,:) !\nabla |\nabla n|
 REAL(q),ALLOCATABLE        :: CE_ggn2(:) !| \nabla |\nabla n| |^2
 REAL(q),ALLOCATABLE        :: CE_divn(:) ! laplacian of n
 REAL(q),ALLOCATABLE        :: CE_term1(:) ! grad shape func wrt n
! REAL(q), ALLOCATABLE       :: RWORK(:)
 RGRID, ALLOCATABLE       :: RWORK(:)

 REAL(q), EXTERNAL      :: ERRFC !compimentary error function implemented in errf.F in vasp.5.lib
 INTEGER(q)             :: I, J, NPC, NPR, RLNP, RCNP
 REAL(q)                :: Acav,CE_qs,dv,tmp,tmp1,tmp2

 NPC =  CE_GRIDC%MPLWV
 RLNP = CE_GRIDC%RL%NP
 RCNP = CE_GRIDC%RC%NP
 NPR =  DIMREAL(CE_GRIDC%MPLWV) !NPR = NPC, in complex mode

 IF (SP_IO%IU0>=0) print *, 'NPLWV, MPLWV, RL%NP, RC%NP', CE_GRIDC%NPLWV, CE_GRIDC%MPLWV, CE_GRIDC%RL%NP, CE_GRIDC%RC%NP

 ALLOCATE(CE_gn2(NPC)) !complex
 ALLOCATE(CE_S(RLNP))
 ALLOCATE(CE_gn2_tmp(NPR), CE_gn(3, RLNP) )
 ALLOCATE(CE_ggn2(NPR), CE_ggn(3, RLNP) )
 ALLOCATE(CE_divn(NPR), CE_term1(RLNP))
 ALLOCATE(RWORK(NPR), VWORK(NPR))

!volume element
 dv = CE_LATT_CUR%OMEGA / CE_GRIDC%NPLWV

!remove the scaling in the reciprocal charge density
 CE_PS_TOTN = CE_PS_TOTN / CE_LATT_CUR%OMEGA

!charge density in rec space, returns \nabla n and |\nabla n|^2 in real space
 CALL NABLAN(CE_PS_TOTN, CE_GRIDC, CE_LATT_CUR, CE_gn2_tmp, CE_gn)

!charge density in rec space, returns laplacian of n, Lap(n)
!Lap(n) = (\nabla \dot \nabla n) in real space
 CALL LAPLN(CE_PS_TOTN, CE_GRIDC, CE_LATT_CUR, CE_divn)

!cavitaion pot in real space
 Vcav = 0._q
 VWORK = 0._q

!add -Lap(n)/|\nabla n| to Vcav and copy \sqrt(|\nabla n|^2) to CE_gn2_tmp
 DO I=1,DIMREAL(CE_GRIDC%MPLWV)
   tmp = SQRT(MAX(REAL(CE_gn2_tmp(I),KIND=q),CHGMIN))
   VWORK(I) = VWORK(I) - (CE_divn(I))/tmp
   CE_gn2_tmp(I) = tmp
 ENDDO

!CE_PS_TOTN to the real space
 CALL FFT3D(CE_PS_TOTN,CE_GRIDC,1)

!unpack complex array to a real array
!in complex mode, just copies CE_PS_TOTN to RWORK
 CALL UNPACK_C2R(CE_PS_TOTN, 1.0_q, CE_PS_TOTN, 0.0_q, RWORK, CE_GRIDC)

!charge density must be in real space(unscaled), 
!returns derivative of shape function wrt the charge density in the real space,
! CE_term1 = gamma*A / n, where
! gamma * A = exp( -( log(n/n_c) )^2 / (2 sigma^2) ) / ( sigma * sqrt(2*pi) )
 CALL SHAPE_GRADN(CE_GRIDC, CE_T_INFO, CE_P, CE_LATT_CUR, RWORK, CE_term1)

!-------------------------------------------------------------
! epsilon  = 1.0 + (eb_k - 1) * shape function
!-------------------------------------------------------------
 DO I=1, CE_GRIDC%RL%NP
   CE_S(I) = ERRFC( log(MAX(REAL(RWORK(I),KIND=q),CHGMIN)/nc_k) / SQRT(2._q)/sigma_k ) / 2._q
   CE_d_eps(I)= 1._q + (eb_k - 1._q) * CE_S(I)
 END DO

!quantum surface area
 CE_qs = 0._q

!quantum surface area, integral of (gamma*A / n) * |\nabla n| 
!= CE_term1 * \sqrt(CE_gn2_tmp)
 DO I=1, CE_GRIDC%RL%NP
   CE_qs = CE_qs + (CE_term1(I)) * (CE_gn2_tmp(I))
 ENDDO

 CE_qs = CE_qs * dv

 CALLMPI( M_sum_d(CE_GRIDC%COMM,CE_qs,1))

!cavitation energy
 Acav = tau * CE_qs

#ifdef debugsol
 IF (SP_IO%IU0>=0) THEN
   print *,'tau ',tau
   print *,'surface area ',CE_qs
   print *,'Acav ',Acav
 ENDIF
#endif

!1/ |nabla n|
 DO I=1, DIMREAL(CE_GRIDC%MPLWV)
   CE_gn2_tmp(I) = 1._q / ( MAX(CE_gn2_tmp(I), CHGMIN) ) 
 ENDDO

!-------------------------------------------------------------
!Compute cavitation potential, Vcav
!-------------------------------------------------------------
!packs the real array into a complex one
!in complex mode, same as  CE_gn2 = CMPLX(CE_gn2_tmp)
 CALL PACK_R2C(CE_gn2_tmp,CE_gn2,CE_GRIDC)

!1/ |nabla n| to reciprocal space to take the gradient
 CALL FFT3D(CE_gn2,CE_GRIDC,-1)
 CALL SETUNB(CE_gn2,CE_GRIDC)
 CE_gn2 = CE_gn2/CE_GRIDC%NPLWV

!\nabla(1 / |\nabla n|), CE_ggn in real space
 CALL NABLAN(CE_gn2, CE_GRIDC, CE_LATT_CUR, CE_ggn2, CE_ggn)

!add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
!and multiply by CE_term1 = gamma*A/n in real space
 DO I=1, CE_GRIDC%RL%NP
  tmp = ( CE_gn(1,I)*CE_ggn(1,I) + CE_gn(2,I)*CE_ggn(2,I) + CE_gn(3,I)*CE_ggn(3,I) ) * (CE_term1(I))
  VWORK(I) = VWORK(I) - tmp
 ENDDO

 CALL PACK_R2C(VWORK, Vcav, CE_GRIDC)

!Vcav to reciprocal space
 CALL FFT3D(Vcav,CE_GRIDC,-1)
 CALL SETUNB(Vcav,CE_GRIDC)
 Vcav = Vcav/CE_GRIDC%NPLWV

!multiply by the surface tension
 Vcav = Vcav * tau

 DEALLOCATE(CE_S, RWORK, VWORK)
 DEALLOCATE(CE_gn2, CE_gn2_tmp)
 DEALLOCATE(CE_gn, CE_ggn2)
 DEALLOCATE(CE_ggn, CE_divn, CE_term1)

 RETURN

END SUBROUTINE CREATECAVITY

!>=====================SUBROUTINE MINIMIZE==============================
!!
!!routine for solving poissons eqn using conjugate gradient or 
!!steepest descent
!!
!!======================================================================

SUBROUTINE MINIMIZE(MN_GRIDC,MN_LATT_CUR,MN_P,MN_T_INFO,MN_d_eps,MN_TOTN,&
                    &MN_phi,NCGSOL)

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
 COMPLEX(q),ALLOCATABLE :: MN_GSQU(:)

 COMPLEX(q),ALLOCATABLE :: gradphi_x_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_y_in(:)
 COMPLEX(q),ALLOCATABLE :: gradphi_z_in(:)

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
 ALLOCATE(gradphi_z_in(MN_GRIDC%MPLWV))


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
   GX=(MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(1,1)+MN_GRIDC%LPCTY(N2)*MN_LATT&
       &_CUR%B(1,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(1,3))
   GY=(MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(2,1)+MN_GRIDC%LPCTY(N2)*MN_LATT&
       &_CUR%B(2,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(2,3))
   GZ=(MN_GRIDC%LPCTX(N1)*MN_LATT_CUR%B(3,1)+MN_GRIDC%LPCTY(N2)*MN_LATT&
       &_CUR%B(3,2)+MN_GRIDC%LPCTZ(N3)*MN_LATT_CUR%B(3,3))
!set G=0 component 0
   IF ((MN_GRIDC%LPCTX(N1)==0).AND.(MN_GRIDC%LPCTY(N2)==0).AND.&
      &(MN_GRIDC%LPCTZ(N3)==0)) THEN
    MN_GSQU(I)=(0.0_q,0.0_q)
   ELSE
! 1/G^2
    MN_GSQU(I)=1/(GX*GX+GY*GY+GZ*GZ)/(TPI**2)
   ENDIF
 
 ENDDO row
ENDDO col

!*********************************************************
! start of Pre-conditioned Conjugate Gradient
!**********************************************************

!-------------------------------------------------------------
!intial guess for phi
!-------------------------------------------------------------
 IF (PK_counter ==1) THEN
   DO I=1,MN_GRIDC%RC%NP
     MN_phi(I) = REAL(MN_GSQU(I))*(MN_TOTN(I))
   ENDDO
 ELSE
   MN_phi = Sol_phi !phi from previous vasp iteration, unscaled
 ENDIF

!-------------------------------------------------------------
!Leps returns Laplasian Lp=div ( epsilon * grad ) phi in k space
!epsilon is in real space
!Lp = - A*phi
!-------------------------------------------------------------
 CALL Leps(MN_GRIDC,MN_LATT_CUR,MN_T_INFO,MN_phi,MN_d_eps,gradphi_x_in,gradphi_y_in,gradphi_z_in,MN_Lp)

!-------------------------------------------------------------
! the residue
!r=A*phi+CHTOT+N
!-------------------------------------------------------------
 DO I=1,MN_GRIDC%RC%NP
   MN_RESID(I)=MN_Lp(I) + MN_TOTN(I)
 ENDDO

!-------------------------------------------------------------
!preconditioning of the residue, MN_Z = invL(r)
!z=r/G^2
!-------------------------------------------------------------
 MN_Z=0._q

 DO I=1,MN_GRIDC%RC%NP
   MN_Z(I)=REAL(MN_GSQU(I))*MN_RESID(I)
 ENDDO

!-------------------------------------------------------------
!without preconditioning
!-------------------------------------------------------------
!MN_Z=MN_RESID

! rinvLr = DOT_PRODUCT(MN_RESID,MN_Z)
 CALL MY_D_PROD(rinvLr, MN_RESID, MN_Z, MN_GRIDC)

 CALLMPI( M_sum_d(MN_GRIDC%COMM,rinvLr,1))

!-------------------------------------------------------------
!square of the residue
!r^2
!-------------------------------------------------------------
! r2 = DOT_PRODUCT(MN_RESID,MN_RESID)
 CALL MY_D_PROD(r2, MN_RESID, MN_RESID, MN_GRIDC)
 CALLMPI( M_sum_d(MN_GRIDC%COMM,r2,1))

 r20 = r2

 MN_d = MN_Z

!-------------------------------------------------------------
!CG loop
!-------------------------------------------------------------

 DO WHILE ( (count < 20000) .AND. (sqrt(r2) > EDIFFSOL).AND. (sqrt(rinvLr) > 1e-10) )

   CALL Leps(MN_GRIDC,MN_LATT_CUR,MN_T_INFO,MN_d,MN_d_eps,gradphi_x_in,gradphi_y_in,gradphi_z_in,MN_Lp)

!   MN_alpha = DOT_PRODUCT(MN_d,MN_Lp)
   CALL MY_D_PROD(MN_alpha, MN_d, MN_Lp, MN_GRIDC)

   CALLMPI( M_sum_d(MN_GRIDC%COMM,MN_alpha,1))

   MN_alpha = 1._q/MN_alpha

   MN_alpha = -rinvLr*MN_alpha

   MN_phi = MN_phi + MN_alpha * MN_d

   MN_RESID = MN_RESID +  (MN_alpha*MN_Lp)

   DO I=1,MN_GRIDC%RC%NP
     MN_Z(I) = (REAL(MN_GSQU(I)))*MN_RESID(I)
   ENDDO

   MN_beta = 1/rinvLr

   CALL MY_D_PROD(rinvLr, MN_RESID, MN_Z, MN_GRIDC)

   CALLMPI( M_sum_d(MN_GRIDC%COMM,rinvLr,1))

   MN_beta = MN_beta * rinvLr

   MN_d = MN_Z + MN_beta*MN_d

   r2 = 0._q

   CALL MY_D_PROD(r2, MN_RESID, MN_RESID, MN_GRIDC)
!   DO I=1,MN_GRIDC%RC%NP
!     r2 = r2 + MN_RESID(I)*CONJG(MN_RESID(I))
!   ENDDO

   CALLMPI( M_sum_d(MN_GRIDC%COMM,r2,1))

#ifdef debugsol
   IF (SP_IO%IU0>=0)   WRITE(SP_IO%IU0,*) SQRT(rinvLr), SQRT(r2)
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
!multiply by 1/epsilon_0and divide by volume for CHTOT
!-------------------------------------------------------------
 MN_phi = MN_phi * EDEPS / MN_LATT_CUR%OMEGA

!-------------------------------------------------------------
!dealloacte all the work arrays
!-------------------------------------------------------------
 DEALLOCATE(MN_RESID,MN_Z,MN_Lp,MN_GSQU)
 DEALLOCATE(MN_d)
 DEALLOCATE(gradphi_x_in,gradphi_y_in,gradphi_z_in)

 RETURN
END SUBROUTINE MINIMIZE

!>======================SUBROUTINE Leps==================================
!!
!!subroutine that evaluates the fourier transform of div(epsilon*grad phi)
!!Ax in Ax=b
!!
!!=======================================================================

SUBROUTINE Leps(LPS_GRIDC,LPS_LATT_CUR,LPS_T_INFO,LPS_phi,LPS_d_eps,&
                &gradphi_x_in,gradphi_y_in,gradphi_z_in,LPS_Lp)

 IMPLICIT NONE

 TYPE (grid_3d) LPS_GRIDC
 TYPE (type_info) LPS_T_INFO
 TYPE (latt) LPS_LATT_CUR

 REAL(q),INTENT(IN) ::  LPS_d_eps(DIMREAL(LPS_GRIDC%MPLWV))
 COMPLEX(q),INTENT(IN) ::  LPS_phi(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(OUT) :: LPS_Lp(LPS_GRIDC%MPLWV)
 REAL(q) RINPL,GX,GY,GZ

 COMPLEX(q),INTENT(INOUT) :: gradphi_x_in(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(INOUT) :: gradphi_y_in(LPS_GRIDC%MPLWV)
 COMPLEX(q),INTENT(INOUT) :: gradphi_z_in(LPS_GRIDC%MPLWV)

 INTEGER(q) N
 INTEGER I,J,k,N1,N2,N3,NC
      
 N=LPS_GRIDC%NPLWV


!-------------------------------------------------------------------------
!Fourier transform of grad phi
!-------------------------------------------------------------------------
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

!-------------------------------------------------------------------------
!gradphi to real space
!-------------------------------------------------------------------------
 CALL FFT3D(gradphi_x_in,LPS_GRIDC,1) 
 CALL FFT3D(gradphi_y_in,LPS_GRIDC,1) 
 CALL FFT3D(gradphi_z_in,LPS_GRIDC,1) 

!-------------------------------------------------------------------------
!multiply grad phi with epsilon in real space
!-------------------------------------------------------------------------
#ifdef realmode
k =1
 DO j=1,DIMREAL(LPS_GRIDC%MPLWV),2
   gradphi_x_in(k) = CMPLX(REAL(gradphi_x_in(k),KIND=q)*LPS_d_eps(j),AIMAG(gradphi_x_in(k))*LPS_d_eps(j+1),KIND=q)
   gradphi_y_in(k) = CMPLX(REAL(gradphi_y_in(k),KIND=q)*LPS_d_eps(j),AIMAG(gradphi_y_in(k))*LPS_d_eps(j+1),KIND=q)
   gradphi_z_in(k) = CMPLX(REAL(gradphi_z_in(k),KIND=q)*LPS_d_eps(j),AIMAG(gradphi_z_in(k))*LPS_d_eps(j+1),KIND=q)
k = k+1
 ENDDO
#else
 DO j=1,LPS_GRIDC%RL%NP
   gradphi_x_in(j) = gradphi_x_in(j)*LPS_d_eps(j)
   gradphi_y_in(j) = gradphi_y_in(j)*LPS_d_eps(j)
   gradphi_z_in(j) = gradphi_z_in(j)*LPS_d_eps(j)
 ENDDO
#endif
!-------------------------------------------------------------------------
!epsilon*grad phi in kspace
!-------------------------------------------------------------------------
 CALL FFT3D(gradphi_x_in,LPS_GRIDC,-1)   
 CALL FFT3D(gradphi_y_in,LPS_GRIDC,-1)   
 CALL FFT3D(gradphi_z_in,LPS_GRIDC,-1)   

 CALL SETUNB(gradphi_x_in,LPS_GRIDC)   
 CALL SETUNB(gradphi_y_in,LPS_GRIDC)   
 CALL SETUNB(gradphi_z_in,LPS_GRIDC)   

!-------------------------------------------------------------------------
!scaling
!-------------------------------------------------------------------------
 DO I=1,LPS_GRIDC%RC%NP
   gradphi_x_in(I) = gradphi_x_in(I)/N
   gradphi_y_in(I) = gradphi_y_in(I)/N
   gradphi_z_in(I) = gradphi_z_in(I)/N
 ENDDO

!-------------------------------------------------------------------------
!div(epsilon*grad phi) in kspace
!-------------------------------------------------------------------------
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
              &GZ*gradphi_z_in(I))

 ENDDO
 
 END SUBROUTINE Leps

!>===========================SUBROUTINE  EPS_POT=======================
!!
!!subroutine that evaluates the V_epsilon
!!returns -(depsilon(n)/dn)*|grad phi|^2 in fourier space
!!though the multiplication is done in real space
!!
!!====================================================================

SUBROUTINE  EPS_POT(EPP_GRIDC, EPP_LATT_CUR,EPP_T_INFO,EPP_P,EPP_PS_TOTN,EPP_phi,EPP_deps,EPP_CWORK)

 IMPLICIT NONE

 TYPE (grid_3d),INTENT(IN)   :: EPP_GRIDC
 TYPE (latt),INTENT(IN)      :: EPP_LATT_CUR
 TYPE (type_info),INTENT(IN) :: EPP_T_INFO
 TYPE(potcar),INTENT(IN)     :: EPP_P (EPP_T_INFO%NTYP)

 COMPLEX(q),INTENT(OUT) :: EPP_CWORK(EPP_GRIDC%MPLWV)
 RGRID,INTENT(IN)       :: EPP_PS_TOTN(DIMREAL(EPP_GRIDC%MPLWV))
 COMPLEX(q),INTENT(IN)  :: EPP_phi(EPP_GRIDC%MPLWV)
 REAL(q),INTENT(IN)     :: EPP_deps(EPP_GRIDC%MPLWV)

 REAL(q), ALLOCATABLE   :: eprime(:)
 REAL(q),ALLOCATABLE    :: phisq(:), RWORK(:)
 REAL(q)                :: EPP_CWORK_temp
 INTEGER(q)             :: I,J  

 ALLOCATE(eprime(EPP_GRIDC%RL%NP), phisq(DIMREAL(EPP_GRIDC%MPLWV)), RWORK(DIMREAL(EPP_GRIDC%MPLWV)))

!-------------------------------------------------------------------------
!SHAPE_GRADN returns eprime, the derivative of shape function wrt the charge density in real space
!-------------------------------------------------------------------------
   CALL SHAPE_GRADN(EPP_GRIDC, EPP_T_INFO, EPP_P, EPP_LATT_CUR, EPP_PS_TOTN, eprime)
   eprime = eprime * (eb_k - 1._q)

!-------------------------------------------------------------------------
!EPP_phi,the solution to the generalized poisson equation,
! is in reciprocal space, 
!NABLAN returns |\nabla phi| ^2 in real space
!-------------------------------------------------------------------------
   CALL NABLAN(EPP_phi,EPP_GRIDC,EPP_LATT_CUR,phisq)

!-------------------------------------------------------------------------
! eps potential in real space
!-------------------------------------------------------------------------
   EPP_CWORK_temp = 0._q

!-------------------------------------------------------------------------
!multiply eprime and phisq in real space
!-------------------------------------------------------------------------
   RWORK = 0._q

   DO I=1, EPP_GRIDC%RL%NP
      EPP_CWORK_temp = (eprime(I)) * (phisq(I))
      RWORK(I) = EPP_CWORK_temp
   ENDDO

 CALL PACK_R2C(RWORK, EPP_CWORK, EPP_GRIDC)
   
!-------------------------------------------------------------------------
!take the eps potential to reciprocal space
!-------------------------------------------------------------------------
   CALL FFT3D(EPP_CWORK,EPP_GRIDC,-1)
   CALL SETUNB(EPP_CWORK,EPP_GRIDC)
   EPP_CWORK = EPP_CWORK/EPP_GRIDC%NPLWV

 DEALLOCATE(eprime,phisq,RWORK)
 
 RETURN

END SUBROUTINE EPS_POT

!>======================SUBROUTINE SHAPE_GRADN============================
!!
!!calculates first derivative of the shape function  wrt CHTOT in realspace
!!
!!=====================================================================

SUBROUTINE SHAPE_GRADN(EPR_GRIDC, EPR_T_INFO, EPR_P, EPR_LATT_CUR, EPR_PS_TOTN, eprime)

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

!>==========================SUBROUTINE NABLAN===============================
!!
!! input : D_PHI , in reciprocal space, unscaled
!! calculate |\nabla D_PHI |^2 = D_PHI2 in real space
!! optional arg, D_grad = \nabla D_PHI, real space vector
!!
!!======================================================================

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

!-------------------------------------------------------------------------
! x-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GX = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(1,1)+D_GRIDC%LPCTY(N2)*D_LATT_CUR%B(1,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(1,3))

   CR(I) = D_PHI(I) * GX * CITPI
 ENDDO
!-------------------------------------------------------------------------
! gradphi_x to real space:
!-------------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
   tmp = REAL(RWORK(I),KIND=q)
   D_PHI2(I) = D_PHI2(I) + tmp * tmp
   IF (PRESENT(D_grad)) D_grad(1,I) = tmp
  ENDDO
!-------------------------------------------------------------------------
! y-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GY = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(2,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(2,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(2,3))

   CR(I) = D_PHI(I) * GY * CITPI
 ENDDO

!-------------------------------------------------------------------------
! grad_y to real space:
!-------------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
   tmp = REAL(RWORK(I),KIND=q)
   D_PHI2(I) = D_PHI2(I) + tmp * tmp
   IF (PRESENT(D_grad)) D_grad(2,I) = tmp
  ENDDO
!-------------------------------------------------------------------------
! z-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GZ = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(3,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(3,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(3,3))

   CR(I)=D_PHI(I) * GZ * CITPI
 ENDDO

!-------------------------------------------------------------------------
! grad_z to real space:
!-------------------------------------------------------------------------
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

!>==========================SUBROUTINE LAPLN===============================
!!
!! calculate the Laplacian of D_PHI in real space
!! D_PHI, the input in reciprocal space
!! D_PHI2, the output in real space
!!======================================================================

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

!-------------------------------------------------------------------------
! x-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GX = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(1,1)+D_GRIDC%LPCTY(N2)*&
       &D_LATT_CUR%B(1,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(1,3))

   CR(I) = D_PHI(I) * (GX*CITPI)**2.0
 ENDDO

!-------------------------------------------------------------------------
! gradphi to real space:
!-------------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
    D_PHI2(I) = D_PHI2(I) + REAL(RWORK(I),KIND=q)
 ENDDO

!-------------------------------------------------------------------------
! y-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GY = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(2,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(2,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(2,3))

   CR(I) = D_PHI(I) * (GY*CITPI)**2
 ENDDO

!-------------------------------------------------------------------------
! grad_y to real space:
!-------------------------------------------------------------------------
 CALL SETUNB(CR,D_GRIDC)
 CALL FFT3D(CR,D_GRIDC,1)

 CALL UNPACK_C2R(CR, 1.0_q, CR, 0.0_q, RWORK, D_GRIDC)

 DO I=1, D_GRIDC%RL%NP
    D_PHI2(I) = D_PHI2(I) + REAL(RWORK(I),KIND=q)
 ENDDO
!-------------------------------------------------------------------------
! z-component:
!-------------------------------------------------------------------------
 DO I=1,D_GRIDC%RC%NP
   N1= MOD((I-1),D_GRIDC%RC%NROW) +1
   NC= (I-1)/D_GRIDC%RC%NROW+1
   N2= D_GRIDC%RC%I2(NC)
   N3= D_GRIDC%RC%I3(NC)

   GZ = (D_GRIDC%LPCTX(N1)*D_LATT_CUR%B(3,1)+D_GRIDC%LPCTY(N2)*&
         &D_LATT_CUR%B(3,2)+D_GRIDC%LPCTZ(N3)*D_LATT_CUR%B(3,3))

   CR(I)=D_PHI(I) * (GZ*CITPI)**2.0
 ENDDO

!-------------------------------------------------------------------------
! grad_z to real space:
!-------------------------------------------------------------------------
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

!>==========================SUBROUTINE MY_D_PROD=================
!!
!! helper routine for dot product
!! also used in pot_k.F
!!
!!==============================================================

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

!>==========================SUBROUTINE PACK_R2C=================
!!
!!  helper routine to pack real array(RA) to a complex array(CB)
!!
!!==============================================================

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

!>==========================SUBROUTINE UNPACK_C2R=================
!!
!!  helper routine to unpack and maybe scale
!!  complex array(A & B) to real array(C)
!!
!!==============================================================

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

END MODULE POT_LPCM_K
!*******************************************************************************************
!*****************************************END OF MODULE POT_LPCM_K ******************************
!*******************************************************************************************
