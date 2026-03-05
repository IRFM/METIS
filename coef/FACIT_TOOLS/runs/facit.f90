! Software name : FACIT
! Authors : P. Maget and D. Fajardo, C. Angioni, P. Manas
! Copyright holders : Commissariat à l’Energie Atomique et aux Energies Alternatives (CEA), France, Max-Planck Institut für Plasmaphysik, Germany
! CEA and IPP authorize the use of the FACIT software under the CeCILL-C open source license https://cecill.info/licences/Licence_CeCILL-C_V1-en.html  
! The terms and conditions of the CeCILL-C license are deemed to be accepted upon downloading the software and/or exercising any of the rights granted under the CeCILL-C license.
! 
subroutine FACIT(nx, nth, xn, theta, &                                ! grid parameters
                 Za, Aa, Zi, Ai, &                                    ! impurity and main ion charge and mass
                 Te, Ti, Ne, Ni, Na, Machi, &                         ! plasma profiles
                 gradTi, gradNi, gradNa, &                            ! gradients
                 invaspct, B0, R0, qmag, dpsidx, FV, BV, RV, jacob, & ! equilibrium
                 AsymPhi, AsymN, &                                    ! input asymmetries
                 pol_asym, full_geom, rotation, regulopt, &           ! different options in the model
                 Flux_imp, Da, Vconv, nn, dmin, dmaj, &               ! main outputs: flux, transport coefficients and poloidal asymmetry
                 Da_PS, Da_BP, Da_CL, Ka_PS, Ka_BP, Ka_CL, &          ! supplementary outputs: PS, BP, CL & centrifugal components
                 Ha_PS, Ha_BP, Ha_CL,Va_PS,Va_BP,Va_CL)               ! supplementary outputs: PS, BP, CL

!*******************************************************************************
! Self-consistent calculation of collisional impurity transport and poloidal
! distribution of the impurity density
!------------
! References:
! Plasma Phys. Control. Fusion (2020) 62 105001 (doi: 10.1088/1361-6587/aba7f9)
! Modifications : 
! * 2022/01/10 : 1. correction to rotation term (P. Maget, corrigendum submitted to PPCF)
! *              2. last three outputs are Va_PS,Va_BP,Va_CL (one more than before)
! *              3. Use of a metric qmag_metrics consistent with full equilibrium when dpsidx available
! * 2022/01/21 : 1. correction of Cgeo_G for full_geom=.t. (line 288)
! *              2. correction of GGG for full_geom=.t. (line 778)
!* 2025/07/04 : 1. qmag input is the good quantity, compute dpsidx from it when not given in input
!*                          2. Correct mu_ie expression to be as in D. Fajardo's paper
!*******************************************************************************
! INPUTS:
! -----------> description [unit] {variable type, shape}
! ------------
! nx --------> size of radial arrays [-] {int}
! nth -------> size of poloidal arrays [-] {int}
! xn --------> radial coordinate [-] {arr, nx}
! theta -----> poloidal coordinate [-] {arr, nth}
! Za --------> impurity charge number [-] {arr, nx}
! Aa --------> impurity mass number [-] {float}
! Zi --------> main ion charge number [-] {float}
! Ai --------> main ion mass number [-] {float}
! Te --------> electron temperature [eV] {arr, nx}
! Ti --------> main ion temperature [eV] {arr, nx}
! Ne --------> electron density [1/m^3] {arr, nx}
! Ni --------> main ion denisty [1/m^3] {arr, nx}
! Na --------> impurity density [1/m^3] {arr, nx}
! Machi -----> Mach number of main ion [-] {arr, nx}
! gradTi ----> main ion temperature gradient [eV/-] {arr, nx}
! gradNi ----> main ion density gradient [1/m^3/-] {arr, nx}
! gradNa ----> impurity density gradient [1/m^3/-] {arr, nx}
! invaspct --> inverse aspect ratio [-] {float}
! B0 --------> magnetic field at magnetic axis [T] {float}
! R0 --------> major radius at magnetic axis [m] {float}
! qmag ------> safety factor [-] {arr, nx}
! dpsidx ----> radial derivative of poloidal flux [V*s/-] {arr, nx}
! FV --------> poloidal current flux function [T*m] {arr, nx}
! BV --------> magnetic field [T] {arr, (nx,nth)}
! RV --------> major radius [m] {arr, (nx,nth)}
! jacob -----> Jacobian of the coordinate system [m/T] {arr, (nx,nth)}
! AsymPhi ---> poloidal asymmetry of electrostatic potential [-] {arr, (nx,2)}
! AsymN -----> poloidal asymmetry of main ion density [-] {arr, (nx,2)}
! pol_asym --> poloidally symmetric (False) or asymmetric (True) system {boolean}
! full_geom -> analytical (False) or iterative (True) geometry {boolean}
! rotation --> consider toroidal rotation {boolean}
! regulopt --> options for iterative calculations [-] {arr, 4}
!*******************************************************************************
! OUTPUTS:
! --------
! Flux_imp --> surface-averaged impurity flux [1/(m^2 s)] {arr, nx}
! Da --------> diffusion coefficient [m^2/s] {arr, nx}
! Vconv -----> convective velocity [m/s] {arr, nx}
! nn --------> poloidal asymmetry of the impurity density Na/<Na> [-] {arr, (nx,nth)}
! dmin ------> horizontal asymmetry of the impurity density [-] {arr, nx}
! dmaj ------> vertical asymmetry of the impurity density [-] {arr, nx}
! Da_* ------> Pfirsch-Schlüter, Banana-Plateau and classical components of
!              the diffusion coefficient [m^2/s] {arr, nx}
! Ka_* ------> components of the coefficient of the main ion density gradient [m^2/s] {arr, nx}
! Ha_* ------> components of the coefficient of the main ion temperature gradient [m^2/s] {arr, nx}
! Va_* ------> Pfirsch-Schlüter, Banana-Plateau and classical components of
!              the pinch velocity [m/s] {arr, nx}!*******************************************************************************

  use constants, only: rkind, mp, me, q_e, eps_pi_fac, sqrt2
  implicit none

  !------------------------------ declarations --------------------------------

  ! INPUTS
  integer :: nx, nth
  real(rkind), dimension(nth) :: theta
  real(rkind), dimension(nx)  :: Te, Ti, Ne, Ni, Na, gradTi, gradNi, gradNa
  real(rkind), dimension(nx)  :: qmag, xn, dpsidx, FV, Za, Machi
  real(rkind) :: B0, R0, invaspct, Ai, Aa, Zi
  real(rkind), dimension(nx,nth) :: jacob, BV, RV, PhiV, NV
  real(rkind), dimension(nx,2) :: AsymPhi, AsymN
  real(rkind), dimension(4) :: regulopt
  logical :: pol_asym, full_geom, rotation

  ! OUTPUTS
  real(rkind), dimension(nx) :: Flux_imp, Da, Vconv, dmin, dmaj
  real(rkind), dimension(nx, nth) :: nn
  real(rkind), dimension(nx) :: Da_BP, Da_PS, Da_CL, Ka_BP, Ka_PS, Ka_CL, Ha_BP, Ha_PS, Ha_CL, Va_BP, Va_PS, Va_CL

  ! OTHER
  integer :: i, p, n, it, ix!, info, ierr, ierrmax
  real(rkind) :: amin, mi, ma, ftrap, ki_Redl, C2
  real(rkind), dimension(nx) :: grad_ln_ni, grad_ln_Ti, grad_ln_na, Ti15, Te15
  real(rkind), dimension(nx) :: epsk, eps15, eps2, ft, C0a, wca, deltaM, ki, g, dD2
  real(rkind), dimension(nx) :: f1, f2, f3, y1, y2, y3, y4, adps
  real(rkind), dimension(nx) :: LneeNRL, LneiNRL, LneimpNRL
  real(rkind), dimension(nx) :: LniiNRL, LniimpNRL, LnimpeNRL, LnimpimpNRL
  real(rkind), dimension(nx) :: Tauee, Tauei, Taueimp
  real(rkind), dimension(nx) :: Tauie, Tauii, Tauiimp
  real(rkind), dimension(nx) :: Tauimpe, Tauimpi, Tauimpimp
  real(rkind), dimension(nx) :: wee, wii, wimpimp, nuestar, nuistar, nuimpstar, rhoLimp2
  real(rkind), dimension(nx) :: L11impi, nuswca, mu_ie, alpha, Zeff
  real(rkind), dimension(nx) :: UU, GG, B2avg, dNH, dNV, dminphia, dmajphia
  real(rkind), dimension(nx) :: Cgeo_G, Cgeo_U, Cgeo_Gcl
  real(rkind), dimension(nx, nth) :: b2
  real(rkind), dimension(nx) :: b2navg, nb2avg, nNVavg, b2NVavg
  real(rkind), dimension(nx) :: K11a, K12a, K22a, K11i, K12i, K22i
  real(rkind), dimension(nx) :: Vra_BP, Vra_PS, Vra_CL
  real(rkind), dimension(nx) :: Ka, Ha

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI
  !(previous two lines here because "A specification statement cannot appear in the executable section")

  !---------------------------------------------------------------------------

  mi = Ai*mp ! main ion mass
  ma = Aa*mp ! impurity mass

  epsk  = xn*invaspct ! local inverse aspect ratio
  eps15 = (epsk + 1.0e-33)**1.5
  eps2  = (epsk + 1.0e-33)**2
  amin  = invaspct*R0 ! minor radius

  if (maxval(abs(dpsidx)).lt.1.0e-33_rkind) then
	dpsidx = amin**2*B0*xn/qmag
  endif

  Zeff = (Za**2*Na + Zi**2*Ni)/(Ne) ! effective charge

  do i = 1, nx
    ! logarithmic gradients
    grad_ln_ni(i) = gradNi(i)/(Ni(i) + 1.e-33)
    grad_ln_Ti(i) = gradTi(i)/(Ti(i) + 1.e-33)
    grad_ln_na(i) = gradNa(i)/(Na(i) + 1.e-33)
    ! trapped particle fraction
    ft(i) = ftrap(epsk(i))
  enddo


  if (.not.rotation) then
    Machi = 0.0_rkind
  endif

  deltaM = 2*(Aa/Ai)*Machi**2*epsk ! rotation strength parameter


  ! Coulomb Logarithms (from NRL formulary)
  LneeNRL = 23.5 - 0.5*log(Ne/1e6) + 1.25*log(Te) - sqrt(1e-5 + (1.0/16.0)*(log(Te) - 2)**2)

  do p = 1, nx

   if ((Ti(p)*me/(mi)<Te(p)).and.(Te(p)<10*Zi**2)) then
     LneiNRL(p) = 23 - 0.5*log(Zi**2*Ne(p)/1e6) + 1.5*log(Te(p))
   elseif ((Ti(p)*me/(mi)<10*Zi**2).and.(Te(p)>10*Zi**2)) then
     LneiNRL(p) = 24 - 0.5*log(Ne(p)/1e6) + log(Te(p))
   else
     LneiNRL(p) = 30 - log(Zi**2/Ai*(Ni(p)/1e6)**0.5) + 1.5*log(Ti(p))
   endif

   if ((Ti(p)*me/(Aa*mp)<Te(p)).and.(Te(p)<10*Za(p)**2)) then
     LneimpNRL(p) = 23. - log(Za(p)*(Ne(p)/1e6)**0.5) + 1.5*log(Te(p))
   elseif ((Ti(p)*me/(Aa*mp)<10.*Za(p)**2).and.(Te(p)>10.*Za(p)**2)) then
     LneimpNRL(p) = 24. - 0.5*log(Ne(p)/1e6) + log(Te(p))
   else
     LneimpNRL(p) = 30. - log(Za(p)**2/Aa*(Na(p)/1e6)**0.5) + 1.5*log(Ti(p))
   endif

  enddo

  LniiNRL     = 23. - log(Zi*Zi*sqrt(2*(Ni/1e6)*Zi**2)) + 1.5*log(Ti)
  LniimpNRL   = 23. - log(Zi*Za*sqrt((Ni/1e6)*Zi**2+(Na/1e6)*Za**2)) + 1.5*log(Ti)
  LnimpimpNRL = 23. - log(Za*Za*sqrt((Na/1e6)*Za**2+(Na/1e6)*Za**2)) + 1.5*log(Ti)


  ! Collision times (Braginskii)
  Ti15 = Ti**1.5
  Te15 = Te**1.5

  Tauee     = (eps_pi_fac*sqrt(me)*Te15)/(Ne*LneeNRL)
  Tauei     = (eps_pi_fac*sqrt(me)*Te15)/(Zi**2*Ni*LneiNRL)
  Taueimp   = (eps_pi_fac*sqrt(me)*Te15)/(Zi**2*Za**2*Na*LneimpNRL)

  Tauie     = (eps_pi_fac*sqrt(mi)*Ti15)/(Zi**2*Ne*LneiNRL)
  Tauii     = (eps_pi_fac*sqrt(mi)*Ti15)/(Zi**4*Ni*LniiNRL)
  Tauiimp   = (eps_pi_fac*sqrt(mi)*Ti15)/(Zi**2*Za**2*Na*LniimpNRL)

  Tauimpe   = (eps_pi_fac*sqrt(ma)*Ti15)/(Za**2*Ne*LneimpNRL)
  Tauimpi   = (eps_pi_fac*sqrt(ma)*Ti15)/(Zi**2*Za**2*Ni*LniimpNRL)
  Tauimpimp = (eps_pi_fac*sqrt(ma)*Ti15)/(Za**4*Na*LnimpimpNRL)


  ! impurity collision frequency
  L11impi = 1.0/(sqrt(1.0 + Aa/Ai)*Tauimpi)

  ! transit frequencies
  wee     = (2.0*q_e*Te/me)**0.5/(R0*qmag)
  wii     = (2.0*q_e*Ti/mi)**0.5/(R0*qmag)
  wimpimp = (2.0*q_e*Ti/ma)**0.5/(R0*qmag)

  ! collisionalities
  nuestar   = 1.0/((eps15 + 1.e-33)*wee*Tauee)
  nuistar   = 1.0/((eps15 + 1.e-33)*wii*Tauii)
  nuimpstar = 1.0/((eps15 + 1.e-33)*wimpimp*Tauimpimp)

  wca    = q_e*Za*B0/ma ! impurity cyclotron frequency
  nuswca = L11impi/wca ! ratio of impurity collision frequency to cyclotron frequency

  ! fitted factors
  call facs(nx, Za, ft, f1, f2, f3, y1, y2, y3, y4, adps)

  g = nuistar*eps15 ! collisionality parameter
  mu_ie = (96.0*sqrt2/125.0)*(1/Zi**2)*sqrt(me/mi)*(Ti15/Te15) ! ion-electron heat exchange term (Fülöp-Helander PoP '01)
  alpha = Na*Za**2/(Ni*Zi**2) ! impurity strength parameter

  do i = 1, nx
    C0a(i) = C2(alpha(i), g(i), f1(i), f2(i), Aa, Ai)/(1 + f3(i)*mu_ie(i)*g(i)**2) ! coefficient of ion heat flux in impurity-ion friction
    ki(i)  = ki_Redl(nuistar(i), ft(i), Zeff(i)) ! neoclassical ion flow coefficient
  enddo

  rhoLimp2 = (2.0*q_e*Ti/ma)/wca**2 ! Impurity Larmor radius (squared)

  ! thermodynamic gradients (for asymmetry calculations)
  UU  = -(Za/Zi)*(C0a + ki)*grad_ln_Ti
  GG  = grad_ln_na - (Za/Zi)*grad_ln_ni + (1 + (Za/Zi)*(C0a - 1))*grad_ln_Ti


  !---------------------------------------------------------------------------
  !-----------------------  Poloidal asymmetry  ------------------------------
  !---------------------------------------------------------------------------


  if (full_geom) then

    ! Full geometry calculation

    call fluxavg(nx, nth, theta, BV**2, jacob, B2avg)

    do ix = 1, nx
      b2(ix,:) = BV(ix,:)**2/B2avg(ix)
    enddo


    if (pol_asym) then

      ! Poloidally asymmetric case


      do i = 1, nx
        PhiV(i,:) = AsymPhi(i,1)*cos(theta) + AsymPhi(i,2)*sin(theta)
       	NV(i,:) = 1. + AsymN(i,1)*cos(theta) + AsymN(i,2)*sin(theta)
      enddo

      call asymmetry_fg(nx, nth, theta, BV, RV, jacob, FV, dpsidx, Machi, L11impi, &
                        R0, Ai, Aa, Zi, Za, B2avg, UU, GG, PhiV, NV, Te, Ti, regulopt, &
                        dmin, dmaj, nn)

    else

      ! Poloidally symmetric case

      dmin = 0.0_rkind
      dmaj = 0.0_rkind
      nn   = 1.0_rkind
      NV   = 1.0_rkind

    endif


    call fluxavg(nx, nth, theta, b2/nn, jacob, b2navg)
    call fluxavg(nx, nth, theta, nn/b2, jacob, nb2avg)
    call fluxavg(nx, nth, theta, nn/NV, jacob, nNVavg)
    call fluxavg(nx, nth, theta, b2/NV, jacob, b2NVavg)

    ! geometric coefficients in the equation for the flux
    Cgeo_G = nb2avg - 1._rkind/b2navg
    Cgeo_U = b2NVavg/b2navg - nNVavg

    Cgeo_Gcl = nb2avg

  else

    ! Simplified geometry

    B2avg = B0**2*(1.0 + 0.5*epsk**2)


    if (pol_asym) then

      !asymmetries of ES potential and main ion density
      dminphia = Za*(Te/Ti)*AsymPhi(:,1)
      dmajphia = Za*(Te/Ti)*AsymPhi(:,2)
      dNH = AsymN(:,1)
      dNV = AsymN(:,2)

      call asymmetry_an(nx, xn, UU, GG, epsk, invaspct, qmag, nuswca, deltaM, Ai, Aa, Zi, Za, &
                        dNH, dNV, dminphia, dmajphia, dmin, dmaj)

      do i = 1, nx
        nn(i,:) = 1.0 + dmin(i)*cos(theta) + dmaj(i)*sin(theta)
      enddo


    else

      ! Poloidally symmetric case

      dminphia = 0.0_rkind
      dmajphia = 0.0_rkind
      dNH = 0.0_rkind
      dNV = 0.0_rkind

      dmin = 0.0_rkind
      dmaj = 0.0_rkind
      nn   = 1.0_rkind

    endif


    dD2 = 0.5*(dmin**2 + dmaj**2)

    ! geometric coefficients in the equation for the flux

    Cgeo_G = 2.0*epsk*dmin + 2.0*eps2 + dD2
    Cgeo_U = epsk*(dNH - dmin) - dD2 + 0.5*(dmin*dNH + dmaj*dNV)
 
    Cgeo_Gcl = 1.0 + epsk*dmin + 2*eps2  

  endif


  !---------------------------------------------------------------------------
  !---------------------------  Impurity flux  -------------------------------
  !---------------------------------------------------------------------------

  ! General form: Va = -D*grad_ln_na + (K*grad_ln_ni + H*grad_ln_Ti + Vrot)
  !                  = -D*grad_ln_na + Vconv


  ! Pfirsch-Schlüter flux

  !Da_PS   = adps*ma*L11impi*FV**2*q_e*Ti*Cgeo_G*amin**2/(Za**2*q_e**2*B2avg*(dpsidx**2 + 1.e-33))
  Da_PS   = adps*qmag**2*rhoLimp2*L11impi*(Cgeo_G/(2.0*eps2))
  !Da_PS   = qmag**2*rhoLimp2*adps*L11impi*FV**2/(R0**2*B2avg)*(Cgeo_G/(2.0*eps2))
  Ka_PS   = (Za/Zi)*Da_PS
  Ha_PS   = -((1.0 + (Za/Zi)*(C0a - 1.0)) + (Cgeo_U/Cgeo_G)*(Za/Zi)*(C0a + ki))*Da_PS

  Vra_PS  = -Da_PS*grad_ln_na/amin + Ka_PS*grad_ln_ni/amin + Ha_PS*grad_ln_Ti/amin 
  Va_PS  = Ka_PS*grad_ln_ni/amin + Ha_PS*grad_ln_Ti/amin 


  ! Classical flux

  Da_CL   = (2.0*eps2*Cgeo_Gcl/Cgeo_G)*Da_PS/(2*qmag**2)
  Ka_CL   = (Za/Zi)*Da_CL
  Ha_CL   = -(1.0 + (Za/Zi)*(C0a - 1.0))*Da_CL

  Vra_CL  = -Da_CL*grad_ln_na/amin + Ka_CL*grad_ln_ni/amin + Ha_CL*grad_ln_Ti/amin
  Va_CL = Ka_CL*grad_ln_ni/amin + Ha_CL*grad_ln_Ti/amin


  ! Banana-Plateau flux

  call K_VISC(nx, Ni, Na, Ti, wii, wimpimp, Zi, Za, Ai, Aa, Tauii, &
              Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, y1, y2, y3, y4, &
              K11a, K12a, K22a, K11i, K12i, K22i)


  Da_BP = 1.5*q_e*Ti*(1.0/(1.0/K11a + 1.0/K11i))/(Za**2*q_e**2*FV**2*Na)
  Ka_BP = (Za/Zi)*Da_BP
  Ha_BP = ((Za/Zi)*(K12i/K11i - 1.5) - (K12a/K11a - 1.5))*Da_BP

  Vra_BP = -Da_BP*grad_ln_na/amin + Ka_BP*grad_ln_ni/amin + Ha_BP*grad_ln_Ti/amin
  Va_BP = Ka_BP*grad_ln_ni/amin + Ha_BP*grad_ln_Ti/amin
 

  ! Total transport coefficients

  Da = Da_PS + Da_BP + Da_CL ! total diffusion coefficient
  Ka = Ka_PS + Ka_BP + Ka_CL
  Ha = Ha_PS + Ha_BP + Ha_CL

  !Vra = Vra_PS + Vra_BP + Vra_CL

  ! total convective velocity
  Vconv = Ka*grad_ln_ni/amin + Ha*grad_ln_Ti/amin

  ! Total surface-averaged flux
  Flux_imp = -Da*gradNa + Na*Vconv


  return
end subroutine FACIT



!*******************************************************************************
!*********************** Complementary subroutines *****************************
!*******************************************************************************



subroutine fluxavg(nx,nth,thetay,AF,JJ,Aavg)

!*******************************************************************************
! Calculates the flux surface average (FSA) of a function AF(x,theta)
!*******************************************************************************
! INPUTS:
! -------
! nx -----> size of radial arrays [-] {int}
! nth ----> size of poloidal arrays [-] {int}
! thetay -> poloidal grid [-] {arr, nth}
! AF -----> function to average [-] {arr, (nx,nth)}
! JJ -----> Jacobian of the coordinate system [-] {arr, (nx,nth)}
!*******************************************************************************
! OUTPUT:
! ------
! Aavg ---> FSA of the AF function [-] {arr, nx}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), dimension(nx,nth) :: AF, JJ
  real(rkind), dimension(nth) :: thetay
  real(rkind), dimension(nx) :: Aavg
  real(rkind) :: denom
  integer   :: nx, nth, ix, ith

  Aavg = 0._rkind

  do ix = 1,nx
    denom = 0._rkind

    do ith=2,nth
      denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ix,ith) + JJ(ix,ith-1))
    enddo

    do ith = 2,nth
      Aavg(ix) = Aavg(ix)+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ix,ith)*JJ(ix,ith) + AF(ix,ith-1)*JJ(ix,ith-1))
    enddo

    if (denom.eq.0.) then
      Aavg(ix) = sum(AF(ix,:))/nth
    else
      Aavg(ix) = Aavg(ix)/denom
    endif

  enddo

  return
end subroutine fluxavg




subroutine fluxavgscal(nth,thetay,AF,JJ,Aavg)

!*******************************************************************************
! Calculates the flux surface average (FSA) of a function AF(x0, theta) at a
! given radial point x0
!*******************************************************************************
! INPUTS:
! -------
! nth ----> size of poloidal arrays [-] {int}
! thetay -> poloidal grid [-] {arr, nth}
! AF -----> function to average [-] {arr, nth}
! JJ -----> Jacobian of the coordinate system at x0 [-] {arr, nth}
!*******************************************************************************
! OUTPUT:
! ------
! Aavg ---> FSA of the AF function [-] {float}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), dimension(nth) :: AF, JJ
  real(rkind), dimension(nth) :: thetay
  real(rkind) :: denom, Aavg
  integer   :: nth, ith

  Aavg = 0._rkind
  denom = 0._rkind

  do ith=2,nth
    denom = denom + 0.5*(thetay(ith)-thetay(ith-1))*(JJ(ith) + JJ(ith-1))
  enddo

  do ith=2,nth
    Aavg = Aavg+ 0.5*(thetay(ith)-thetay(ith-1))*(AF(ith)*JJ(ith) + AF(ith-1)*JJ(ith-1))
  enddo

  if (abs(denom).gt.0.) then
    Aavg = Aavg/denom
  else
    Aavg = sum(AF)/nth
  endif

  return
end subroutine fluxavgscal



function ftrap(epsK)
!*******************************************************************************
! Trapped particle fraction as a function of the local inverse aspect ratio
!*******************************************************************************
! INPUT:
! ------
! epsK --> local inverse aspect ratio [-] {float}
!*******************************************************************************
! OUTPUT:
! ------
! ftrap -> trapped particle fraction [-] {float}
!*******************************************************************************
  use constants, only: rkind
  implicit none

  real(rkind) :: ftrap
  real(rkind), intent(inout) :: epsK

  ftrap = 1. - (1. - epsK)**2/(sqrt(1. - epsK**2)*(1+1.46*sqrt(epsK)))

  return
end function ftrap



subroutine asymmetry_an(nx, xn, UU, GG, epsk, invaspct, qmag, nuswca, deltaM, Ai, Aa, Zi, Za, &
                        dNH, dNV, dminphia, dmajphia, dmin, dmaj)

!*******************************************************************************
! Poloidal asymmetry of the impurity density distribution, analytical
! calculation in circular geometry
!*******************************************************************************
! INPUTS:
! -------
! nx -------> size of radial arrays [-] {int}
! xn -------> radial grid [-] {arr, nx}
! UU -------> thermodynamic gradient U [-] {arr, nx}
! GG -------> thermodynamic gradient G [-] {arr, nx}
! epsK -----> local inverse aspect ratio [-] {arr, nx}
! invaspct -> inverse aspect ratio [-] {float}
! qmag -----> safety factor [-] {arr, nx}
! nuswca ---> ratio of impurity coll. freq. to cyclotron freq. [-] {arr, nx}
! deltaM ---> rotation strength parameter [-] {arr, nx}
! Ai -------> main ion mass number [-] {float}
! Aa -------> impurity mass number [-] {float}
! Zi -------> main ion charge number [-] {float}
! Za -------> impurity charge number [-] {arr, nx}
! dNH ------> horizontal asymmetry of main ion density [-] {arr, nx}
! dNV ------> vertical asymmetry of main ion density [-] {arr, nx}
! dminphia -> horizontal asymmetry of electrostatic potential [-] {arr, nx}
! dmajphia -> vertical asymmetry of electrostatic potential [-] {arr, nx}
!*******************************************************************************
! OUTPUTS:
! --------
! dmin -----> horizontal asymmetry of impurity density [-] {arr, nx}
! dmaj -----> vertical asymmetry of impurity density [-] {arr, nx}
!*******************************************************************************

  use constants, only: rkind
  implicit none

  integer :: nx
  real(rkind), dimension(nx) :: UU, GG, epsK, qmag, nuswca, deltaM, Za, dNH, dNV
  real(rkind), dimension(nx) :: dminphia, dmajphia, xn
  real(rkind) :: Ai, Aa, Zi, invaspct
  real(rkind), dimension(nx) :: dmin, dmaj
  real(rkind), dimension(nx) :: RR, UG, Ae, AGe, CD0, HH, QQ, FF, KK, CD, CDV
  real(rkind), dimension(nx) :: RD, DD, num, cosa, sina



  !RR = 0.5*(1-(Ai*Za)/(Aa*Zi))/(epsk + 1.e-33)
  !RR = 0.5*(1-Ai*Za/(Aa*Zi))/(xn+1.e-33)
  RR = 0.0_rkind

  UG = 1 + UU/GG
  !Ae = nuswca*qmag**2/(epsk + 1.e-33)
  Ae = nuswca*qmag**2/invaspct
  AGe = Ae*GG
  CD0 = -(epsk + 1.e-33)/UG
  !HH = 1.0 + deltaM*CD0*RR/GG
  HH = 1.0_rkind
  !QQ = CD0*(dNV/(epsk + 1.e-33))*UU/GG
  QQ = CD0*(dNV/(epsk + 1.e-33))*(UG-1.0)
  !FF = CD0*(1-0.5*dNH/(epsk + 1.e-33)*UU/GG - deltaM*(RR/GG)/(epsk + 1.e-33))
  !FF = CD0*(1-0.5*dNH*(UG-1.0)/(epsk + 1.e-33) - deltaM*(RR/GG)/(epsk + 1.e-33))
  FF = CD0*(1-0.5*dNH*(UG-1.0)/(epsk + 1.e-33) )
  KK = 1.0_rkind

  CD = FF - 0.5*(dminphia - deltaM)
  CDV = -0.5*(dmajphia + QQ)
  RD = sqrt((FF + 0.5*(dminphia-deltaM))**2 + 0.25*(dmajphia - QQ)**2)
  DD = RD**2 + AGe**2*(RD/CD0)**2

  num  = ((AGe/CD0)**2 - 1)*(FF/(CD0) + 0.5*(dminphia-deltaM)/CD0) + &
         (AGe/CD0)*(0.5*dNV*(UG-1.0)/(epsk + 1.e-33) - 0.5*dmajphia/CD0)
  cosa = RD*CD0*num/DD

  num  = 2*AGe*(FF/CD0 + 0.5*(dminphia-deltaM)/CD0)+((AGe/CD0)**2-1)*&
         (0.5*dmajphia - 0.5*dNV*CD0*(UG-1.0)/((epsk + 1.e-33)))
  sina = RD*num/DD

  dmin = CD + RD*cosa
  dmaj = CDV + RD*sina

end subroutine asymmetry_an



subroutine asymmetry_fg(nx, nth, theta, BV, RV, jacob, FV, dpsidx, Machi, L11impi, &
                        R0, Ai, Aa, Zi, Za, B2avg, UU, GG, PhiV, NV, Te, Ti, regulopt, &
                        dmin, dmaj, nn)
!*******************************************************************************
! Poloidal asymmetry of the impurity density distribution, iterative
! calculation in full geometry
!*******************************************************************************
! INPUTS:
! -------
! nx -------> size of radial arrays [-] {int}
! nth ------> size of poloidal arrays [-] {int}
! theta ----> poloidal grid [-] {arr, nth}
! BV -------> magnetic field [T] {arr, (nx,nth)}
! RV -------> major radius [m] {arr, (nx,nth)}
! jacob ----> Jacobian of the coordinate system [m/T] {arr, (nx,nth)}
! FV -------> poloidal current flux function [T*m] {arr, nx}
! dpsidx ---> radial derivative of poloidal flux [V*s/-] {arr, nx}
! Machi ----> Mach number of main ion [-] {arr, nx}
! L11impi --> impurity - main ion collision frequency [1/s] {arr, nx}
! R0 -------> major radius at magnetic axis [m] {float}
! Ai -------> main ion mass number [-] {float}
! Aa -------> impurity mass number [-] {float}
! Zi -------> main ion charge number [-] {float}
! Za -------> impurity charge number [-] {arr, nx}
! B2avg ----> FSA of magnetic field squared <BV**2> [T^2] {arr, nx}
! UU -------> thermodynamic gradient U [-] {arr, nx}
! GG -------> thermodynamic gradient G [-] {arr, nx}
! PhiV -----> poloidal asymmetry of electrostatic potential [-] {arr, (nx,nth)}
! NV -------> poloidal asymmetry of main ion density [-] {arr, (nx,nth)}
! Te -------> electron temperature [eV] {arr, nx}
! Ti -------> main ion temperature [eV] {arr, nx}
! regulopt -> options for iterative calculations [-] {arr, 4}
!*******************************************************************************
! OUTPUTS:
! --------
! dmin -----> horizontal asymmetry of impurity density [-] {arr, nx}
! dmaj -----> vertical asymmetry of impurity density [-] {arr, nx}
! nn -------> poloidal asymmetry of the impurity density Na/<Na> [-] {arr, (nx,nth)}
!*******************************************************************************

  use constants, only: rkind, mp, q_e, pi
  implicit none

  integer :: nx, nth
  real(rkind), dimension(nth) :: theta
  real(rkind), dimension(nx) :: FV, dpsidx, Machi, L11impi, Za, UU, GG, Te, Ti, B2avg
  real(rkind) :: R0, Ai, Aa, Zi
  real(rkind), dimension(nx, nth) :: BV, RV, jacob, PhiV, NV
  real(rkind), dimension(4) :: regulopt

  real(rkind), dimension(nx) :: dmin, dmaj
  real(rkind), dimension(nx, nth) :: nn


  real(rkind), dimension(nx) :: asym_error, Factrot0
  real(rkind), dimension(1:nth+1,1:nth+1) :: AAA, LL, CC, DD
  real(rkind), dimension(1:nth+1) :: BB, nnya
  real(rkind), dimension(1:nth) :: Apsi, dum, nnx, nny, b2, nnp, TermGG, TermUU, FFF, GGG, HHH
  real(rkind), dimension(size(AAA,1)) :: work  ! work array for LAPACK
  real(rkind)  :: b2snavg, nsb2avg, b2sNNavg, FactV
  !real(rkind)  :: TermUUneo, TermRRneo, TermRRneoa, TermRRneob, TermGGcl
  !real(rkind)  :: TermGGneo, TermUUneoa, TermUUneob, TermRRcla, TermRRclb, TermRRcl
  real(rkind)  :: err, prog, regulweight, dtheta, Erreur, dumscal
  integer, dimension(size(AAA,1)) :: ipiv   ! pivot indices
  integer :: n, info, ix, it, ierr, ierrmax



  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  !allocate(AA(nth+1,nth+1),LL(nth+1,nth+1),CC(nth+1,nth+1),DD(nth+1,nth+1),BB(nth))

  err         = regulopt(1)
  prog        = regulopt(2)
  regulweight = regulopt(3)
  ierrmax     = int(regulopt(4))


  ! Poloidal asymmetry
  
  nn = 1.0_rkind

  do ix = 1,nx

 	Factrot0(ix) = (Aa/Ai)*(Machi(ix)/R0)**2

	do it=1,nth
	    b2(it) = BV(ix,it)**2/B2avg(ix)
	    Apsi(it) = jacob(ix,it)*FV(ix)*(Aa*mp)*L11impi(ix)/(q_e*Za(ix))/((dpsidx(ix))**2+1.e-33)
	enddo

   	call fluxavgscal(nth, theta , b2/NV(ix,:), jacob(ix,:), b2sNNavg)

  	nnp = 2
  	nnx = nn(ix,:)
  	ierr = 1
  	Erreur = 2*err

  	do while ( Erreur>err .and. ierr<ierrmax)

      AAA = 0._rkind
  		LL  = 0._rkind
  		BB  = 0._rkind

  		call fluxavgscal(nth, theta, b2/nnx, jacob(ix,:), b2snavg)

  		TermGG = 1. - b2/nnx/b2snavg
  		TermUU = b2/NV(ix,:) - b2sNNavg*b2/nnx/b2snavg
 
  		FFF  = Apsi*( GG(ix) + b2/NV(ix,:)*UU(ix) )
  		GGG = -Za(ix)*Te(ix)/Ti(ix)*(PhiV(ix,:)-PhiV(ix,1)) + Factrot0(ix)*(RV(ix,:)**2-RV(ix,1)**2)
  		HHH = Apsi*b2/b2snavg*(GG(ix) + b2sNNavg*UU(ix) )

      ! AA*n = BB
      ! Remplissage matrices AA, BB & LL

  		do it = 2,nth-1

  			dtheta = 0.5*(theta(it+1) - theta(it-1))
  			AAA(it,it-1) = -0.5/dtheta
  			AAA(it,it) = -FFF(it)-0.5/dtheta*(GGG(it+1)-GGG(it-1))
  			AAA(it,it+1) = 0.5/dtheta

  			LL(it, it-1) = 1.
  			LL(it, it) = -2.
  			LL(it, it+1) = 1.

  			BB(it) = -HHH(it)

  		enddo

  		it = 1
  		dtheta = 0.5*(theta(2)-theta(nth) + 2*PI)
  		AAA(it,it)   = -FFF(it) - (0.5/dtheta)*(GGG(it+1) - GGG(nth-1))
  		AAA(it,it+1) = 0.5/dtheta
  		AAA(it,nth)  = -0.5/dtheta
  		LL(it, it)   = -2.
  		LL(it, it+1) = 1.
  		LL(it, nth)  = 1.
  		BB(it) = -HHH(it)


  		it = nth
  		dtheta = 0.5*(theta(1)-theta(nth-1) + 2*PI)
  		AAA(it,it-1) = -0.5/dtheta
  		AAA(it,it)   = -FFF(it)-0.5/dtheta*(GGG(1)-GGG(it-1))
  		AAA(it,it+1) = 0.5/dtheta
  		LL(it, it-1) = 1.
  		LL(it, it)   = -2.
  		LL(it, it+1) = 1.
  		BB(it) = -HHH(it)

  		it = nth+1
  		dtheta = 0.5*(theta(2)-theta(nth) + 2*PI)
      		AAA(it,it) = -1.
  		AAA(it,1)  = 1
  		BB(it) = 0.

  		! CC = transpose(AA)*AA+regulweight*transpose(LL)*LL
  		CC = transpose(AAA)
  		CC = matmul(CC,AAA)
  		DD = transpose(LL)
  		DD = matmul(DD,LL)
  		CC = CC + regulweight*DD

  		n  = size(CC,1)
  		DD = CC
  		call DGETRF(n, n, DD, n, ipiv, info)

   		if (info /= 0) then
			   stop 'Matrix is numerically singular!'
  		end if

  		! DGETRI computes the inverse of a matrix using the LU factorization
  		! computed by DGETRF.
   		call DGETRI(n, DD, n, ipiv, work, n, info)

      if (info /= 0) then
			   stop 'Matrix inversion failed!'
  		end if

      CC = matmul(DD,transpose(AAA))
  		nnya = matmul(CC,BB)
      !nnya = 1.

  		call fluxavgscal(nth, theta, nnya(1:nth), jacob(ix,:), dumscal)

  		nny = nnya(1:nth)/dumscal
  		nny = max(nny, 1.e-5)

  		nnp = nnx
  		nnx = prog*nnp + (1-prog)*nny

  		do it = 2, nth-1
			   dum(it) = (log(nnx(it+1))-log(nnx(it-1))-GGG(it+1)+GGG(it-1))/(theta(it+1)-theta(it-1))
  		enddo

  		it = 1
  		dum(it) = (log(nnx(it+1))-log(nnx(it))-GGG(it+1)+GGG(it))/(theta(it+1)-theta(it))

  		it = nth
  		dum(it) = (log(nnx(it))-log(nnx(it-1))-GGG(it)+GGG(it-1))/(theta(it)-theta(it-1))

  		Erreur = maxval(abs(dum-FFF+HHH/nnx))

  		ierr = ierr+1

    end do

  	asym_error(ix) = Erreur

  	nn(ix,:) = nnx
  	dmin(ix) = 2.*sum((nnx-1.)*cos(theta))/float(nth)
  	dmaj(ix) = 2.*sum((nnx-1.)*sin(theta))/float(nth)

  enddo

end subroutine asymmetry_fg







function ki_Redl(nuistar, ft, Zeff)

!*******************************************************************************
! Neoclassical ion flow coefficient, fitted w.r.t. NEO in Redl PoP (2021)
!*******************************************************************************
! INPUTS:
! -------
! nuistar -> main ion collisionality [-]
! ft ------> trapped particle fraction [-]
! Zeff ----> effective charge [-]
!*******************************************************************************
! OUTPUT:
! ------
! ki ------> main ion flow coefficient [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), intent(in) :: nuistar, ft, Zeff
  real(rkind) :: alpha0
  real(rkind) :: ki_Redl

  alpha0 = -(0.62 + 0.055*(Zeff - 1.0))*(1.0 - ft)/((0.53 + 0.17*(Zeff - 1.0))*&
            (1.0 - (0.31 - 0.065*(Zeff - 1.0))*ft - 0.25*ft**2))

  ki_Redl = ((alpha0 + 0.7*Zeff*SQRT(ft*nuistar))/(1.0 + 0.18*SQRT(nuistar))-&
              0.002*nuistar**2*ft**6)/(1.0 + 0.004*nuistar**2*ft**6)

  return
end function ki_Redl




function C2(alpha, g, f1, f2, Aimp, Ai)

!*******************************************************************************
! new parametrized version of the C2 function
! see Hirshman-Sigmar NF (1981), page 61 (1138)
!*******************************************************************************
! INPUTS:
! -------
! alpha -> impurity strength parameter [-]
! g -----> main ion collisionality parameter [-]
! f1 ----> f1(Zimp) factor [-]
! f2 ----> f2(Zimp) factor [-]
!*******************************************************************************
! OUTPUT:
! ------
! C2 ----> term in impurity-ion friction C0a [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none

  real(rkind), intent(in) :: alpha, g, f1, f2, Aimp, Ai
  real(rkind) :: C2

  C2 = 1.5/(1.0 + (Ai/Aimp)*f1) - (0.29 + 0.68*alpha)/(0.59 + alpha + (1.34 + f2)/g**2)

  return
end function C2



subroutine facs(nx, Zimp, ft, f1, f2, f3, y1, y2, y3, y4, adps)

!*******************************************************************************
! set of factors fitted with respect to NEO
!*******************************************************************************
! INPUTS:
! -------
! nx ---> size of radial arrays [-]
! Zimp -> impurity charge [-]
! ft ---> trapped particle fraction [-]
!*******************************************************************************
! OUTPUTS:
! --------
! f1 --> low collisionality saturation of PS C0a coefficient [-]
! f2 --> correction to low-Z impurity-ion friction in C0a coefficient [-]
! f3 --> factor in i-e heat exchange contribution to C0a [-]
! y1 --> factor of the (1,1) banana viscosity coefficient of the impurity [-]
! y2 --> factor of the (1,2) plateau viscosity coefficient of the impurity [-]
! y3 --> factor of the (1,2) Pfirsch-Schlüter viscosity coefficient of the impurity [-]
! y4 --> factor of the (1,2) banana viscosity coefficient of the main ion [-]
!*******************************************************************************

  use constants, only: rkind
  implicit none


  integer :: nx
  !real(dp) :: Zimp
  real(rkind), dimension(nx) :: Zimp
  real(rkind), dimension(nx) :: ft

  !real(dp) :: f1, f2, f3
  real(rkind), dimension(nx) :: f1, f2, f3
  real(rkind), dimension(nx) :: adps, y2, y3
  real(rkind), dimension(nx) :: w11, w12, w13, w14, y1
  real(rkind), dimension(nx) :: wp1, wp2, wp3, wp4, yp, y4


  ! f's:
  f1 =  (-6.83808564e5 + 2.46855534e6*Zimp)/( 1 + 6.04692708e5*Zimp**1.61470425)
  f2 = (88.28389935 + 10.50852772*Zimp)/( 1 + 0.2157175*Zimp**2.57338463)
  f3 = (-3.15171054e6 + 1.92908543e6*Zimp)/(1+5.26716920e6*Zimp**8.33610108e-01)


  ! y's:
  w11 = 1.23805214e-05*Zimp**3 - 1.03611576e-03*Zimp**2 + 1.85221287e-02*Zimp + 1.29758029
  w12 = -5.84996779e-05*Zimp**3 + 4.79623045e-03*Zimp**2 -1.01924030e-01*Zimp - 5.81816797
  w13 = 5.98152997e-05*Zimp**3 - 4.84255587e-03*Zimp**2 + 1.03585208e-01*Zimp + 5.71139227
  w14 = -1.33253954e-05*Zimp**3 + 1.04416741e-03*Zimp**2 - 1.94973421e-02*Zimp - 1.10061492

  y1 = w11*ft**2 + w12*ft + w13*ft**0.5 + w14

  wp1 = ( 0.11603574 + 0.47297835*Zimp**0.94456671)/( 1 + 0.1245426*Zimp**1.20221441)
  wp2 = ( 0.6327602 - 2.92116611*Zimp**1.06310182)/( 1 + 0.30469549*Zimp**1.16495283)
  wp3 = ( -0.69318477 + 2.85619511*Zimp**1.11765611)/( 1 + 0.34127361*Zimp**1.1952383)
  wp4 = ( 1.80217558 - 1.0080554*Zimp**0.96345934)/( 1 + 0.51677339*Zimp**1.11131896)

  yp = wp1*ft**2 + wp2*ft + wp3*ft**0.5 + wp4
  y4 = yp/y1

  y2 = 21.31*ft**3 - 21.88*ft**2 + 7.316*ft + 0.6264
!  y3 = ((4.2857e5 - 4.4978e5*Zimp)/(1 - 1.3557e5*Zimp**1.9))*((-9.09204093e6 + 8.15802759e6*Zimp)/(1.0 + 3.25263641e6*Zimp**1.07640221))
  y3 = 8.85/Zimp**2.98 - 7.96/Zimp**1.82 - 9.27/Zimp**1.98 + 8.34/Zimp**0.82

!  adps = (-110378.3491 + 753838.926571*Zimp**1.06107833841)/(1+ 1007273.22737*Zimp)
  adps = -0.109/Zimp + 0.743*Zimp**0.06


end subroutine facs


subroutine K_VISC(nx, ni, nimp, Ti, wii, wimpimp, Zi, Zimp, Ai, Aimp, Tauii, &
                  Tauimpi, Tauiimp, Tauimpimp, eps2, R0, qmag, ft, y1, y2, y3, y4, &
                  K11a, K12a, K22a, K11i, K12i, K22i)

!*******************************************************************************
! Calculates the positive definite matrix of neoclassical viscosity coefficients
! in a fully analytical way by solving for the coefficients in the individual
! banana, plateau and Pfirsch-Schlüter collisionality regimes and using a
! rational approximation to interpolate between them, as well as evaluating
! kinetic integrals for a Maxwellian distribution
!*******************************************************************************
! INPUTS:
! -------
! nx --------> size of radial arrays [-]
! ni --------> main ion density [1/m^3]
! nimp ------> impurity density [1/m^3]
! Ti --------> main ion temperature [eV]
! wii -------> ion transit frequency [1/s]
! wimpimp ---> impurity transit frequecy [1/s]
! Zi --------> main ion charge [-]
! Zimp ------> impurity charge [-]
! Ai --------> main ion mass [-]
! Aimp ------> impurity mass [-]
! Tauii -----> ion-ion collision time [s]
! Tauimpi ---> impurity-ion collision time [s]
! Tauiimp ---> ion-impurity collision time [s]
! Tauimpimp -> impurity impurity collision time [s]
! eps2 ------> local inverse aspect ratio squared [-]
! R0 --------> major radius at magnetic axis [m]
! qmag ------> safety factor [-]
! ft --------> trapped particle fraction [-]
! y1 --------> fitted factor from facs subroutine [-]
! y2 --------> fitted factor from facs subroutine [-]
! y3 --------> fitted factor from facs subroutine [-]
! y4 --------> fitted factor from facs subroutine [-]
!*******************************************************************************
! OUTPUTS:
! --------
! K11a ------> (1,1) impurity viscosity coefficient [kg/(m*s)]
! K12a ------> (1,2)=(2,1) impurity viscosity coefficient [kg/(m*s)]
! K22a ------> (2,2) impurity viscosity coefficient [kg/(m*s)]
! K11i ------> (1,1) main ion viscosity coefficient [kg/(m*s)]
! K12i ------> (1,2)=(2,1) main ion viscosity coefficient [kg/(m*s)]
! K22i ------> (2,2) main ion viscosity coefficient [kg/(m*s)]
!*******************************************************************************

  use constants, only: rkind, mp, q_e, pi, sqrt2
  implicit none

  integer :: nx
  real(rkind), dimension(nx) :: ni, nimp, Ti, wii, wimpimp, Tauii, Tauimpi, Tauiimp, Tauimpimp, eps2, qmag, ft
  real(rkind) :: Zi, Ai, Aimp, R0, mimp, mi
  !real(dp) :: Zimp
  real(rkind), dimension(nx) :: Zimp
  real(rkind), dimension(nx) :: K11a, K12a, K22a, K11i, K12i, K22i
  !real(dp) :: f1, f2, f3
  real(rkind), dimension(nx) :: y1, y2, y3, y4

  real(rkind), dimension(nx) :: fac_a_P, fac_i_P, K11aP, K12aP, K22aP, K11iP, K12iP, K22iP
  real(rkind) :: r00, r01, r11, xai, xia, x2ai, x2ia, xfac_ai, xfac_ia
  real(rkind) :: qaa00, qaa01, qaa11, qai00, qia00, qai01, qia01, qai11, qia11
  real(rkind), dimension(nx) :: fac_qai_PS, fac_qia_PS
  real(rkind), dimension(nx) :: qa00, qi00, qa01, qi01, qa11, qi11, Qa, Qi
  real(rkind), dimension(nx) :: la11, la12, la22, li11, li12, li22
  real(rkind), dimension(nx) :: fac_imp_PS, fac_ion_PS
  real(rkind), dimension(nx) :: K11aPS, K12aPS, K22aPS, K11iPS, K12iPS, K22iPS
  real(rkind), dimension(nx) :: fac_B, nuDai, nuD2ai, nuD4ai
  real(rkind), dimension(nx) :: nuDia, nuD2ia, nuD4ia
  real(rkind), dimension(nx) :: nuDaa, nuD2aa, nuD4aa
  real(rkind), dimension(nx) :: nuDii, nuD2ii, nuD4ii
  real(rkind), dimension(nx) :: K11aB, K12aB, K22aB, K11iB, K12iB, K22iB


  mimp = Aimp*mp
  mi = Ai*mp


  ! Plateau regime

  fac_a_P = nimp*(q_e*Ti)*SQRT(pi)/(3.0*wimpimp)
  fac_i_P = ni*(q_e*Ti)*SQRT(pi)/(3.0*wii)

  K11aP = fac_a_P*2.0
  K12aP = fac_a_P*2.0*3.0
  K22aP = fac_a_P*2.0*3.0*4.0

  K11iP = fac_i_P*2.0
  K12iP = fac_i_P*2.0*3.0
  K22iP = fac_i_P*2.0*3.0*4.0

  ! Pfirsch-Schlüter regime

  r00 = 1.0/sqrt2
  r01 = 1.5/sqrt2
  r11 = 3.75/sqrt2

  xai = SQRT(Aimp/Ai)
  xia = 1.0/xai

  x2ai = xai**2
  x2ia = xia**2

  xfac_ai = (1+x2ai)**0.5
  xfac_ia = (1+x2ia)**0.5

  qaa00 = 8.0/2**1.5
  qaa01 = 15.0/2**2.5
  qaa11 = 132.5/2**3.5

  qai00 = (3.0+5.0*x2ai)/xfac_ai**3
  qia00 = (3.0+5.0*x2ia)/xfac_ia**3

  qai01 = 1.5*(3.0+7.0*x2ai)/xfac_ai**5
  qia01 = 1.5*(3.0+7.0*x2ia)/xfac_ia**5

  qai11 = (35.0*x2ai**3 + 38.5*x2ai**2 + 46.25*x2ai + 12.75)/xfac_ai**7
  qia11 = (35.0*x2ia**3 + 38.5*x2ia**2 + 46.25*x2ia + 12.75)/xfac_ia**7

  fac_qai_PS = (ni*Zi**2/(nimp*Zimp**2))
  fac_qia_PS = (nimp*Zimp**2/(ni*Zi**2))


  qa00 = fac_qai_PS*qai00 +  qaa00 - r00
  qi00 = fac_qia_PS*qia00 + qaa00 - r00

  qa01 = fac_qai_PS*qai01 + qaa01 - r01
  qi01 = fac_qia_PS*qia01 + qaa01 - r01

  qa11 = fac_qai_PS*qai11 + qaa11 - r11
  qi11 = fac_qia_PS*qia11 + qaa11 - r11


  Qa = 0.4*(qa00*qa11-qa01*qa01)
  Qi = 0.4*(qi00*qi11-qi01*qi01)


  la11 = qa11/Qa
  la12 = 3.5*(qa11+qa01)/Qa
  la22 = 12.25*(qa11+qa00+2*qa01)/Qa

  li11 = qi11/Qi
  li12 = 3.5*(qi11 + qi01)/Qi
  li22 = 12.25*(qi11+qi00 + 2.0*qi01)/Qi

  fac_imp_PS = nimp*(q_e*Ti)*Tauimpimp
  fac_ion_PS = ni*(q_e*Ti)*Tauii


  K11aPS = fac_imp_PS*la11
  K12aPS = fac_imp_PS*la12
  K22aPS = fac_imp_PS*la22

  K11iPS = fac_ion_PS*li11
  K12iPS = fac_ion_PS*li12
  K22iPS = fac_ion_PS*li22



  ! Banana regime

  fac_B = (ft/(1.0 - ft))*(2.0*R0**2*qmag**2/(3.0*eps2))

  ! Maxwellian integrals

  nuDai  = (xfac_ai + x2ai*LOG(xai/(1 + xfac_ai)))/Tauimpi
  nuD2ai = 1.0/(xfac_ai*Tauimpi)
  nuD4ai = 2.0*(1.0 + 1.25*x2ai)/(xfac_ai**3*Tauimpi)

  nuDia  = (xfac_ia + x2ia*LOG(xia/(1 + xfac_ia)))/Tauiimp
  nuD2ia = 1.0/(xfac_ia*Tauiimp)
  nuD4ia = 2.0*(1.0 + 1.25*x2ia)/(xfac_ia**3*Tauiimp)

  nuDaa  = (sqrt2 + LOG(1/(1 + sqrt2)))/Tauimpimp
  nuD2aa = 1.0/(sqrt2*Tauimpimp)
  nuD4aa = 4.5/(2.0**1.5*Tauimpimp)

  nuDii  = (sqrt2 + LOG(1/(1 + sqrt2)))/Tauii
  nuD2ii = 1.0/(sqrt2*Tauii)
  nuD4ii = 4.5/(2.0**1.5*Tauii)

  ! Banana regime viscosity coefficient

  K11aB = fac_B*nimp*mimp*(nuDai + nuDaa)
  K12aB = fac_B*nimp*mimp*(nuD2ai + nuD2aa)
  K22aB = fac_B*nimp*mimp*(nuD4ai + nuD4aa)

  K11iB = fac_B*ni*mi*(nuDia + nuDii)
  K12iB = fac_B*ni*mi*(nuD2ia + nuD2ii)
  K22iB = fac_B*ni*mi*(nuD4ia + nuD4ii)

  ! total viscosity coefficients, rational approximation for interpolation

  K11a = y1*K11aB/((1.0 + y1*K11aB/(K11aP))*(1.0 + K11aP/(K11aPS)))
  K12a = K12aB/((1.0 + K12aB/(y2*K12aP))*(1.0 + y2*K12aP/(y3*K12aPS)))
  K22a = K22aB/((1.0 + K22aB/(K22aP))*(1.0 + K22aP/(K22aPS)))

  K11i = K11iB/((1.0 + K11iB/(K11iP))*(1.0 + K11iP/(K11iPS)))
  K12i = y4*K12iB/((1.0 + y4*K12iB/(K12iP))*(1.0 + K12iP/(K12iPS)))
  K22i = K22iB/((1.0 + K22iB/(K22iP))*(1.0 + K22iP/(K22iPS)))



  return
end subroutine K_VISC
