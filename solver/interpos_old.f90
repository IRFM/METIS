







#ifdef MATLAB73 
#include "fintrf.h"
#endif

MODULE prec_rkind
!
!   Precision for real and complex
!
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,300)
  INTEGER, PARAMETER :: CKIND = RKIND
!
END MODULE prec_rkind
MODULE interpos_module
  !
  ! Examples of calls to this generic routine:
  !
  ! call interpos(xin,yin,nin,xout,yout,nout)
  ! call interpos(xin,yin,nin,xout,yout,nout,tension)
  ! call interpos(xin,yin,nin,xout,yout,nout,tension,youtp=yp)
  ! call interpos(xin,yin,nin,xout,yout,nout,tension,nbc=bc_type,ybc=bc_values)
  ! call interpos(xin,yin,nin,xout,yout,nout,tension,nbc=bc_type,ybc=bc_values,option=32)
  !
  ! The full possibilities can be expressed as follows:
  !
  ! call interpos(xin,yin[[,yinnew],yinpp,]nin,xout,yout,nout,[tension,youtp,youtpp,youtint,nbc,ybc,sigma,option])
  ! with nout optional if equals one (xout, yout, etc are scalars)
  !
  ! These are grouped into 3 main routines:
  !
  ! main:
  ! call interpos(xin,yin,nin,[xout,yout,nout,tension,youtp,youtpp,youtint,nbc,ybc,sigma,option])
  ! These two to compute 1st yinpp and then to just evaluate for different xout values with second calls
  ! call interpos(xin,yin,yinnew,yinpp,nin,[tension,nbc,ybc,sigma])                 ! only computes yinpp
  ! call interpos(xin,yinnew,yinpp,nin,nout,xout,yout,youtp,youtpp,youtint,option])    ! only computes yout's
  !
  ! Each of these allow for scalar xout, yout, etc which means that nout is also optional in this case
  !
  !
  implicit none
  public :: interpos
  interface interpos
     module procedure interpos_def, interpos_defper, interpos_yinpp, interpos_interp &
       & ,interpos_defxscal, interpos_defperxscal, interpos_interpxscal
  end interface

contains

  SUBROUTINE interpos_def(XIN,YIN,NIN,NOUT,tension,xout,yout,youtp,youtpp,youtint,nbc,ybc,sigma,option,info)
    !
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: NIN
    INTEGER, intent(in),optional :: NOUT
    REAL(RKIND), intent(in) ::  XIN(NIN), YIN(NIN)
    REAL(RKIND), intent(in), optional ::  tension
    REAL(RKIND), intent(in), optional ::  xout(:)
    REAL(RKIND), intent(out), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    REAL(RKIND), intent(in), optional ::  ybc(6), sigma(:)
    INTEGER, intent(in), optional :: nbc(2), option
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:)
    REAL(RKIND) ::  sigma2(NIN)
    integer :: info2, nout2
    !
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, inttype, optabs
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    !write(*,*) 'debut interpos'
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin
!      if (present(nout)) nout = nin
      allocate(xout2(nout2))
      xout2 = xin(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
       if (tension .ge. 0._rkind) then
          sigma2=tension*sigma2
       else
          ! tension=-1 or negative, uses default value |tension| * Delta_x**3
          zz=minval(xin(2:nin)-xin(1:nin-1))
          ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
          sigma2=abs(tension)*zz**3*sigma2
       end if
    else
       zz=minval(xin(2:nin)-xin(1:nin-1))
       sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbc given, ybc should also be given)
    if (present(nbc)) then
      if (present(ybc)) then
        NBCLFT = nbc(1)
        NBCRGT = nbc(2)
        YBCLFT = ybc(1)
        YBCRGT = ybc(2)
        if (NBCLFT .ge. 10) XBCLFT = ybc(3)
        if (NBCRGT .ge. 10) XBCRGT = ybc(4)
        PXEXP0 = ybc(5)
        PXEXPDL= ybc(6)
      else
        print *,' problems in interpos, nbc is given but not ybc'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = xin(1)
      XBCRGT = xin(nin)
      PXEXP0 = xin(1)-xin(nin)
      PXEXPDL= 1._rkind
    end if
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      inttype=3 ! cubic
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        inttype = 1
        ! extrapolation
        if (optabs .eq. 1)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 11) iextrapo = sign(1,option) ! linear extrapolation with edge-deriv for 11; edge-val for -11
        if (optabs .eq. 21) iextrapo = sign(10,option) ! constant outside: y=yedge for option=21, y=0. for -21
      case(2,12,22,32,42)
        ! quadratic interpolation
        inttype = 2
        ! extrapolation
        if (optabs .eq. 2)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 12) iextrapo = sign(21,option) ! quadratic on dx then linear extrapol. with edge-deriv for 12; edge-val for -12
        if (optabs .eq. 22) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 32) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 32; edge-val for -32
        if (optabs .eq. 42) iextrapo = sign(10,option) ! constant outside: y=yedge for option=42, y=0. for -42
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          iextrapo = 1
        case (2)
          iextrapo = 21
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    !write(*,*) 'type interpolation ',inttype
    select case (inttype)
    case (1)
      ! linear
      call intlinear(XIN,YIN,NIN,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo)
    case (2)
      ! quadratic
      call intquadratic(XIN,YIN,NIN,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo,NBCLFT)
    case (3)
      ! cubic
      ! general cubic spline routine with non-periodic boundary conditions

      CALL CBSPLGEN(XIN,YIN,YNEW,YINPP,Nin,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
        &    nout2,ioptder,sigma2,NBCLFT, &
        &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL,info2)
      !
    case default
      print *,' inttype=',inttype,' ; error in this routine, this should not be possible'
      if (present(info)) info = 4
      return
    end select
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_def

  SUBROUTINE interpos_defper(XIN,YIN,NIN,NOUT,tension,xout,yout,youtp,youtpp,youtint,NBC,YBC,sigma,info)
    ! for periodic boundary conditions, nbc=-1 and ybc=period should be given
    ! option not needed since no extrapolation and only valid for cubic spline here
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: NIN, NBC
    INTEGER, intent(in),optional :: NOUT
    REAL(RKIND), intent(in) ::  XIN(NIN), YIN(NIN)
    REAL(RKIND), intent(in), optional ::  tension
    REAL(RKIND), intent(in), optional ::  xout(:)
    REAL(RKIND), intent(out), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    REAL(RKIND), intent(in), optional ::   ybc, sigma(:)
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:)
    REAL(RKIND) ::  sigma2(NIN)
    integer :: info2, nout2, nin_eff
    !
    INTEGER ioptder
    REAL(RKIND) :: period
    REAL(RKIND) ::  ynew(NIN), yinpp(NIN), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! period
    nin_eff = nin
    if (present(ybc)) then
      period = ybc
      ! check if extra point is given, if yes skip it with warning:
      if (abs(xin(nin_eff)-xin(1)-period) .lt. 1.0e-10_rkind*(xin(2)-xin(1))) then
        print *,'% it seems periodic point x(n)=x(1)+period is included in input, use nin-1 points'
        nin_eff = nin_eff-1
      end if
    else
      ! Assumes last x point corresponds to x(1)+period and is redundant
      period=xin(nin_eff)-xin(1)
      nin_eff = nin_eff-1
    end if
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin_eff
!      if (present(nout)) nout = nin_eff
      allocate(xout2(nout2))
      xout2 = xin(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      if (tension .ge. 0._rkind) then
        sigma2=tension*sigma2
      else
        ! tension=-1 or negative, uses default value |tension| * Delta_x**3
        zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
        ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
        sigma2=abs(tension)*zz**3*sigma2
      end if
    else
      zz=minval(xin(2:nin_eff)-xin(1:nin_eff-1))
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions
    if (nbc .ne. -1) then
      print *,'nbc when scalar should be -1 for periodic boundary conditions, otherwise vector of size 2'
      if (present(info)) info = 2
      return
    end if

    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic
    ! general cubic spline routine with periodic boundary conditions
    CALL CBSPLGNP(XIN,YIN,YNEW,YINPP,Nin_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
      &    nout2,ioptder,sigma2,period,PXEXP0,PXEXPDL,info2)
      !
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    ! if xout = xin, could add extra periodic point if nin_eff was set to nin-1
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_defper

  SUBROUTINE interpos_yinpp(XIN,YIN,YINNEW,YINPP,NIN,tension,nbc,ybc,sigma,info)
    !
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: NIN
    REAL(RKIND), intent(in) ::  XIN(NIN), YIN(NIN), YINNEW(NIN), YINPP(NIN)
    REAL(RKIND), intent(in), optional ::  tension
    REAL(RKIND), intent(in), optional ::  ybc(6), sigma(:)
    INTEGER, intent(in), optional :: nbc(2)
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  sigma2(NIN)
    integer :: info2
    !
    INTEGER NBCLFT,NBCRGT
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    ! 1. Deals with various optional arguments
    !
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that taus=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
       if (tension .ge. 0._rkind) then
          sigma2=tension*sigma2
       else
          zz=minval(xin(2:nin)-xin(1:nin-1))
          ! tension=-1 or negative, uses default value |tension| * Delta_x**3
          ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
          sigma2=abs(tension)*zz**3*sigma2
       end if
    else
       zz=minval(xin(2:nin)-xin(1:nin-1))
       sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbc given, ybc should also be given)
    if (present(nbc)) then
      if (present(ybc)) then
        NBCLFT = nbc(1)
        NBCRGT = nbc(2)
        YBCLFT = ybc(1)
        YBCRGT = ybc(2)
        if (NBCLFT .ge. 10) XBCLFT = ybc(3)
        if (NBCRGT .ge. 10) XBCRGT = ybc(4)
        PXEXP0 = ybc(5)
        PXEXPDL= ybc(6)
      else
        print *,' problems in interpos, nbc is given but not ybc'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = xin(1)
      XBCRGT = xin(nin)
      PXEXP0 = xin(1)-xin(nin)
      PXEXPDL= 1._rkind
    end if
    !
    ! 2. call relevant interpolation routine
    !
    CALL CBFITBND(XIN,YIN,YINNEW,NIN,YINPP,SIGMA2,NBCLFT,NBCRGT, &
      &  XBCLFT,XBCRGT,YBCLFT,YBCRGT,PXEXP0,PXEXPDL)
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_yinpp

  SUBROUTINE interpos_interp(XIN,YIN,YPP,NIN,NOUT,xout,yout,youtp,youtpp,youtint,option,info)
    ! Note: calls arg ypp instead of yinpp otherwise there is ambiguity with interpos_yinpp
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: NIN
    INTEGER, intent(in),optional :: NOUT
    REAL(RKIND), intent(in) ::  XIN(NIN), YIN(NIN), YPP(NIN)
    REAL(RKIND), intent(in), optional ::  xout(:)
    REAL(RKIND), intent(out), optional ::  yout(:), youtp(:), youtpp(:), youtint(:)
    INTEGER, intent(in), optional :: option
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND), allocatable ::  xout2(:), yout2(:), youtp2(:), youtpp2(:), youtint2(:)
    integer :: info2, nout2, optabs
    !
    INTEGER :: ioptder, iextrapo, inttype
    !
    ! 1. Deals with various optional arguments
    !
    ! xout
    if (present(xout)) then
      if (present(nout)) then
        nout2=nout
      else
        print *,'% nout not present, use size(xout) instead'
        nout2 = size(xout)
      end if
      allocate(xout2(nout2))
      xout2 = xout(1:nout2)
    else
      nout2 = nin
!      if (present(nout)) nout = nin
      allocate(xout2(nout2))
      xout2 = xin(1:nout2)
    end if
    allocate(yout2(nout2))
    allocate(youtp2(nout2))
    allocate(youtpp2(nout2))
    allocate(youtint2(nout2))
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(youtp)) ioptder=1
    if (present(youtpp)) ioptder=2
    if (present(youtint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 3
        return
      case(2,12,22,32,42)
        ! quadratic interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 4
        return
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 3
          return
        case (2)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 4
          return
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic interpolation
    CALL SPLIBNDA(XIN,YIN,YPP,NIN,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2,NOUT2, &
      & iextrapo, ioptder)
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yout)) yout(1:nout2) = yout2
    if (present(youtp)) youtp(1:nout2) = youtp2
    if (present(youtpp)) youtpp(1:nout2) = youtpp2
    if (present(youtint)) youtint(1:nout2) = youtint2
    !
    if (present(info)) info = 0
    !
    deallocate(xout2)
    deallocate(yout2)
    deallocate(youtp2)
    deallocate(youtpp2)
    deallocate(youtint2)
    return
  END SUBROUTINE interpos_interp

  SUBROUTINE interpos_defxscal(X,Y,N,xscal,tension,yscal,yscalp,yscalpp,yscalint,nbcscal,ybcscal,sigma,option,info)
    !
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: N
    REAL(RKIND), intent(in) ::  X(N), Y(N)
    REAL(RKIND), intent(in), optional ::  tension
    REAL(RKIND), intent(in) ::  xscal
    REAL(RKIND), intent(out), optional ::  yscal, yscalp, yscalpp, yscalint
    REAL(RKIND), intent(in), optional ::  ybcscal(6), sigma(:)
    INTEGER, intent(in), optional :: nbcscal(2), option
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2, yout2, youtp2, youtpp2, youtint2
    REAL(RKIND) ::  sigma2(N)
    integer :: info2, nout2
    !
    INTEGER ioptder, NBCLFT,NBCRGT, iextrapo, inttype, optabs
    REAL(RKIND) :: XBCLFT,XBCRGT, YBCLFT,YBCRGT
    REAL(RKIND) ::  ynew(N), ypp(N), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! xscal
    nout2 = 1
    xout2 = xscal
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that tension=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      if (tension .ge. 0._rkind) then
        sigma2=tension*sigma2
      else
        ! tension=-1 or negative, uses default value |tension| * Delta_x**3
        zz=minval(x(2:n)-x(1:n-1))
        ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
        sigma2=abs(tension)*zz**3*sigma2
      end if
    else
      zz=minval(x(2:n)-x(1:n-1))
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions (if nbcscal given, ybcscal should also be given)
    if (present(nbcscal)) then
      if (present(ybcscal)) then
        NBCLFT = nbcscal(1)
        NBCRGT = nbcscal(2)
        YBCLFT = ybcscal(1)
        YBCRGT = ybcscal(2)
        if (NBCLFT .ge. 10) XBCLFT = ybcscal(3)
        if (NBCRGT .ge. 10) XBCRGT = ybcscal(4)
        PXEXP0 = ybcscal(5)
        PXEXPDL= ybcscal(6)
      else
        print *,' problems in interpos, nbcscal is given but not ybcscal'
        if (present(info)) info = 2
        return
      end if
    else
      ! default is 2nd derivative=0 at each end
      ! For quadratic interpolation, default is NBCLFT=0 is used for quadratic defined with points and discontinuous 1st derivative
      NBCLFT = 0
      NBCRGT = 0
      YBCLFT = 0._rkind
      YBCRGT = 0._rkind
      ! not used but sets value to make sure:
      XBCLFT = x(1)
      XBCRGT = x(n)
      PXEXP0 = x(1)-x(n)
      PXEXPDL= 1._rkind
    end if
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      inttype=3 ! cubic
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        inttype = 1
        ! extrapolation
        if (optabs .eq. 1)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 11) iextrapo = sign(1,option) ! linear extrapolation with edge-deriv for 11; edge-val for -11
        if (optabs .eq. 21) iextrapo = sign(10,option) ! constant outside: y=yedge for option=21, y=0. for -21
      case(2,12,22,32,42)
        ! quadratic interpolation
        inttype = 2
        ! extrapolation
        if (optabs .eq. 2)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 12) iextrapo = sign(21,option) ! quadratic on dx then linear extrapol. with edge-deriv for 12; edge-val for -12
        if (optabs .eq. 22) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 32) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 32; edge-val for -32
        if (optabs .eq. 42) iextrapo = sign(10,option) ! constant outside: y=yedge for option=42, y=0. for -42
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          iextrapo = 1
        case (2)
          iextrapo = 21
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    select case (inttype)
    case (1)
      ! linear
      call intlinear(X,Y,N,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo)
    case (2)
      ! quadratic
      call intquadratic(X,Y,N,xout2,yout2,youtp2,youtpp2,youtint2,NOUT2,ioptder,iextrapo,NBCLFT)
    case (3)
      ! cubic
      ! general cubic spline routine with non-periodic boundary conditions
      CALL CBSPLGEN(X,Y,YNEW,YPP,N,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
        &    nout2,ioptder,sigma2,NBCLFT, &
        &    NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,iextrapo,PXEXP0,PXEXPDL,info2)
      !
    case default
      print *,' inttype=',inttype,' ; error in this routine, this should not be possible'
      if (present(info)) info = 4
      return
    end select
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yscal)) yscal = yout2
    if (present(yscalp)) yscalp = youtp2
    if (present(yscalpp)) yscalpp = youtpp2
    if (present(yscalint)) yscalint = youtint2
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_defxscal

  SUBROUTINE interpos_defperxscal(X,Y,N,tension,xscal,yscal,yscalp,yscalpp,yscalint,NBCSCAL,YBCSCAL,sigma,info)
    ! for periodic boundary conditions, nbcscal=-1 and ybcscal=period should be given
    ! option not needed since no extrapolation and only valid for cubic spline here
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: N, NBCSCAL
    REAL(RKIND), intent(in) ::  X(N), Y(N)
    REAL(RKIND), intent(in), optional ::  tension
    REAL(RKIND), intent(in) ::  xscal
    REAL(RKIND), intent(out), optional ::  yscal, yscalp, yscalpp, yscalint
    REAL(RKIND), intent(in), optional ::   ybcscal, sigma(:)
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2, yout2, youtp2, youtpp2, youtint2
    REAL(RKIND) ::  sigma2(N)
    integer :: info2, nout2, n_eff
    !
    INTEGER ioptder
    REAL(RKIND) :: period
    REAL(RKIND) ::  ynew(N), ypp(N), sigmamin
    REAL(RKIND) :: PXEXP0,PXEXPDL, zz
    !
    info2 = -1
    !
    ! 1. Deals with various optional arguments
    !
    ! period
    n_eff = n
    if (present(ybcscal)) then
      period = ybcscal
      ! check if extra point is given, if yes skip it with warning:
      if (abs(x(n_eff)-x(1)-period) .lt. 1.0e-10_rkind*(x(2)-x(1))) then
        print *,'% it seems periodic point x(n)=x(1)+period is included in input, use n-1 points'
        n_eff = n_eff-1
      end if
    else
      ! Assumes last x point corresponds to x(1)+period and is redundant
      period=x(n_eff)-x(1)
      n_eff = n_eff-1
    end if
    ! xscal
    nout2 = 1
    xout2 = xscal
    ! sigma
    if (present(sigma)) then
      ! normalize to minimum value, so that tension=tension at this position
      sigmamin = minval(sigma)
      if (sigmamin .le. 0._rkind) then
        print *,' ERROR: min(sigma)=',sigmamin,' is .le. 0.'
        if (present(info)) info = 1
        return
      end if
      sigma2 = sigma/sigmamin
    else
      sigma2 = 1._rkind
    end if
    ! tension
    if (present(tension)) then
      if (tension .ge. 0._rkind) then
        sigma2=tension*sigma2
      else
        ! tension=-1 or negative, uses default value |tension| * Delta_x**3
        zz=minval(x(2:n_eff)-x(1:n_eff-1))
        ! tension=abs(tension)*zz**3 ! if intent(inout) but does not allow call: tension=-1.
        sigma2=abs(tension)*zz**3*sigma2
      end if
    else
      zz=minval(x(2:n_eff)-x(1:n_eff-1))
      sigma2=zz**3*sigma2
    end if
    !
    ! boundary conditions
    if (nbcscal .ne. -1) then
      print *,'nbcscal when scalar should be -1 for periodic boundary conditions, otherwise vector of size 2'
      if (present(info)) info = 2
      return
    end if

    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic
    ! general cubic spline routine with periodic boundary conditions
    CALL CBSPLGNP(X,Y,YNEW,YPP,N_eff,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2, &
      &    nout2,ioptder,sigma2,period,PXEXP0,PXEXPDL,info2)
      !
    !
    ! 3. Fill in output values
    !
    ! yscal, p, pp, int
    if (present(yscal)) yscal = yout2
    if (present(yscalp)) yscalp = youtp2
    if (present(yscalpp)) yscalpp = youtpp2
    if (present(yscalint)) yscalint = youtint2
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_defperxscal

  SUBROUTINE interpos_interpxscal(X,Y,YPP,N,xscal,yscal,yscalp,yscalpp,yscalint,option,info)
    ! Note: calls arg ypp instead of yinpp otherwise there is ambiguity with interpos_yinpp
    USE prec_rkind
    implicit none
    INTEGER, intent(in) :: N
    REAL(RKIND), intent(in) ::  X(N), Y(N), YPP(N)
    REAL(RKIND), intent(in) ::  XSCAL
    REAL(RKIND), intent(out), optional ::  yscal, yscalp, yscalpp, yscalint
    INTEGER, intent(in), optional :: option
    INTEGER, intent(out), optional :: info
    ! need variables computed even if not given as argument
    REAL(RKIND) ::  xout2, yout2, youtp2, youtpp2, youtint2
    integer :: info2, nout2, optabs
    !
    INTEGER :: ioptder, iextrapo, inttype
    !
    ! 1. Deals with various optional arguments
    !
    ! xscal
    nout2 = 1
    xout2 = xscal
    ! compute y (0), also yp (1), also ypp (2) or also yint(3)
    ioptder = 0
    if (present(yscalp)) ioptder=1
    if (present(yscalpp)) ioptder=2
    if (present(yscalint)) ioptder=3
    ! option: interpolation and extrapolation choice
    if (.NOT. present(option)) then
      iextrapo = 32 !extrapolation cubic on 1st dx outside and then quadratic
    else
      optabs = abs(option)
      select case (optabs)
      case(1, 11, 21)
        ! linear interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 3
        return
      case(2,12,22,32,42)
        ! quadratic interpolation
        print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
        if (present(info)) info = 4
        return
      case(3,13,23,33,43,53,63)
        ! cubic interpolation
        inttype = 3
        ! extrapolation
        if (optabs .eq. 3)  iextrapo = 0 ! no extrapolation
        if (optabs .eq. 13) iextrapo = sign(32,option) ! cubic on dx then quadratic extrapol. with edge-deriv for 13; edge-val for -13
        if (optabs .eq. 23) iextrapo = sign(1,option) ! linear extrapol. with edge-deriv for 22; edge-val for -22
        if (optabs .eq. 33) iextrapo = sign(2,option) ! quadratic extrapol. with edge-deriv for 33; edge-val for -33
        if (optabs .eq. 43) iextrapo = sign(3,option) ! cubic extrapol. with edge-deriv for 43; edge-val for -43
        if (optabs .eq. 53) iextrapo = sign(31,option) ! cubic on dx then linear extrapol. with edge-deriv for 53; edge-val for -53
        if (optabs .eq. 63) iextrapo = sign(10,option) ! constant outside: y=yedge for option=63, y=0. for -63
      case default
        inttype=mod(optabs,10)
        select case (inttype)
        case (1)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 3
          return
        case (2)
          print *,' should call it only for cubic since uses ypp already calculated by cubic interpolation assumption'
          if (present(info)) info = 4
          return
        case (3)
          iextrapo = 32
        case default
          print *,'option= ',option, &
            & ' is not valid. The last digit should be 1, 2 or 3 for linear, quadratic or cubic interpolation'
          if (present(info)) info = 3
          return
        end select
      end select
    end if
    !
    ! 2. call relevant interpolation routine
    !
    ! cubic interpolation
    CALL SPLIBNDA(X,Y,YPP,N,XOUT2,YOUT2,YOUTP2,YOUTPP2,YOUTINT2,NOUT2, &
      & iextrapo, ioptder)
    !
    ! 3. Fill in output values
    !
    ! yout, p, pp, int
    if (present(yscal)) yscal = yout2
    if (present(yscalp)) yscalp = youtp2
    if (present(yscalpp)) yscalpp = youtpp2
    if (present(yscalint)) yscalint = youtint2
    !
    if (present(info)) info = 0
    !
    return
  END SUBROUTINE interpos_interpxscal

end MODULE interpos_module
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
  USE prec_rkind
  USE interpos_module
  implicit none
  !   the mex function should have the name of the file. In this case interpos
  !   Then typically in matlab one calls:
  !   >> [a,b,c,..]=interpos(arg1,arg2,...)
  !
  !   nlhs, nrhs: number of input arguments on left and right-hand side respectively
  !   thus one can test and depending on number of arguments choose different options
  !
  !   plhs, prhs: pointer to the different arguments: plhs(i), i=1,nlhs
  !   (plhs and prhs are integer arrays)
  !
  !   fmex function: mxGetPr gets the pointer value in prhs, plhs
  !   .              mxCopyPtrToReal8: copy pointed values to local array
  !   .              mxCopyReal8ToPtr: copy local array values to pointed array
  !   .              mxGetM: get nb of rows in pointed array
  !   .              mxGetN: get nb of columns in pointed array
  !   .              mxcreatedoublematrix: creates a matrix for the return argument
  !
  !   In this example, matlab call expected:
  !
  !   >> [yout{,youtp,youtpp,youtint}] = interpos({kopt,}xin,yin{,xout}{,taus{,nbc,ybc}{,sigma}})
  !
  !   where {..} means facultative, thus one can have for example the following calls:
  !
  !   >> [yout]=interpos(xin,yin,yout)
  !   >> [yout]=interpos(kopt,xin,yin,xout)
  !   >> [yout]=interpos(xin,yin,taus)
  !   >> [yout,youtp]=interpos(xin,yin)
  !   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,-1.)
  !   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,taus,nbc,ybc)
  !   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,xout,taus,nbc,ybc,sigma)
  !
  !   NOTE: If xout is a single point, taus should also be given otherwise it assumes it is a taus value:
  !          interpos(xin,yin,xxx) means taus=xxx if xxx size of 1
  !          interpos(xin,yin,xxx,0.) means xout=xxx, taus=0. if xxx size of 1
  !
  !   The arguments have the following meaning with respect to the function y(x)
  !   that one wants to interpolate/extrapolate or compute its derivatives and integrals:
  !
  !   kopt: option for interpolation and extrapolation method (default is kopt=13):
  !   .     = 1, 11 or 21 : linear interpolation
  !   .     = 2, 12 or 22 : quadratic "
  !   .     = 3, 13 or 23 : cubic spline interpolation, with tension=taus
  !   .     < 10 : send warning if need to extrapolate
  !   .     in [11,19] : extrapolate using a lower order than for interpolation
  !   .     in [21,29] : extrapolate using same order as interpolation
  !
  !   xin : array giving the input values of x
  !   yin : array giving the input values of y(x=xin(i))
  !
  !   xout: array of values of x at which the function and/or its derivatives
  !   .     are to be computed. If xout is not given, assumes xout=xin
  !   yout: interpolated values of y at x=xout(i).
  !   youtp: 1st derivative of y at x=xout(i)
  !   youtpp: 2nd derivative of y at x=xout(i)
  !   youtint: Integral of (y dx) from x=xin(1) to x=xout(i)
  !
  !   taus  : tension value for cubic spline interpolation. If not given, uses taus=0
  !           if taus < 0 (typically -1) uses a default taus value: taus=abs(taus)*default_taus
  !           (the default_taus is typically min(delta_x)**3)
  !
  !   nbc(2): [NBCLFT NBCRGT] (default: [0 0])
  !     BOUNDARY CONDITIONS, 4 TYPES DETERMINED BY THE VALUE OF (nbc(1)=NBCLFT
  !     and nbc(2) for left and right-hand side BC.):
  !
  !     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
  !     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
  !     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
  !     if nbc(1)=-1, then assumes periodic boundary condition and uses ybc(1) for the period
  !
  !     The value of nbc(1 or 2) should be >= 10 if the BC is not at an end point:
  !     XBCLFT~=xin(1) or XBCRGT~=xin(end)
  !
  !     Examples:
  !       A good value for radial profiles in rho between [0,1] is to specify
  !       the first derivative = 0 at left and second derivative=0 at right:
  !       => nbc = [1 0] and ybc=[0. 0.]
  !
  !       or to obtained the first derivative at right-hand side from a
  !       lagrangian interpolation of the last points:
  !       => nbc=[1 1] and ybc=[0. 1e32]
  !
  !       Note: The B.C. can only be given anywhere but within the interval [xin(1),xin(end)]
  !       Therefore, say rho=[0.1 ... 1] and one wants to impose zero
  !       derivative at rho=0, one should add a point in the input before. One uses sigma to avoid forcing the value
  !       calling interpos:
  !             rho_eff(1)=0.;
  !             rho_eff(2:length(rho)+1) = rho;
  !       then  [yout]=interpos([0.;rho],[yin(1);yin],..,taus,[1 0],[0. 0.],[1000;ones(size(rho))]);
  !
  !   ybc(2 or 4 or 6) (default: [0 0]): [YBCLFT YBCRGT] or [YBCLFT YBCRGT XBCLFT XBCRGT] with
  !     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
  !
  !     FOR nbc type 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
  !     FOR nbc type 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
  !
  !     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
  !     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
  !
  !   ybc(6): [YBCLFT YBCRGT XBCLFT XBCRGT PXEXP0 PXEXPDL]
  !        Enables one to specify a gaussian wight to taus with PXEXP0 and PXEXPDL:
  !              taus_eff =  taus * EXP(-((XIN-PXEXP0)/PXEXPDL)**2)
  !        IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND taus_eff=cst= taus
  !        if ybc(5:6) not given, then  PXEXP0=xin(1)-1. and pxexpdl=1. to get constant taus
  !
  !   ybc(1)=period for periodic boundary condition, if nbc(1)=-1. Defines the condition that y(x+period)=y(x)
  !          In this case xin(end) should not be equal to xin(1)+period, since it is redundant. However interpos just
  !          does not use the end point in this case.
  !        if period=-1, assumes period=xin(end)-xin(1). Useful if xin goes from 0 to 2pi for example
  !
  !    sigma: error_bar at each (x,y) point used for the fit. The effective taus value will then
  !           be taus(i)=taus .* sigma(i) ./ min(sigma). So uses the relative values of sigma.
  !
  !
  CRONOSINT plhs(*), prhs(*)
  
      MATLABINOUTDIM nlhs,nrhs
  !
  CRONOSINT mxGetPr
  !   pointers for matlab input arguments
  CRONOSINT pkopt, pxin, pyin, pxout, ptaus, pnbc, pybc, psig
  REAL(RKIND),pointer :: xin(:), yin(:), xout(:), sig(:)
  REAL(RKIND) :: ybc(6)
  integer kopt, nbc(2)
  REAL(RKIND) :: taus, zkopt, znbc(2)
  !   pointers for matlab output arguments
  CRONOSINT pyout, pyoutp, pyoutpp, pyoutint
  REAL(RKIND),pointer :: yout(:), youtp(:), youtpp(:), youtint(:)
  REAL(RKIND) :: XX, ZTAUEFF, ZCOFEXP, ZXEXP0, ZXEXPDL
  CRONOSINT mxgetm, mxgetn, mxcreatedoublematrix
  CRONOSINT malloc_f
  integer nbyt_tot, infomalloc
  common /memoryuse/ nbyt_tot
  REAL(RKIND) :: sigmin, fun_sigma
  integer kopt_sign, inttype, iextrapo, option, ninrow, nincol, nin, &
    &  nout, ixoutxin, i1len, i4len, i, inextrhs, idoexp_error, iybclen, &
    &  noutrow, noutcol, ioptder, iflag, i5len, i4isxout, i4arg
  integer nel
  CRONOSINT ninint8, noutint8, nelint8, iybclenint8, noutrowint8, &
    &  noutcolint8, ninrowint8, nincolint8
!!$      integer :: inan
!!$      real :: znan
!!$      data inan/B'01111111100000100000000000000000'/ ! ok on lac g95
!!$!      data inan/B'01111111100100010001001010101010'/  ! ok also (0 as first)
!!$!       data inan/B'01111111100000000000000000000000'/ ! +Inf
  !
  fun_sigma(XX)= ZTAUEFF*EXP(-ZCOFEXP*(XX-ZXEXP0)**2/ZXEXPDL**2)
  !
  !.......................................................................
  !
  !   1. defaults values
  !
  !%OS      
  !write(0,*) ' nrhs= ',nrhs
  !%OS      
  !write(0,*) 'prhs= ',prhs(1),prhs(2),prhs(3)
  nbyt_tot=0
  taus = 0._RKIND
  !%OS      print *,' nbyt_tot= ',nbyt_tot
  !
  !   2. Input arguments
  !
  !%OS      print *,' nb of return arguments= ',nlhs
  !%OS      print *,' nb of input arguments= ',nrhs
  !
  !   2.1 kopt
  !
  inextrhs=1
  i1len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
  if (i1len .eq. 1) then
    ! 1st input is kopt
    pkopt=mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    nelint8 = int(1,8)
    call mxCopyPtrToReal8(pkopt, zkopt, nelint8)
    kopt = int(abs(zkopt))
    kopt_sign=sign(1,int(zkopt))
  else
    ! 1st input is xin, thus define default for kopt
    kopt = 33 ! 3 for cubic, large tenth for default extrapolation
    kopt_sign = 1
    inextrhs = 1
  end if
  !
  !%OS      
  !write(0,*) ' kopt= ',kopt,zkopt
  !
  !   define interpolation type and extrapolation type
  !
  nel=10
  inttype = mod(kopt,10)
  iextrapo = kopt/nel
  !  print *,' kopt, iextrapo, inttype = ',kopt,iextrapo,inttype 
  !
  !   2.2 get xin, first length and allocate space to yin as well
  !
  pxin=mxGetPr(prhs(inextrhs))
  ninrowint8 = mxGetM(prhs(inextrhs))
  nincolint8 = mxGetN(prhs(inextrhs))
  inextrhs = inextrhs + 1
  ninrow = int(ninrowint8)
  nincol = int(nincolint8)
  !%OS      print *,'nincol, ninrow= ',nincol, ninrow
  !   input 1D arrays in row or column
  nin = max(ninrow,nincol)
  allocate(xin(nin))
  allocate(yin(nin))
  ninint8 = int(nin,8)
  call mxCopyPtrToReal8(pxin, xin, ninint8)
  !%OS      print *,xin(1:10)
  !%OS      
  !write(0,*) ' nin= ',nin
  !%OS      print *,' xin= ',(xin(i),i=1,nin)
  !
  !   2.3 get yin
  !
  pyin=mxGetPr(prhs(inextrhs))     ! Get yin
  inextrhs = inextrhs + 1
  call mxCopyPtrToReal8(pyin, yin, ninint8)
  !%OS      print *,' yin= ',(yin(i),i=1,nin)
  !
  !   2.4 check for 4th input
  !
  !%OS      print *,' nrhs= ',nrhs
  if (nrhs .lt. inextrhs) then
    nout = nin
    !   flag for xout=xin set to 1
    ixoutxin = 1
    !%OS        print *,' ixoutxin= ',ixoutxin
  else
    !
    ! check if next argument is an array, then it is xout.
    ! If it is a scalar and next one as well, then xout scalar and taus given, otherwise xout=xin and taus given
    !
    i4len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
    i4arg = inextrhs
    inextrhs = inextrhs + 1
    !%OS      print *,'i4len= ',i4len
    i4isxout = 0
    if (i4len.gt.1) i4isxout = 1
    if (nrhs.ge.inextrhs) then
      i5len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
      if (i4len.eq.1 .and. i5len.eq.1) then
        i4isxout = 1
      endif
    endif
    ! print *,'i4len, i5len, i4isxout= ',i4len, i5len, i4isxout
    if (i4isxout .eq. 0) then
      !   in this case prhs(4) gives the tau value  
      ixoutxin = 1
      ptaus = mxGetPr(prhs(i4arg))
      nelint8=int(1,8)
      call mxCopyPtrToReal8(ptaus, taus, nelint8)
      !%OS       print *,'salut11'
    else
      ixoutxin = 0
      nout = i4len
      allocate(xout(nout))
      pxout = mxGetPr(prhs(i4arg))
      noutint8 = int(nout,8)
      call mxCopyPtrToReal8(pxout, xout, noutint8)
    endif
  endif
  if (ixoutxin .eq. 1) then
    nout = nin
    noutint8 = int(nout,8)
    allocate(xout(nout))
    do i=1,nout
      xout(i) = xin(i)
    end do
  endif
  ! print *,'ixoutxin, nout= ',ixoutxin, nout
  !%OS      print *,' xout= ',(xout(i),i=1,nout)
  !
  !   2.5 5th argument is taus if xout .ne. xin
  !
  if (nrhs.ge.inextrhs .and. ixoutxin.eq.0) then
    ptaus = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    nelint8=int(1,8)
    call mxCopyPtrToReal8(ptaus, taus, nelint8)
  endif
  !%OS      
  !write(0,*) 'taus= ',taus
  !
  !   2.6 5,6th or 6,7th arguments are nbc, ybc
  !
  nbc(1) = 0
  nbc(2) = 0
  ybc(1) = 0._RKIND
  ybc(2) = 0._RKIND
  idoexp_error=0
  if (nrhs .ge. inextrhs) then
    i4len = int(max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs))))
!!$         print *,'i4len= ',i4len
    pnbc = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    if (i4len .ne. 0) then
      if (i4len .eq. 1) then
        nelint8=int(1,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(1)
        ybc(1) = -1.
      else
        nelint8=int(2,8)
        call mxCopyPtrToReal8(pnbc, znbc, nelint8)
        nbc(1) = znbc(1)
        nbc(2) = znbc(2)
      end if
    end if
!!$ ! to fill in with NaNs when iextrapol=0 and xout outside interval, could use: (but may not work on all compiler)
!!$         taus=-1._rkind
!!$         taus=real(sqrt(taus),rkind)
    !
!!! can use this instead with data inan/B'01111111100100010001001010101010'/ or data inan/B'01111111100000100000000000000000'/
!!$         znan=transfer(inan,znan)
!!$         print *,'znan= ',znan,' inan= ',inan

!!$         print *,'taus,nbc= ',taus,nbc,' znbc= ',znbc
    !
    if (nrhs .ge. inextrhs) then
      iybclenint8 = max(mxGetM(prhs(inextrhs)),mxGetN(prhs(inextrhs)))
!!$            print *,' iybclenint8= ',iybclenint8
      pybc = mxGetPr(prhs(inextrhs))
      inextrhs = inextrhs + 1
      if (iybclenint8 .eq. 0) then
        ybc(1) = 0._rkind
        ybc(2) = 0._rkind
        idoexp_error = 0
      else
        call mxCopyPtrToReal8(pybc, ybc, iybclenint8)
        iybclen=int(iybclenint8)
        if (iybclen .le. 4) then
          ybc(5) = xin(1) - xin(nin)
          ybc(6) = 1._RKIND
        else
          idoexp_error = 1
        endif
      end if
    end if
  else
    ybc(5) = xin(1) - xin(nin)
    ybc(6) = 1._RKIND
  endif
!!$      
   !write(0,*) 'ybc=',ybc
  !
  !   2.7 7th or 8th is sigmain
  allocate(sig(nin))
  if (nrhs .ge. inextrhs) then
    psig = mxGetPr(prhs(inextrhs))
    inextrhs = inextrhs + 1
    call mxCopyPtrToReal8(psig, sig, ninint8)
!!$    sigmin=sig(1)
!!$    do i=2,nin
!!$      sigmin=min(sigmin,sig(i))
!!$    end do
!!$    do i=1,nin
!!$      sig(i) = sig(i)/sigmin*abs(taus)
!!$    end do
  else
    ! now multiplication by taus done in fortran interpos
    ZTAUEFF = 1._rkind
    ZXEXP0 = ybc(5)
    ZXEXPDL = ybc(6)
    ZCOFEXP = 0._RKIND
    if (idoexp_error.EQ.1) THEN
      ZCOFEXP = 1.0_RKIND
    endif
    do i=1,nin
      sig(i)=fun_sigma(xin(i))
    end do
  endif
  ! print *,'sig= ',sig
!!$      if (taus.LT.0._RKIND) then
!!$        sig(1)=10._RKIND*sig(1)
!!$      endif
  !
  !   3. Allocate return arguments arrays (in same geometry as xin)
  !
  if (ninrow .eq. 1) then
    noutrow = 1
    noutcol = nout
  else
    noutrow = nout
    noutcol = 1
  endif
  !   pointers for yout, youtp and youtpp and youtint respectively
  nelint8=int(0,8)
  noutrowint8 = int(noutrow,8)
  noutcolint8 = int(noutcol,8)
  !%OS      print *,' noutrowint8, noutcolint8= ',noutrowint8, noutcolint8
  !%OS      print *,' noutrow, noutcol= ',noutrow, noutcol
  plhs(1) = mxcreatedoublematrix(noutrowint8,noutcolint8,nelint8)
  !%OS      print *,plhs(1)
  pyout = mxGetPr(plhs(1))
  allocate(yout(nout))
  yout=0._RKIND
  ioptder = 0
  if (nlhs .ge. 2) then
    nelint8=int(0,8)
    plhs(2) = mxcreatedoublematrix(noutrowint8,noutcolint8,nelint8)
    pyoutp = mxGetPr(plhs(2))
    allocate(youtp(nout))
    youtp=0._RKIND
    ioptder = 1
  endif
  if (nlhs .ge. 3) then
    nelint8=int(0,8)
    plhs(3) = mxcreatedoublematrix(noutrowint8,noutcolint8,nelint8)
    pyoutpp = mxGetPr(plhs(3))
    allocate(youtpp(nout))
    youtpp=0._RKIND
    ioptder = 2
  endif
  if (nlhs .ge. 4) then
    nelint8=int(0,8)
    plhs(4) = mxcreatedoublematrix(noutrowint8,noutcolint8,nelint8)
    pyoutint = mxGetPr(plhs(4))
    allocate(youtint(nout))
    youtint=0._RKIND
    ioptder = 3
  endif
  !
  !   4. compute interpolated function and its derivatives
  !
  !  print *,' inttype, ioptder,iextrapo= ',inttype, ioptder,iextrapo
  if (inttype .eq. 1) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 1
    else if (iextrapo .eq. 2) then
      iextrapo = 10
    else
      iextrapo = 1
    end if
    iextrapo = kopt_sign * iextrapo
    call intlinear(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout,ioptder,iextrapo)
  end if
  !
  if (inttype .eq. 2) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 1
    else if (iextrapo .eq. 2) then
      iextrapo = 2
    else if (iextrapo .eq. 3) then
      iextrapo = 21
    else if (iextrapo .eq. 4) then
      iextrapo = 10
    else
      iextrapo = 2
    end if
    iextrapo = kopt_sign * iextrapo
    call intquadratic(xin,yin,nin,xout,yout,youtp,youtpp,youtint,nout,ioptder,iextrapo,nbc(1))
  end if
  !
  if (inttype .eq. 3) then
    if (iextrapo .eq. 0) then
      iextrapo = 0
    else if (iextrapo .eq. 1) then
      iextrapo = 3
    else if (iextrapo .eq. 2) then
      iextrapo = 31
    else if (iextrapo .eq. 3) then
      iextrapo = 32
    else if (iextrapo .eq. 4) then
      iextrapo = 1
    else if (iextrapo .eq. 5) then
      iextrapo = 2
    else if (iextrapo .eq. 6) then
      iextrapo = 10
    else
      iextrapo = 32
    endif
    !        print *,'iextrapo= ',iextrapo
    !
    !   add sign of kopt for extrapol
    !%OS        print *,' ioptder,iextrapo,inttype= ',ioptder,iextrapo,inttype
    !%OS        print *,'xin,yin,nin,xout',xin,yin,nin,xout
    if (nbc(1).lt. 0) then
      ! periodic boundary conditions
      if (ioptder .eq. 0) then
        call interpos(xin,yin,nin,nout,taus,xout,yout &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else if (ioptder .eq. 1) then
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else if (ioptder .eq. 2) then
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      else
        call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp,youtint &
          &    ,nbc=nbc(1),ybc=ybc(1),sigma=sig,info=iflag)
      end if
    else
      if (ioptder .eq. 0) then
        call interpos(xin,yin,nin,nout,taus,xout,yout &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else if (ioptder .eq. 1) then
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else if (ioptder .eq. 2) then
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      else
      call interpos(xin,yin,nin,nout,taus,xout,yout,youtp,youtpp,youtint &
          &    ,nbc=nbc,ybc=ybc,sigma=sig,option=kopt*kopt_sign,info=iflag)
      end if
    endif
    if (iflag .ne. 0) then
      print *
      print *,' problem in cbsplgen'
      print *,' iflag = ',iflag
      return
    end if
  end if
  !   
  !   5. copy arrays to return value pointers
  !

  !%OS      print *,'yout(1:10)',yout(1:10)
  !%OS      print *,'pyout= ',pyout
  !%OS      print *,'noutint8= ',noutint8
  call mxCopyReal8ToPtr(yout, pyout, noutint8)
  if (ioptder .ge. 1) then
    call mxCopyReal8ToPtr(youtp,pyoutp,noutint8)
  endif
  if (ioptder .ge. 2) then
    call mxCopyReal8ToPtr(youtpp,pyoutpp,noutint8)
  endif
  if (ioptder .ge. 3) then
    call mxCopyReal8ToPtr(youtint,pyoutint,noutint8)
  endif
  !
  deallocate(xin)
  deallocate(yin)
  deallocate(sig)
  !
  return
end subroutine mexfunction
SUBROUTINE CBSPLGEN(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
  &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PSIG,NBCLFT &
  &  ,NBCRGT,XBCLFT,XBCRGT,YBCLFT,YBCRGT,KEXTRPOL,PXEXP0,PXEXPDL &
  &  ,KFLAG)
  !     =================================================================
  !
  !     NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PTAUS=0 (i.e. psig=0):
  !           THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
  !           (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
  !
  !     Interpolate (pxin,pyin) on (pxout,pyout) using
  !     Hirshman fitted cubic spline with ptaus value or
  !     standard cubic spline if PTAUS=0 (psig=ptaus*sigma/sigmamin)
  !
  !     KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
  !     KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
  !     KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
  !     KOPT = 3: AS KOPT=2 PLUS INTEGRAL FROM PXIN(1) UP TO PXOUT(J)
  !
  !     SEE COMMENTS FOR ROUTINE CBFITBND FOR MORE INFORMATION
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     PXIN    : INPUT ABSCISSA (GIVEN DATA)
  !     PYIN    : INPUT VALUES OF FUNCTION AT PXIN(I),I=1,KNIN
  !     PYINNEW : IF PTAUS.NE.0, THEN PYNEW CONTAINS ON OUTPUT THE NEW VALUES
  !     .         OF THE FUNCTION AT PXIN(I) FOR THE CUBIC SPLINE FIT
  !     PYINPP  : SECOND DER. OF THE CUBIC SPLINE FIT FOR (PXIN,PYIN) IF PTAUS=0 OR
  !     .         ON (PXIN,PYNEW) OTHERWISE
  !     KNIN    : NUMBER OF INPUT POINTS
  !     PXOUT   : X VALUES AT WHICH THE FUNCTION HAS TO BE INTERPOLATED (INPUT)
  !     PYOUT   : INTERPOLATED VALUES AT PXOUT(I),I=1,KNOUT (OUTPUT)
  !     PYOUTP  : INTERPOLATED VALUES OF 1ST DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.1)
  !     PYOUTPP : INTERPOLATED VALUES OF 2ND DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.EQ.2)
  !     PYOUTINT: INTEGRAL OF Y FROM PXIN(1) TO PXOUT(J)
  !     KNOUT   : NUMBER OF POINTS FOR OUTPUT
  !     KOPT    : SEE ABOVE
  !     PSIG    : SIGMA at each point, normalized by minimum sigma, times PTAUS
  !     PTAUS   : WEIGHT OF SECOND DERIVATIVE IN THE CHI**2 TO BE MINIMIZED. PTAUS=0 GIVES THE
  !     .         STANDARD CUBIC SPLINE. LARGER VALUES OF PTAUS WILL SMOOTH MORE THE 2ND DER.
  !     NBCLFT  : FOR LEFT B.C., VALUE SHOULD BE 0,1,2,10,11 OR 12. (SEE ROUTINE CBFITBND BELOW)
  !     NBCRGT  : FOR RIGHT B.C. (SEE ROUTINE CBFITBND BELOW)
  !     XBCLFT  : FOR LEFT B.C., USED ONLY IF NBCLFT.GE.10 (SEE ROUTINE CBFITBND BELOW)
  !     XBCRGT  : FOR RIGHT B.C., USED ONLY IF NBCRGT.GE.10 (SEE ROUTINE CBFITBND BELOW)
  !     YBCLFT  : VALUE OF LEFT B.C.
  !     YBCRGT  : VALUE OF RIGHT B.C.
  !     .         STANDARD B.C. (SECOND DER. = 0) IS OBTAINED WITH:
  !     .         NBCLFT = NBCRGT = 0 AND YBCLFT = YBCRGT = 0.
  !     KEXTRPOL: OPTION ON HOW TO EXTRAPOLATE THE FUNCTION IF PXOUT(I) IS OUTSIDE [PXIN(1),PXIN(KNIN)]
  !     .       = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
  !     .       = 1: LINEAR EXTRAPOLATION
  !     .       = 2: USE QUADRATIC INTERPOLATION IF X OUT OF BOUND
  !     .       = 3: USE CUBIC INTERPOLATION IF X OUT OF BOUND
  !     .       = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
  !     .       = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
  !     .       = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
  !     PXEXP0  : PTAUS IS WEIGHTED BY AN EXP(-((X-PXEXP0)/PXEXPDL)**2)
  !     PXEXPDL : IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND PTAUS=CST
  !     .         (SEE ROUTINE CBFITBND BELOW)
  !     .         GIVE PXEXP0=PXIN(1)-1. AND PXEXPDL=1. TO GET CONSTANT PTAUS
  !     KFLAG   : ERROR FLAG: IF NOT 0, THERE IS A PROBLEM
  !
  !-----------------------------------------------------------------------
  USE prec_rkind
  implicit none
  REAL(RKIND) :: EPSILON
  PARAMETER(EPSILON = 1.0E-10_RKIND)
  !
  integer knin, knout, kopt, mdmatot, nbclft, nbcrgt, kextrpol, &
    &  kflag, i, ioptmono
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), PYINPP(KNIN), &
    &  PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), &
    &  PYOUTINT(KNOUT), PSIG(KNIN)
  !   
  REAL(RKIND) :: xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl, zy, zyp, &
    &  zypp
  !
  !-----------------------------------------------------------------------
  !     0. CHECK INPUT CONSISTENCY
  !
  KFLAG = 0
  IF (PSIG(1) .EQ. 0._RKIND) THEN
    IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
      &  .OR. NBCRGT.EQ.2) THEN
      PRINT *,' PTAUS=0, BUT NEED SMOOTHING, WHEN'
      PRINT *,'     NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR.', &
        &      ' NBCLFT.EQ.2 .OR. NBCRGT.EQ.2'
      PRINT *,' NBCLFT = ',NBCLFT
      PRINT *,' NBCRGT = ',NBCRGT
      !%OS          STOP 'TAU=0'
      KFLAG = 1
      return
    ENDIF
  ENDIF
  !
  !   PXIN in ASCENDING ORDER
  !
  DO i=1,KNIN-1
    if (PXIN(i) .GE. PXIN(i+1)) then
      print *,' xin not in ascending order:'
      print *,' xin(',i,')= ',PXIN(i),'   >=   xin(',i+1,')= ', &
        &      PXIN(i+1)
      KFLAG = 2
      RETURN
    endif
  END DO
  !
  CALL CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG,NBCLFT,NBCRGT, &
    &  XBCLFT,XBCRGT,YBCLFT,YBCRGT,PXEXP0,PXEXPDL)
  !
  !L    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT
  !
  IF (PSIG(1) .EQ. 0.0_RKIND) THEN
    CALL SPLIBNDA(PXIN,PYIN   ,PYINPP,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT, &
      &      KEXTRPOL,KOPT)
  ELSE
    CALL SPLIBNDA(PXIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT, &
      &      KEXTRPOL,KOPT)
  ENDIF
  !
  RETURN
END SUBROUTINE CBSPLGEN
SUBROUTINE CBSPLGNP(PXIN,PYIN,PYINNEW,PYINPP,KNIN,PXOUT,PYOUT, &
  &  PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPT,PTAUS,PERIOD,PXEXP0,PXEXPDL,KFLAG)
  !     =================================================================
  !
  !   ASSUME PERIODIC BOUNDARY CONDITION WITH Y(PXIN(1)+PERIOD) = PYIN(1)
  !   BUT PYIN(KNIN).NE. PYIN(1), SO DON-T INCLUDE EXTRA PERIODIC POINT
  !
  !     NOTE: THIS ROUTINE INCLUDES THE STANDARD CUBIC SPLINE IF PTAUS=0:
  !           THEN PYINNEW IS NOT USED AND MDAMAT=3 IS SUFFICIENT
  !           (=> PYINNEW(1) OR PYINNEW=PYIN IS OK)
  !
  !     Interpolate (pxin,pyin) on (pxout,pyout) using
  !     Hirshman fitted cubic spline with ptaus value or
  !     standard cubic spline if PTAUS=0
  !
  !     KOPT = 0: ONLY INTERPOLATE FUNCTION INTO PYOUT
  !     KOPT = 1: INTERPOLATE FUNCTION INTO PYOUT AND 1ST DER. INTO PYOUTP
  !     KOPT = 2: AS KOPT=1 PLUS 2ND DER. INTO PYOUTPP
  !     KOPT = 3: AS KOPT=2 PLUS PRIMITIVE(Y SPLINE) INTO PYOUTINT
  !
  !     SEE COMMENTS FOR ROUTINE CBFITBND FOR MORE INFORMATION
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     PXIN(KNIN)    : INPUT ABSCISSA (GIVEN DATA)
  !     PYIN(KNIN)    : INPUT VALUES OF FUNCTION AT PXIN(I),I=1,KNIN
  !     PYINNEW(KNIN)   : IF PTAUS.NE.0, THEN PYINNEW CONTAINS ON OUTPUT THE NEW VALUES
  !     .         OF THE FUNCTION AT PXIN(I) FOR THE CUBIC SPLINE FIT
  !     PYINPP(KNIN+1)  : SECOND DER. OF THE CUBIC SPLINE FIT FOR (PXIN,PYIN) IF PTAUS=0 OR
  !     .         ON (PXIN,PYINNEW) OTHERWISE
  !     KNIN    : NUMBER OF INPUT POINTS
  !     PXOUT(KNOUT)   : X VALUES AT WHICH THE FUNCTION HAS TO BE INTERPOLATED (INPUT)
  !     PYOUT(KNOUT)   : INTERPOLATED VALUES AT PXOUT(I),I=1,KNOUT (OUTPUT)
  !     PYOUTP(KNOUT)  : INTERPOLATED VALUES OF 1ST DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.1)
  !     PYOUTPP(KNOUT) : INTERPOLATED VALUES OF 2ND DER. OF FUNCTIONS AT PXOUT(I) (OUTPUT, IF KOPT.GE.2)
  !     PYOUTINT(KNOUT): INTEGRAL OF YOUT AT PXOUT(I) (OUTPUT, IF KOPT.EQ.3)
  !     KNOUT   : NUMBER OF POINTS FOR OUTPUT
  !     KOPT    : SEE ABOVE
  !     PTAUS   : WEIGHT OF SECOND DERIVATIVE IN THE CHI**2 TO BE MINIMIZED. PTAUS=0 GIVES THE
  !     .         STANDARD CUBIC SPLINE. LARGER VALUES OF PTAUS WILL SMOOTH MORE THE 2ND DER.
  !     PXEXP0  : PTAUS IS WEIGHTED BY AN EXP(-((X-PXEXP0)/PXEXPDL)**2)
  !     PXEXPDL : IF PXEXP0 NOT IN [PXIN(1),PXIN(1)+PERIOD], EXP() IGNORED AND PTAUS=CST
  !     .         (SEE ROUTINE CBFITBND BELOW)
  !     KFLAG   : ERROR FLAG: IF NOT 0, THERE IS A PROBLEM
  !
  !-----------------------------------------------------------------------
  USE prec_rkind
  implicit none
  REAL(RKIND) :: EPSILON
  PARAMETER(EPSILON = 1.0E-10)
  !
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), PYINPP(KNIN), &
    &  PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
  !
  REAL(RKIND) :: PTAUS, PERIOD, PXEXP0, PXEXPDL, ZY, ZYP, ZYPP, ZYINT, ZDXMIN
  INTEGER :: I, KOPT, KNIN, KNOUT, KFLAG
  !-----------------------------------------------------------------------
  !
  !     0. CHECK INPUT CONSISTENCY
  !
  KFLAG=1
  IF (PTAUS .LT. 0) THEN
    ZDXMIN = MINVAL(PXIN(2:KNIN)-PXIN(1:KNIN-1))
    PTAUS = ABS(PTAUS) * ZDXMIN**3
    PRINT *,'% TAUS CHANGED TO DEFAULT VALUE = ',PTAUS
  END IF
  !
  CALL CBFITPER(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PTAUS,PERIOD,PXEXP0,PXEXPDL)
  !
  !L    2. COMPUTE INTERPOLATED VALUE AT EACH PXOUT
  !
  KFLAG=2
  IF (PTAUS .EQ. 0.0) THEN
    CALL SPLIPERA(PXIN,PYIN   ,PYINPP,KNIN,PXOUT,KNOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT, &
      &      PERIOD,KOPT)
  ELSE
    CALL SPLIPERA(PXIN,PYINNEW,PYINPP,KNIN,PXOUT,KNOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT, &
      &      PERIOD,KOPT)
  ENDIF
  !
  KFLAG=0
  RETURN
END SUBROUTINE CBSPLGNP
!-----------------------------------------------------------------------
SUBROUTINE CBFITBND(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PSIG, &
  &  NBCLFT,NBCRGT,XBCLFT,XBCRGT,YBCLFT, &
  &  YBCRGT,PXEXP0,PXEXPDL)
  !
  !     NON-PERIODIC B.C
  !
  !     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
  !     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
  !     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
  !     SIG = TAU*SIGMA_K/min(SIGMA_K) OF PAPER.
  !
  !     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
  !     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE spgbtrf_s.f
  !     (LAPACK SOURCE COPIED FROM NETLIB.ORG)
  !
  !     IF TAUS=0, PYINNEW NOT USED => PYINNEW(1) OR PYINNEW=PYIN IS OK
  !
  !     BOUNDARY CONDITIONS, 3 TYPES DETERMINED BY THE VALUE OF (NBCLFT/RGT):
  !
  !     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
  !     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
  !     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
  !
  !     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
  !
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
  !     FOR TYPE 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
  !
  !     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
  !     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
  !
  !     IF END POINTS ARE USED FOR THE B.C. AND TYPE 0 OR 1, THEN USE SYMMETRY
  !     OF MATRIX
  !
  !     IF XBCLFT OR XBCRGT ARE USED, IMPOSE B.C. ON NODE CLOSEST TO XBCLFT OR XBCRGT
  !
  !     TENSION TAUS(K) IS GIVEN WITH AN EXPONENTIAL FORM TO BE ABLE TO LOCALIZE
  !     IT:
  !     .     TAU_K = PTAUS * EXP( -COF * ((X-X0)/DX)**2)
  !
  !     WHERE X0 = PXEXP0 AND DX = PXEXPDL, AND:
  !     COF = 1. IF PXEXP0 IN [PXIN(1),PXIN(KNIN)], 0 OTHERWISE
  !     THUS SETTING PXEXP0 OUTSIDE DOMAIN GIVES A CST TAU_K VALUE
  !
  !     TOTAL DIMENSION OF PAMAT. THE REQUIRED SPACE DEPENDS IF
  !     .         PTAUS IS ZERO OR NOT AND ON THE B.C. (SYMMETRIC OR NOT)
  !     THUS, MDMATOT CAN VARY BETWEEN 2*KNIN AND 10*KNIN (MAX. VALUE NEEDED)
  !
  !-----------------------------------------------------------------------
  !
  USE prec_rkind
  implicit none
  integer :: knin, idamat, nbclft, nbcrgt, kpm2, ixbc, ibctyp, iik, &
    &  itauval, isym, n, i, j, k, iup, idown, idiag, ishift, ieff, ikp2 &
    &  ,ikp1, ikm1, ikm2, jk, jkp1, jkp2, jeff, iii, iupsofar, idwnsofa &
    &  , jbc, ik, idiamik, idiapik, iklft, irhs,  &
    &   jkm1, jkm2
  LAPACKINT INFO2,INFO
  DIMENSION :: KPM2(KNIN,-2:+2), IXBC(2), IBCTYP(2)
  !OS      pointer(iptr_ftauk,ftauk)
  REAL(RKIND), DIMENSION(KNIN) :: ftauk
  !
  REAL(RKIND) :: xbclft, xbcrgt, ybclft, ybcrgt, pxexp0, pxexpdl, pxin, &
    &  pyin, pyinnew, pyinpp, whk, wohk, psig
  DIMENSION :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), &
    &  PYINPP(KNIN), WHK(KNIN), WOHK(KNIN), &
    &  ZYBC(2), PSIG(KNIN)
  REAL(RKIND) :: zybc, xtkm1, xohkm1, xohkm2, xtk, xhkm1, &
    &  xohk, xtkp1, xhk, xohkp1, xykp1, xyk, &
    &  xykm1, ztaueff, zcofexp, zxexp0, zxexpdl, a1, a2, a3, a4, b1, &
    &  b2, b3, b4, px, fakk, fakkp1, fakkp2, frhs, zdelx, zero, zvalue, &
    &  zypeff, ztohkk1, zsign, fa2, fa3, fa0, fa1, &
    &  fakkm1, fakkm2, fun_ftauk, fcccc0, fcccc1, fcccc2, fcccc3
  integer nel
  REAL(RKIND), ALLOCATABLE :: PAMAT(:,:)
  !
  !
  !     FUNCTIONS FOR MATRIX COEFFICIENTS
  !
  REAL(RKIND) :: zsix, zthree, ztwo, zone
  PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
  FAKKM2(XTKM1,XOHKM1,XOHKM2) = XTKM1*XOHKM1*XOHKM2
  FAKKM1(XTK,XTKM1,XHKM1,XOHK,XOHKM1,XOHKM2) = XHKM1/zsix &
    &  - XOHKM1*(XTK*XOHK+(XTK+XTKM1)*XOHKM1 + XTKM1*XOHKM2)
  FAKK(XTKP1,XTK,XTKM1,XHK,XHKM1,XOHK,XOHKM1) = (XHK+XHKM1)/ZTHREE &
    &  + XOHK*XOHK*(XTKP1+XTK) &
    &  + XOHKM1*(ZTWO*XTK*XOHK+(XTK+XTKM1)*XOHKM1)
  FAKKP1(XTKP1,XTK,XHK,XOHKP1,XOHK,XOHKM1) = XHK/zsix &
    &  - XOHK*(XTKP1*XOHKP1+(XTK+XTKP1)*XOHK + XTK*XOHKM1)
  FAKKP2(XTKP1,XOHKP1,XOHK) = XTKP1*XOHKP1*XOHK
  !
  FRHS(XYKP1,XYK,XYKM1,XOHK,XOHKM1) = (XYKP1-XYK)*XOHK &
    &  - (XYK-XYKM1)*XOHKM1
  !
  !     WHEN ONE WANTS AN ARRAY FOR TAU*SIGMA_K**2, THEN ONE SHOULD REPLACE
  !     THE FUNCTION FTAU BY AN ARRAY
  !
  fun_FTAUK(IIK)= ZTAUEFF*EXP(-ZCOFEXP*(PXIN(IIK)-ZXEXP0)**2 &
    &  /ZXEXPDL**2)
  !%OS      FTAUK(IIK) = ZTAUEFF
  !
  !.......................................................................
  !*COMDECK CUCCCC
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         23.04.88            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
  ! ----------------------------------------------------------------------
  !
  FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
    &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
    &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
  FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
    &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
    &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / (B1-B2) - &
    &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
    &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        A1 - &
    &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
  ! -- FCCCC0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  ! ----------------------------------------------------------------------
  ! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  !.......................................................................
  !
  !-----------------------------------------------------------------------
  !
  !     0. INITIALIZATION
  !
  ITAUVAL = 1
  IF (PSIG(1) .EQ. 0._RKIND) ITAUVAL = 0
  ZTAUEFF = abs(PSIG(1))
  ZXEXP0 = PXEXP0
  ZXEXPDL = PXEXPDL
  ZCOFEXP = 1.0_RKIND
  IF (ZXEXP0.LT.PXIN(1) .OR. ZXEXP0.GT.PXIN(KNIN)) ZCOFEXP=0.0_RKIND
  !
  ISYM = 1
  IF (NBCLFT.GE.10 .OR. NBCRGT.GE.10 .OR. NBCLFT.EQ.2 &
    &  .OR. NBCRGT.EQ.2) ISYM = 0
  !
  !
  !     0.3 PREPARE BAND WIDTH
  !
  !     PREPARE MATRIX DIMENSION.
  !     MINIMUM REQUIRED:
  !     IUP = 1 = IDOWN; IF PTAUS.NE.0 => IUP = IDOWN = 2
  !     IF SYMMETRIC, USE ONLY UPPER BAND AND DIAGONAL =>IDAMAT=IUP+1
  !     IF ASYMMETRIC => IDAMAT = 2*IDOWN + IUP + 1
  !     IF B.C. NOT AT END OF INTERVAL => IUP = IUP + 1, AND/OR IDOWN=IDOWN+1
  !
  !     => ALTOGETHER, MINIMUM VALUE: IDOWN=1, IUP=1, SYM. =>IDAMAT_MAX = 2
  !     => ALTOGETHER, MAXIMUM VALUE: IDOWN=3, IUP=3 =>IDAMAT_MAX = 10
  !
  !
  IF (ITAUVAL .EQ. 0) THEN
    IUP   = 1
    IDOWN = 1
    IDAMAT = 2
  ELSE
    IUP   = 2
    IDOWN = 2
    IDAMAT = 3
  END IF
  IDIAG = IUP + 1
  IF (ISYM .EQ. 0) THEN
    IDAMAT = 3*(IDAMAT-1)+1
    IF (NBCLFT .GE. 10) THEN
      IUP = IUP + 1
      IDAMAT = IDAMAT + 1
    END IF
    IF (NBCRGT .GE. 10) THEN
      IDOWN = IDOWN + 1
      IDAMAT = IDAMAT + 2
    END IF
    IDIAG = IUP + IDOWN + 1
  ENDIF
  !
  ALLOCATE(PAMAT(IDAMAT,KNIN))
  !
  KNIN = KNIN
  PYINPP = 0._RKIND
  PAMAT = 0._RKIND
  !
  !     0.2 PRE-COMPUTE H_K AND 1./H_K, and zftauk
  !
  DO K=1,KNIN-1
    WHK(K)  = (PXIN(K+1) - PXIN(K))
  enddo
  DO K=1,KNIN-1
    WHK(K)  = (PXIN(K+1) - PXIN(K))
    WOHK(K) = zone / WHK(K)
    ftauk(k)=psig(k)
    !%OS        ftauk(k)=fun_ftauk(k)
  END DO
  WHK(KNIN) = 0.0_RKIND
  WOHK(KNIN) = 0.0_RKIND
  ftauk(knin)=psig(KNIN)
  !%OS      ftauk(knin)=fun_ftauk(KNIN)
  if (PSIG(1).lt.0._RKIND) ftauk(1)=-10._RKIND*ftauk(1)
  !
  !     0.4 DETERMINE NEIGHBOURS: K-2, K-1, .., K+2
  !     WHEN OUT OF BOUNDS, POINT TO INDEX KNIN, AS WHK, WOHK(KNIN)=0.0
  !
  DO ISHIFT=-2,+2
    DO K=1,KNIN
      KPM2(K,ISHIFT) = K + ISHIFT
    END DO
  END DO
  !     OUT OF INTERVAL: SEND TO KNIN
  KPM2(1,-2)   = KNIN
  KPM2(1,-1)   = KNIN
  KPM2(2,-2)   = KNIN
  KPM2(KNIN-1,+2) = KNIN
  KPM2(KNIN  ,+1) = KNIN
  KPM2(KNIN  ,+2) = KNIN
  !
  !     1. CONSTRUCT MATRIX AND R.H.S
  !     LAPACK SET-UP OF MATRIX:    A(I,J) -> PAMAT(I-J+IDIAG, J)
  !
  !.......................................................................
  !     AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART, THAT IS IF J.GE.I
  !
  DO K=1,KNIN-1
    IEFF = K + IDIAG
    IKP2 = KPM2(K,+2)
    IKP1 = KPM2(K,+1)
    IKM1 = KPM2(K,-1)
    IKM2 = KPM2(K,-2)
    !     A(K,K)
    JK = K
    PAMAT(IEFF-JK,JK) = FAKK(FTAUK(IKP1),FTAUK(K),FTAUK(IKM1),WHK(K) &
      &    ,WHK(IKM1),WOHK(K),WOHK(IKM1))
    !
    !     A(K,K+1)
    JKP1 = K + 1
    PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),FTAUK(K),WHK(K), &
      &    WOHK(IKP1),WOHK(K),WOHK(IKM1))
    !     A(K,K-1)
    !%OS        JKM1 = K - 1
    !%OS        PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(K),FTAUK(IKM1),WHK(IKM1), &
    !%OS     &    WOHK(K),WOHK(IKM1),WOHK(IKM2))
    !
    IF (ITAUVAL .EQ. 1) THEN
      !     A(K,K+2)
      JKP2 = K + 2
      IF (JKP2 .LE. KNIN) &
        &      PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
        &      WOHK(K))
      !     A(K,K-2)
      !%OS          JKM2 = K - 2
      !%OS          PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1),WOHK(IKM1), &
      !%OS     &      WOHK(IKM2))
    ENDIF
    !     B(K)
    PYINPP(K) = FRHS(PYIN(IKP1),PYIN(K),PYIN(IKM1),WOHK(K), &
      &    WOHK(IKM1))
    !
  END DO
  !
  !     2. BOUNDARY CONDITIONS
  !
  !     2.1 IF NON-SYMMETRIC, COPY TOP PART TO BOTTOM BEFORE APPLYING
  !     B.C.
  !
  IF (ISYM .EQ. 0) THEN
    DO I=1,KNIN
      IEFF = I + IDIAG
      DO J=I+1,MIN(I+MIN(IUP,IDOWN),KNIN)
        JEFF = J + IDIAG
        !     A(J,I) = A(I,J)
        PAMAT(JEFF-I,I) = PAMAT(IEFF-J,J)
      END DO
    END DO
  ENDIF
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,knin
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,knin)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,knin)
  !
  !%OS
  !
  !     2.2 B.C. AT TWO LOCATIONS PXIN(IXBC(JBC)), JBC=1,2
  !     IBCTYP(JBC) = 0, 1 OR 2 (TYPE OF B.C, SEE ABOVE).
  !     SO FAR USES NODE CLOSEST TO XBCLFT/RGT FOR LOCATION
  !     OF B.C., INSTEAD OF ACTUAL VALUE OF XBCLFT/RGT
  !
  IXBC(1) = 1
  IXBC(2) = KNIN
  IF (NBCLFT .GE. 10) THEN
    DO I=1,KNIN
      IF (PXIN(I) .GE. XBCLFT) GO TO 220
    END DO
220 CONTINUE
    ZDELX = ABS(PXIN(I)-XBCLFT)
    IXBC(1) = I
    IF (I .GE. KNIN) THEN
      IXBC(1) = KNIN
      PRINT *,' WARNING: LEFT B.C. AT I=KNIN: XBCLFT=',XBCLFT, &
        &      '  PXIN(KNIN)= ',PXIN(KNIN)
    ELSE IF (ABS(PXIN(I-1)-XBCLFT).LE.ZDELX .AND. I.NE.1) THEN
      IXBC(1) = I-1
    ENDIF
  ENDIF
  !
  IF (NBCRGT .GE. 10) THEN
    DO I=1,KNIN
      IF (PXIN(I) .GE. XBCRGT) GO TO 221
    END DO
221 CONTINUE
    ZDELX = ABS(PXIN(I)-XBCRGT)
    IXBC(2) = I
    IF (I .LE. 1) THEN
      IXBC(2) = 1
      PRINT *,' WARNING: RIGHT B.C. AT I=1: XBCRGT=',XBCRGT, &
        &      '  PXIN(1)= ',PXIN(1)
    ELSE IF (I .GT. KNIN) THEN
      IXBC(2) = KNIN
    ELSE IF (ABS(PXIN(I-1)-XBCRGT) .LE. ZDELX) THEN
      IXBC(2) = I-1
    ENDIF
  ENDIF
  !
  ZYBC(1) = YBCLFT
  ZYBC(2) = YBCRGT
  nel=10
  IBCTYP(1) = MOD(NBCLFT,nel)
  IBCTYP(2) = MOD(NBCRGT,nel)
  IF (IXBC(1) .EQ. IXBC(2)) THEN
    PRINT *,' ERROR, B.C. AT SAME LOCATIONS: IXBC(1)=IXBC(2)= ', &
      &    IXBC(1)
    !%OS        STOP '1=2'
    RETURN
  ELSE IF (IXBC(1) .GT. IXBC(2)) THEN
    PRINT *,' WARNING, NEEDED TO SWITCH B.C. POINTS AS IXBC(1)= ', &
      &    IXBC(1),' > IXBC(2)= ',IXBC(2)
    III = IXBC(1)
    IXBC(1) = IXBC(2)
    IXBC(2) = III
    ZYBC(1) = YBCRGT
    ZYBC(2) = YBCLFT
    nel=10
    IBCTYP(1) = MOD(NBCRGT,nel)
    IBCTYP(2) = MOD(NBCLFT,nel)
  ENDIF
  !
  !     2.3 MOVE EQUATIONS UP OR DOWN IF B.C. IS NOT AN END POINT
  !
  IF (IXBC(1) .NE. 1) THEN
    !
    !     MOVE ROW EQ. K=2,..,IXBC(1) UP BY ONE
    IUPSOFAR = IUP - 1
    DO K=2,IXBC(1)
      DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUPSOFAR)
        PAMAT(IDIAG+(K-1)-J,J) = PAMAT(IDIAG+K-J,J)
      END DO
      PYINPP(K-1) = PYINPP(K)
      !     ZERO A((K-1),(K-1)-IDOWN)
      IF (K-1-IDOWN .GE. 1) PAMAT(IDIAG+IDOWN,K-1-IDOWN) = 0.0_RKIND
    END DO
    !     ZERO ROW IXBC(1) AND RHS
    K = IXBC(1)
    DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUP)
      PAMAT(IDIAG+K-J,J) = 0.0_RKIND
    END DO
    PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  IF (IXBC(2) .NE. KNIN) THEN
    !     
    !     MOVE EQ. K=IXBC(2),..,KNIN-1 DOWN BY ONE
    IDWNSOFA = IDOWN - 1
    DO K=KNIN-1,IXBC(2),-1
      DO J=MAX(1,K-IDWNSOFA),MIN(KNIN,K+IUP)
        PAMAT(IDIAG+(K+1)-J,J) = PAMAT(IDIAG+K-J,J)
      END DO
      PYINPP(K+1) = PYINPP(K)
      !     ZERO A((K+1),(K+1)+IUP)
      IF (K+1+IUP .LE. KNIN) PAMAT(IDIAG-IUP,K+1+IUP) = 0.0_RKIND
    END DO
    !     ZERO ROW IXBC(2) AND RHS
    K = IXBC(2)
    DO J=MAX(1,K-IDOWN),MIN(KNIN,K+IUP)
      PAMAT(IDIAG+K-J,J) = 0.0_RKIND
    END DO
    PYINPP(K) = 0.0_RKIND
  ENDIF
  !
  !     2.4 FOR ROW=IXBC(), MODIFY MATRIX AND RHS ACCORDING TO B.C. TYPE
  !
  ZERO = 0.0_RKIND
  DO JBC=1,2
    IK = IXBC(JBC)
    ZVALUE = ZYBC(JBC)
    IEFF = IK + IDIAG
    IKP2 = KPM2(IK,+2)
    IKP1 = KPM2(IK,+1)
    IKM1 = KPM2(IK,-1)
    IKM2 = KPM2(IK,-2)
    IF (IBCTYP(JBC) .EQ. 0) THEN
      !
      !     SYMMETRIZE => COL IK GOES TO RIGHT-HAND SIDE AND THEN ZEROED
      !
      IF (ISYM .EQ. 1) THEN
        IDIAMIK = IDIAG - IK
        DO I=MAX(1,IK-IUP),IK-1
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
          PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
        END DO
        IDIAPIK = IDIAG + IK
        DO I=IK+1,MIN(KNIN,IK+IUP)
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(IDIAPIK-I,I)
          PAMAT(IDIAPIK-I,I) = 0.0_RKIND
        END DO
      ELSE
        IDIAMIK = IDIAG - IK
        DO I=MAX(1,IK-IUP),MIN(KNIN,IK+IDOWN)
          PYINPP(I) = PYINPP(I) - ZVALUE * PAMAT(I+IDIAMIK,IK)
          PAMAT(I+IDIAMIK,IK) = 0.0_RKIND
        END DO
        !     ZERO ROW IK
        DO J=MAX(1,IK-IDOWN),MIN(KNIN,IK+IUP)
          PAMAT(IEFF-J,J) = 0.0_RKIND
        END DO
      ENDIF
      !
      !     REPLACE ROW IK BY EQUATION: G_K = ZVALUE
      PAMAT(IDIAG,IK) = 1.0_RKIND
      PYINPP(IK) = ZVALUE
      !
    ELSE IF (IBCTYP(JBC) .EQ. 1) THEN
      !
      !     1ST DERIVATIVE GIVEN
      !
      ZYPEFF = ZVALUE
      IF (ZVALUE .GT. 1.E+31_RKIND) THEN
        !     FROM LGRANGIAN INTERPOLATION
        IKLFT = IK - 1
        IF (IK .EQ. 1) IKLFT = IK
        IF (IKLFT+3 .GT. KNIN) IKLFT = KNIN - 3
        ZYPEFF = FCCCC1(PYIN(IKLFT),PYIN(IKLFT+1),PYIN(IKLFT+2), &
          &        PYIN(IKLFT+3),PXIN(IKLFT),PXIN(IKLFT+1),PXIN(IKLFT+2), &
          &        PXIN(IKLFT+3),PXIN(IK))
      ELSE IF (ZVALUE .LT. -1.E+31_RKIND) THEN
        IKLFT = IK
        IF (IK .EQ. KNIN) IKLFT = IK - 1
        ZYPEFF = (PYIN(IKLFT+1)-PYIN(IKLFT)) &
          &        / (PXIN(IKLFT+1)-PXIN(IKLFT))
      ENDIF
      ZTOHKK1 = FTAUK(IK)*WOHK(IK)*WOHK(IKM1)
      !     A(IK,IK)
      IF (IK .NE. KNIN) PAMAT(IEFF-IK,IK) = FAKK(FTAUK(IKP1),FTAUK(IK), &
        &      ZERO,WHK(IK),ZERO,WOHK(IK),ZERO) + ZTOHKK1
      IF (IK .EQ. KNIN) PAMAT(IEFF-IK,IK) = FAKK(ZERO,FTAUK(IK), &
        &      FTAUK(IKM1),ZERO,WHK(IKM1),ZERO,WOHK(IKM1)) + ZTOHKK1
      !     A(IK,IK-1)
      JKM1 = IK - 1
      IF (ISYM.EQ.0 .AND. JKM1.GE.1) THEN
        IF (IK .NE. KNIN) PAMAT(IEFF-JKM1,JKM1) = - ZTOHKK1
        IF (IK .EQ. KNIN) PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(IK), &
          &        FTAUK(IKM1),WHK(IKM1),ZERO,WOHK(IKM1),WOHK(IKM2))
      ENDIF
      !     A(IK,IK+1)
      JKP1 = IK + 1
      IF (JKP1 .LE. KNIN) &
        &      PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1), &
        &      FTAUK(IK),WHK(IK),WOHK(IKP1),WOHK(IK),ZERO)
      !
      IF (ITAUVAL .EQ. 1) THEN
        !     A(IK,IK+2)
        JKP2 = IK + 2
        IF (JKP2 .LE. KNIN) &
          &        PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1), &
          &        WOHK(IK))
        !     A(IK,IK-2)
        JKM2 = IK - 2
        IF (ISYM.EQ.0 .AND. JKM2.GE.1) THEN
          IF (IK .NE. KNIN) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
          IF (IK .EQ. KNIN) PAMAT(IEFF-JKM2,JKM2) = FAKKM2(FTAUK(IKM1), &
            &          WOHK(IKM1),WOHK(IKM2))
        ENDIF
      ENDIF
      !     RHS
      ZSIGN = -1._RKIND
      IF (IK .EQ. KNIN) ZSIGN = +1._RKIND
      IF (IK .NE. KNIN) PYINPP(IK) = FRHS(PYIN(IKP1),PYIN(IK),ZERO, &
        &      WOHK(IK),ZERO) - ZYPEFF
      IF (IK .EQ. KNIN) PYINPP(IK) = FRHS(ZERO,PYIN(IK),PYIN(IKM1), &
        &      ZERO,WOHK(IKM1)) + ZYPEFF
      !
    ELSE IF (IBCTYP(JBC) .EQ. 2) THEN
      !
      !     FUNCTION IS GIVEN
      !
      !     A(IK,IK)
      PAMAT(IEFF-IK,IK) = - FTAUK(IK) * (WOHK(IK) + WOHK(IKM1))
      !     A(IK,IK+1)
      JKP1 = IK + 1
      IF (JKP1 .LE. KNIN) &
        &      PAMAT(IEFF-JKP1,JKP1) = FTAUK(IK) * WOHK(IK)
      !     A(IK,IK-1)
      JKM1 = IK - 1
      IF (ISYM.EQ.0 .AND. JKM1.GE.1) &
        &      PAMAT(IEFF-JKM1,JKM1) = FTAUK(IK) * WOHK(IKM1)
      !
      IF (ITAUVAL .EQ. 1) THEN
        !     A(IK,IK+2)
        JKP2 = IK + 2
        IF (JKP2 .LE. KNIN) PAMAT(IEFF-JKP2,JKP2) = 0.0_RKIND
        !     A(IK,IK-2)
        JKM2 = IK - 2
        IF (ISYM.EQ.0 .AND. JKM2.GE.1) PAMAT(IEFF-JKM2,JKM2) = 0.0_RKIND
      ENDIF
      !     RHS
      PYINPP(IK) = PYIN(IK) - ZVALUE
      !
    ENDIF
    !
  END DO
  !
  !     3. SOLVE SYSTEM
  !
  !     USE INTEGER WORK SPACE FROM KPM2(0) ARRAY FOR IPIVOT,
  !     AS KPM2(K,0) NOT NEEDED NOR USED
  !
  !%OS
  !     debug, print matrix and rhs
  !%OS      write(6,'(3a4,a)') 'i ','j1 ','j2 ',' i,j1  i,j1+1,..... i,j2'
  !%OS      do i=1,knin
  !%OS        ieff = idiag + i
  !%OS        j1 = i-idown
  !%OS        j2 = i+iup
  !%OSc%OS        j1 = max(i-idown,1)
  !%OSc%OS        j2 = min(i+iup,knin)
  !%OS        write(6,'(3i4,1p10e13.4)') i,j1,j2,(pamat(ieff-j,j),j=j1,j2)
  !%OS      end do
  !%OS      write(6,'(a4,a12)') 'i','RHS'
  !%OS      write(6,'(i4,1pe13.4)') (i,pyinpp(i),i=1,knin)
  !
  !%OS
  !%OS      print *,'isym= ',isym
  IF (ISYM .EQ. 1) THEN
    !   Single precision (using routines in spgbtrf_s.f compiled with -r8 for example)
    !   Double precision (using the lapack libraries from the compiler)
    CALL DPBTRF('U',KNIN,IUP,PAMAT,IDAMAT,INFO)
  ELSE
    CALL DGBTRF(KNIN,KNIN,IDOWN,IUP,PAMAT,IDAMAT,KPM2(1,0),INFO)
  ENDIF
  IRHS = 1
  IF (INFO .EQ. 0) THEN
    IF (ISYM .EQ. 1) THEN
      CALL DPBTRS('U',KNIN,IUP,IRHS,PAMAT,IDAMAT,PYINPP,KNIN,INFO2)
    ELSE
      CALL DGBTRS('N',KNIN,IDOWN,IUP,IRHS,PAMAT,IDAMAT,KPM2(1,0),PYINPP, &
        &      KNIN,INFO2)
    ENDIF
  ELSE
    PRINT *,' ERROR IN SP/GBTRF: INFO = ',INFO
    !%OS        STOP 'INFO'
    RETURN
  ENDIF
  !
  !     4. COMPUTE NEW VALUES OF Y_K (NON-STANDARD CUBIC SPLINE ONLY)
  !
  IF (ITAUVAL .EQ. 1) THEN
    DO K=1,KNIN
      IKP1 = KPM2(K,+1)
      IKM1 = KPM2(K,-1)
      PYINNEW(K) = PYIN(K) - FTAUK(K) * &
        &      ((PYINPP(IKP1)-PYINPP(K))*WOHK(K) &
        &      - (PYINPP(K)-PYINPP(IKM1))*WOHK(IKM1))
    END DO
    !
  ENDIF

  IF (INFO2 .LT. 0) THEN
    PRINT *,' ERROR IN SP/GBTRS: INFO2 = ',INFO2
    !%OS        STOP 'INFO2'
    RETURN
  ENDIF
  !
  RETURN
END SUBROUTINE CBFITBND
SUBROUTINE CBFITPER(PXIN,PYIN,PYINNEW,KNIN,PYINPP,PTAUS,PERIOD,PXEXP0,PXEXPDL)
  !
  !     PERIODIC B.C WITH Y(PXIN(1)+PERIOD) = PYIN(1)
  !
  !     PREPARE SECOND DERIVATIVE OF CUBIC SPLINE INTERPOLATION AND NEW
  !     VALUES OF Y AT NODES YNEW FITTED SUCH THAT CHI**2 + TAUS*F''**2
  !     IS MINIMIZED ACCORDING TO HIRSHMAN ET AL, PHYS. PLASMAS 1 (1994) 2280.
  !     TAUS = TAU*SIGMA_K OF PAPER ASSUMING SIGMA_K CONSTANT.
  !
  !     SETTING TAUS=0., ONE FINDS THE USUAL CUBIC SPLINE INT. WITH CHI**2=0
  !     TAUS LARGE => FIT CLOSER TO STRAIGHT LINE (SECOND DERIV.=0)
  !
  !     IF LAPACK ROUTINES NOT AVAILABLE, USE NONSYM.F AND REMOVE "c%nonsym"
  !     (COPIED SOURCE FROM NETLIB.COM)
  !
  !     IF TAUS=0, PYINNEW NOT USED => PYINNEW(1) OR PYINNEW=PYIN IS OK
  !
  !
  USE prec_rkind
  implicit none
  !
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN), PYINNEW(KNIN), &
    &  PYINPP(KNIN), WHK(KNIN), WOHK(KNIN)
  REAL(RKIND), ALLOCATABLE :: PAMAT(:,:)

  INTEGER :: KTONUM(KNIN), KPM2(KNIN,-2:+2)
  INTEGER INFO2, N, IDOWN, IUP, IRHS
  !
  REAL(RKIND) :: FAKKM2, FAKKM1, FAKK, FAKKP1, FAKKP2, FRHS, FTAUK, ZTAUEFF
  REAL(RKIND) :: PTAUS, PXEXP0, PXEXPDL, ZXEXP0, ZXEXPDL, ZCOFEXP
  REAL(RKIND) :: XTKM1, XTK, XTKP1, XOHKP1, XOHKM1, XHKM1, XHK, XOHKM2, XOHK, XYK, XYKP1, XYKM1, PERIOD
  !
  INTEGER :: ICUBSTD, KNIN, I, J, K, IBAND, IDIAG, &
    &  IDAMAT, IHALF, ISHIFT, IEFF, IKP2, IKP1, IKM1, IKM2, JKM2, JKM1, JK, &
    &  JKP1, JKP2, II
  !
  !%nonsym
!!$      REAL(RKIND) :: ANONSYM(501*9)
!!$      INTEGER :: IPIVOT(501)
  !%nonsym
  CHARACTER*1 ZCHAR
  !
  !     FUNCTIONS FOR MATRIX COEFFICIENTS
  !
  FAKKM2(XTKM1,XOHKM1,XOHKM2) = XTKM1*XOHKM1*XOHKM2
  FAKKM1(XTK,XTKM1,XHKM1,XOHK,XOHKM1,XOHKM2) = XHKM1/6. &
    &  - XOHKM1*(XTK*XOHK+(XTK+XTKM1)*XOHKM1 + XTKM1*XOHKM2)
  FAKK(XTKP1,XTK,XTKM1,XHK,XHKM1,XOHK,XOHKM1) = (XHK+XHKM1)/3. &
    &  +XOHK*XOHK*(XTKP1+XTK) + XOHKM1*(2.*XTK*XOHK+(XTK+XTKM1)*XOHKM1)
  FAKKP1(XTKP1,XTK,XHK,XOHKP1,XOHK,XOHKM1) = XHK/6. &
    &  - XOHK*(XTKP1*XOHKP1+(XTK+XTKP1)*XOHK + XTK*XOHKM1)
  FAKKP2(XTKP1,XOHKP1,XOHK) = XTKP1*XOHKP1*XOHK
  !
  FRHS(XYKP1,XYK,XYKM1,XOHK,XOHKM1) = (XYKP1-XYK)*XOHK &
    &  - (XYK-XYKM1)*XOHKM1
  !
  !     WHEN ONE WANTS AN ARRAY FOR TAU*SIGMA_K**2, THEN ONE SHOULD REPLACE
  !     THE FUNCTION FTAU BY AN ARRAY
  !
  !%OS      FTAUK(IIK)= ZTAUEFF*EXP(-ZCOFEXP*(PXIN(IIK)-ZXEXP0)**2/ZXEXPDL**2)
  FTAUK(II) = ZTAUEFF
  !
  !-----------------------------------------------------------------------
  !
  !     0. INITIALIZATION
  !
  ICUBSTD = 0
  IF (PTAUS .EQ. 0._RKIND) ICUBSTD = 1
  !
  ZTAUEFF = PTAUS
  ZXEXP0 = PXEXP0
  ZXEXPDL = PXEXPDL
  ZCOFEXP = 1.0_RKIND
  IF (ZXEXP0.LT.PXIN(1) .OR. ZXEXP0.GT.PXIN(1)+PERIOD) ZCOFEXP=0.0_RKIND
  !
  !.......................................................................
  !
  !     0.3 PREPARE BAND WIDTH
  !
  IUP   = 4
  IDOWN = 4
  IDAMAT = 13
  IF (ICUBSTD .EQ. 1) THEN
    IUP   = 2
    IDOWN = 2
    IDAMAT = 7
  END IF
  IBAND = IDOWN + IUP + 1
  IDIAG = IBAND
  ALLOCATE(PAMAT(IDAMAT,KNIN))
  !%nonsym
  IDIAG = IBAND
  !%nonsym
  !
  PYINPP = 0.0_RKIND
  PAMAT = 0.0_RKIND
  !
  !     0.2 PRE-COMPUTE H_K AND 1./H_K
  !
  DO K=1,KNIN-1
    WHK(K)  = (PXIN(K+1) - PXIN(K))
    WOHK(K) = 1. / WHK(K)
  END DO
  WHK(KNIN) = PXIN(1) + PERIOD - PXIN(KNIN)
  IF (WHK(KNIN) .EQ. 0.0_RKIND) THEN
    PRINT *,' PXIN(1) + PERIOD - PXIN(KNIN) = 0. IN CBFITPER'
    PRINT *,' DO NOT INCLUDE EXTRA PERIODIC POINT: KNIN .NE. 1', &
      &    ' => CHANGE KNIN TO KNIN-1 ?'
!!$        STOP 'CBFITPER'
    return
  ELSE
    WOHK(KNIN) = 1. / WHK(KNIN)
  ENDIF
  !
  !   0.4 INDEXING OF UNKNOWNS IN ALTERNATIVE UP,DOWN,UP,... WAY IN ORDER
  !   TO KEEP A BAND MATRIX, EVEN IF IT DOUBLES THE BAND-WIDTH:
  !   K -> UNKNOWN_NUMBER:
  !   . KTONUM(K) = INUM WITH INUM=1 FOR K=1 ; =2*(K-1) FOR K IN [2,IHALF];
  !   . INUM=2*(KNIN-K)+3 FOR K IN [IHALF+1,KNIN]; WITH IHALF=KNIN/2 + 1
  !   WARNING: CARE IF KNIN IS ODD
  !
  IHALF = KNIN/2 + 1
  KTONUM(1) = 1
  DO K=2,IHALF
    KTONUM(K) = 2*K - 2
  END DO
  DO K=IHALF+1,KNIN
    KTONUM(K) = 2*(KNIN-K) + 3
  END DO
  !
  !     0.5 DETERMINE NEIGHBOURS: K-2, K-1, .., K+2
  !
  DO ISHIFT=-2,+2
    DO K=1,KNIN
      KPM2(K,ISHIFT) = K + ISHIFT
    END DO
  END DO
  !     PERIODIC ENDS
  KPM2(1,-2)   = KNIN - 1
  KPM2(1,-1)   = KNIN
  KPM2(2,-2)   = KNIN
  KPM2(KNIN-1,+2) = 1
  KPM2(KNIN  ,+1) = 1
  KPM2(KNIN  ,+2) = 2
  !
  !     1. CONSTRUCT MATRIX AND R.H.S
  !     LAPACK SET-UP OF MATRIX:    A(I,J) -> AMAT(I-J+IDIAG, J)
  !
  !.......................................................................
  !     (AS MATRIX SYMMETRIC, COMPUTE ONLY UPPER PART) NOT YET
  !
  DO K=1,KNIN
    IEFF = KTONUM(K) + IDIAG
    IKP2 = KPM2(K,+2)
    IKP1 = KPM2(K,+1)
    IKM1 = KPM2(K,-1)
    IKM2 = KPM2(K,-2)
    !     A(K,K-2)
    JKM2 = KTONUM(IKM2)
    PAMAT(IEFF-JKM2,JKM2) =FAKKM2(FTAUK(IKM1),WOHK(IKM1),WOHK(IKM2))
    !     A(K,K-1)
    JKM1 = KTONUM(IKM1)
    PAMAT(IEFF-JKM1,JKM1) = FAKKM1(FTAUK(K),FTAUK(IKM1),WHK(IKM1), &
      &    WOHK(K),WOHK(IKM1),WOHK(IKM2))
    !     A(K,K)
    JK = KTONUM(K)
    PAMAT(IEFF-JK,JK) =FAKK(FTAUK(IKP1),FTAUK(K),FTAUK(IKM1),WHK(K), &
      &    WHK(IKM1),WOHK(K),WOHK(IKM1))
    !     A(K,K+1)
    JKP1 = KTONUM(IKP1)
    PAMAT(IEFF-JKP1,JKP1) = FAKKP1(FTAUK(IKP1),FTAUK(K),WHK(K), &
      &    WOHK(IKP1),WOHK(K),WOHK(IKM1))
    !     A(K,K+2)
    JKP2 = KTONUM(IKP2)
    PAMAT(IEFF-JKP2,JKP2) = FAKKP2(FTAUK(IKP1),WOHK(IKP1),WOHK(K))
    !     B(K)
    PYINPP(KTONUM(K)) = FRHS(PYIN(IKP1),PYIN(K),PYIN(IKM1),WOHK(K) &
      &    ,WOHK(IKM1))
    !
  END DO
  !
  !     2. SOLVE SYSTEM
  !
  !     USE INTEGER WORK SPACE FROM KPM2(0) ARRAY FOR IPIVOT,
  !     AS KPM2)K,0) NOT NEEDED NOR USED
  !
  IRHS = 1
  INFO2 = 0
  !   Single precision (using routines in spgbtrf_s.f compiled with -r8 for example)
  !%single      CALL SGBTRF(N,KNIN,IDOWN,IUP,PAMAT,IDAMAT,KPM2(1,0),INFO2)
  !   Double precision (using the lapack libraries from the compiler)
  CALL DGBTRF(KNIN,KNIN,IDOWN,IUP,PAMAT,IDAMAT,KPM2(1,0),INFO2)
  IF (INFO2 .EQ. 0) THEN
    !   single precision with S as first character
    !   double precision with D as first character
    CALL DGBTRS('N',KNIN,IDOWN,IUP,IRHS,PAMAT,IDAMAT,KPM2(1,0),PYINPP, &
      &    KNIN,INFO2)
  ELSE
    PRINT *,' ERROR IN SP/GBTRF: INFO = ',INFO2
    STOP 'INFO'
  ENDIF
  !
  !%nonsym
!!$      DO I=1,IBAND*KNIN
!!$        ANONSYM(I) = 0.0_RKIND
!!$      ENDDO
!!$      DO I=1,KNIN
!!$!     UPPER PART
!!$        DO J=max(1,i),min(i+IUP,n)
!!$          ANONSYM((I-1)*IBAND+J-I+IUP+1) = PAMAT(IDIAG+I-J,J)
!!$        ENDDO
!!$!     LOWER PART
!!$        DO J=max(1,i-IUP),min(i-1,n-1)
!!$          ANONSYM((I-1)*IBAND+J-I+IUP+1) = PAMAT(IDIAG+J-I,I)
!!$        ENDDO
!!$      ENDDO
!!$      CALL NONSYM(ANONSYM,PAMAT,PYINPP,KNIN,IUP,IUP,1.0E-06_RKIND,INFO2)
  !%nonsym
  !
  !     2.2 COPY SOLUTION BACK TO STANDARD NUMBERING, USE PAMAT AS WORK SPACE
  !
  DO K=1,KNIN
    PAMAT(1,K) = PYINPP(KTONUM(K))
  END DO
  DO K=1,KNIN
    PYINPP(K) = PAMAT(1,K)
  END DO
  !
  !     3. COMPUTE NEW VALUES OF Y_K (NON-STANDARD CUBIC SPLINE ONLY)
  !
  IF (ICUBSTD .EQ. 0) THEN
    !
    DO K=1,KNIN
      IKP1 = KPM2(K,+1)
      IKM1 = KPM2(K,-1)
      PYINNEW(K) = PYIN(K) &
        &      - FTAUK(K) * ((PYINPP(IKP1)-PYINPP(K))*WOHK(K) &
        &      - (PYINPP(K)-PYINPP(IKM1))*WOHK(IKM1))
    END DO
    !
  ENDIF

  IF (INFO2 .LT. 0) THEN
    PRINT *,' ERROR IN SP/GBTRS: INFO2 = ',INFO2
  ENDIF
  !
  DEALLOCATE(pamat)
  RETURN
END SUBROUTINE CBFITPER
SUBROUTINE SPLIPERA(PXIN,PYIN,PYINPP,KNIN,PXOUT,KNOUT,PY,PYP,PYPP,PYINT,PERIOD,KOPT)
  USE prec_rkind
  implicit none
  REAL(RKIND) :: PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN)
  REAL(RKIND) :: XA(KNIN+1),YA(KNIN+1),Y2A(KNIN+1), YINT(KNIN+1)
  REAL(RKIND) :: PXOUT(KNOUT), PY(KNOUT), PYP(KNOUT), PYPP(KNOUT), PYINT(KNOUT), PERIOD
  INTEGER :: KNIN, KNOUT, IIN_XOUT(KNOUT), IPER(KNOUT), KOPT
  !
  REAL(RKIND) :: ZXSHIFT(KNOUT), XAKLO, XAKHI, YAKLO, YAKHI, Y2AKLO, Y2AKHI, H, A, B, A2, B2
  INTEGER :: KLO, KHI, K, J, I, IPREVIOUS
  !
  !   with periodic boundary conditions, SET
  !   XA(N+1) = PXIN(1) + PERIOD
  !   YA(N+1) = PYIN(1)
  !   Y2A(N+1) = PYINPP(1)
  !
  !   Calculates Y, Yprime, Yprimeprime, int(y)
  !
  !-----------------------------------------------------------------------
  ! COMPUTE IIN_XOUT,IPER,ZXSHIFT TO LOCATE PXOUT WITHIN PXIN
  !
  !  CALL FINDINDICESPER(PXIN,KNIN,PXOUT,KNOUT,IIN_XOUT,IPER,ZXSHIFT,PERIOD)
  XA(1:KNIN) = PXIN
  XA(KNIN+1) = PXIN(1) + PERIOD
  YA(1:KNIN) = PYIN
  YA(KNIN+1) = PYIN(1)
  Y2A(1:KNIN) = PYINPP
  Y2A(KNIN+1) = PYINPP(1)
  IF (KOPT .GE. 3) THEN
    ! PREPARE INTEGRAL UP TO XIN(I) TO ADD UP TO LOCAL VALUE
    ! SINCE IT TAKES SOME EXTRA TIME DO IT ONLY IF REQUIRED
    YINT(1) = 0._RKIND
    DO I=2,KNIN+1
      H = XA(I) - XA(I-1)
      YINT(I) = YINT(I-1) +  H/2._RKIND*(YA(I-1) + YA(I) &
        & - H**2/12._RKIND*(Y2A(I-1)+Y2A(I)))
    END DO
  END IF
  !
  IPREVIOUS = 1
  DO J=1,KNOUT
    !     SHIFT PXOUT POINTS WITHIN [PXIN(1),PXIN(KNIN)] WITH IPER * PERIOD
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      IPER(J) = INT((PXIN(1)-PXOUT(J))/PERIOD) + 1
      ZXSHIFT(J) = PXOUT(J) + IPER(J)*PERIOD
    ELSE IF (PXOUT(J) .GT. PXIN(1)+PERIOD) THEN
      IPER(J) = - (INT((PXOUT(J)-PXIN(1)-PERIOD)/PERIOD) + 1)
      ZXSHIFT(J) = PXOUT(J) + IPER(J)*PERIOD
    ELSE
      IPER(J) = 0
      ZXSHIFT(J) = PXOUT(J)
    ENDIF
    ! FIND INTERVAL IN PXIN WHERE PXOUT+-IPER(J)*PERIOD BELONGS, START FROM PREVIOUS INTERVAL, EXCEPT IF THERE WAS A JUMP BACK
    IF (ZXSHIFT(J) .LT. PXIN(IPREVIOUS)) IPREVIOUS = 1
    DO I=IPREVIOUS,KNIN
      IF ((ZXSHIFT(J) .GE. XA(I)) .AND. (ZXSHIFT(J) .LE. XA(I+1))) THEN
        IIN_XOUT(J) = I
        IPREVIOUS = I
        EXIT
      END IF
    END DO
  END DO
  !
  !
  DO J=1,KNOUT
    H = XA(IIN_XOUT(J)+1) - XA(IIN_XOUT(J))
    YAKLO = YA(IIN_XOUT(J))
    YAKHI = YA(IIN_XOUT(J)+1)
    Y2AKLO = Y2A(IIN_XOUT(J))
    Y2AKHI = Y2A(IIN_XOUT(J)+1)
    IF (H.EQ.0._RKIND) STOP 'BAD XA INPUT.'
    A=(XA(IIN_XOUT(J)+1) - ZXSHIFT(J)) / H
    B=1._RKIND - A
    A2=A**2
    B2=B**2
    PY(J)=A*YAKLO+B*YAKHI+ &
      &      (A*(A2-1._RKIND)*Y2AKLO+B*(B2-1._RKIND)*Y2AKHI)*(H**2)/6._RKIND
    IF (KOPT .GE. 1) PYP(J)=(YAKHI-YAKLO)/H - &
      &  ( (3._RKIND*A2-1._RKIND)*Y2AKLO - (3._RKIND*B2-1._RKIND)*Y2AKHI )*H/6._RKIND
    IF (KOPT .GE. 2) PYPP(J)=A*Y2AKLO+B*Y2AKHI
    IF (KOPT .GE. 3) PYINT(J)= YINT(IIN_XOUT(J)) + H/2._RKIND*(YAKLO*(1._RKIND-A2) + YAKHI*B2 &
      &  - (Y2AKLO*(1._RKIND-A2)**2 + Y2AKHI*B2*(2._RKIND-B2))  * H**2/12._RKIND) - IPER(J)*YINT(KNIN+1)
  END DO
  !
  RETURN
END SUBROUTINE SPLIPERA
!-----------------------------------------------------------------------
SUBROUTINE SPLIBNDA(PXIN,PYIN,PYINPP,KNIN,PXOUT,PY,PYP,PYPP,PYINT,KNOUT,KOPTXPOL,KOPTDER)
  USE prec_rkind
  implicit none
  REAL(RKIND) :: zsix, zthree, ztwo, zone
  PARAMETER(zsix=6._RKIND, zthree=3._RKIND, ztwo=2._RKIND, zone=1._RKIND)
  REAL(RKIND) :: PXIN(KNIN),PYIN(KNIN),PYINPP(KNIN)
  REAL(RKIND) :: PXOUT(KNOUT), PY(KNOUT), PYP(KNOUT), PYPP(KNOUT), PYINT(KNOUT)
  integer KNIN, KNOUT, KOPTXPOL, KOPTDER
  !
  !     KOPTDER = 0: COMPUTE ONLY FUNCTION PY POINTS PXOUT
  !     KOPTDER = 1: COMPUTE ALSO 1ST DERIVATIVE PYP
  !     KOPTDER = 2: COMPUTE ALSO 2ND DERIVATIVE PYPP
  !     KOPTDER = 3: COMPUTE ALSO INTEGRAL OF Y FROM XIN(1) TO XOUT(J) IN PYINT(J)
  !
  !   ABS(KOPTXPOL):
  !     KOPTXPOL = 0: STOP WITH ERROR MESSAGE IF OUT OF BOUND
  !     KOPTXPOL = 1: LINEAR EXTRAPOLATION: +1 FROM YEDGE,YPEDGE ; -1: FROM Y(LAST TWO POINTS)
  !     KOPTXPOL = 10: CONSTANT EXTRAPOLATION: +10: Y=YEDGE, -10: Y=0. OUTSIDE PXIN INTERVAL
  !     KOPTXPOL = 2: USE QUADRATIC EXTRAPOLATION IF X OUT OF BOUND
  !     KOPTXPOL = 3: USE CUBIC EXTRAPOLATION IF X OUT OF BOUND
  !     KOPTXPOL = 21: USE QUADRATIC WITHIN ALFA*DELTA_X AND LINEAR FURTHER
  !     KOPTXPOL = 31: USE CUBIC WITHIN ALFA*DELTA_X AND LINEAR    FURTHER
  !     KOPTXPOL = 32: USE CUBIC WITHIN ALFA*DELTA_X AND QUADRATIC FURTHER
  !
  !     KOPTXPOL > 0: VALUE AND 1ST DER. CONTINUOUS AT END OF INTERVAL, THUS
  !     .             USES CUBIC SPLINE OF LAST INTERVAL TO CONTINUE
  !     KOPTXPOL < 0: ONLY Y VALUE CONTINUOUS AND USES VALUES AT LAST BUT ONE,
  !     .             TWO, THREE POINTS TO EXTRAPOLATE (BETTER IF DER. AT EDGE
  !     .             IS WILD)
  !
  !-----------------------------------------------------------------------
  ! LOCAL VARIABLES:
  REAL(RKIND) :: ALFA
  PARAMETER(ALFA = 1._RKIND)
  REAL(RKIND) :: H, H2, A, B,  &
    & ZX1, ZY1, ZY1P, zx2, zy2, ZY2P, &
    & ZXN, ZYN, ZYNP, ZXNN, ZYNN, ZYNNP, &
    & ZXLFTDEL, ZYLFTDEL, ZYPLFTDEL, ZYINTLFTDL, &
    & ZXRGTDEL, ZYRGTDEL, ZYPRGTDEL, ZYINTRGTDL, &
    & ZYINT, ZYINT_XIN(KNIN), ZYIN_XPOL
  INTEGER ICONTDER
  INTEGER :: J1ST_XIN(KNIN),JLAST_XIN(KNIN), JLEFT(KNOUT), JRIGHT(KNOUT), &
    & JLEFT_DEL(KNOUT), JRIGHT_DEL(KNOUT), IOPTXPOL, I, J, K, KLO, KHI
  !
  ! VARIABLES RELATED TO FUNCTIONS:
  REAL(RKIND) :: FC3, X1, F1, P1, X2, &
    &  F2, P2, FC2, FC1, FC0, FQQQ0, FQQQ1, FQQQ2, &
    &  FLINEAR, FLINEARP, FCCCC0, FCCCC1, FCCCC2, FCCCC3, FQDQ0, FQDQ1, &
    &  FQDQ2, FCDCD0, FCDCD1, FCDCD2, FCDCD3, FB1, &
    &  FB2, FA2, FA3, FD2, FD1, &
    &  FCCCCM1, FCDCDM1, FQQQM1, FQDQM1, FLINEARM1, FLINXYP, FLINXYPM1
  REAL(RKIND) :: A1, A2, A3, A4, B1, B2, B3, B4, PX
  REAL(RKIND) :: FB0, FD0, FA0, FA1
  !
  !.......................................................................
  !*COMDECK CUCCCC
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         23.04.88            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
  ! ----------------------------------------------------------------------
  !
  FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
    &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
    &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
  FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
    &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
    &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        (A1-A2) / (B1-B2) - &
    &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
    &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
    &        A1 - &
    &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
  ! -- FCCCC0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
  ! ----------------------------------------------------------------------
  ! -- FCCCC1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCCCC1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCCCC1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &              FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &              PX * (ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &                    ZTHREE * PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  ! ----------------------------------------------------------------------
  ! -- FCCCC2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCCCC2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCCCC2(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &             ZTWO * FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &             zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4) * PX
  ! ----------------------------------------------------------------------
  ! -- FCCCC3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:     -
  ! -- FCCCC3(......,PX) = D3F/DX3 (PX)                                  -
  ! ----------------------------------------------------------------------
  FCCCC3(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &                      zsix * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
  ! ----------------------------------------------------------------------
  ! -- FCCCCM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM B1 TO PX:     -
  FCCCCM1(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
    &  (PX-B1)*(FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &  0.5_RKIND*(PX+B1)*FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
    &  FA2(A1,A2,A3,A4,B1,B2,B3,B4)/3._RKIND*(PX*(PX+B1)+B1*B1) + &
    &  0.25_RKIND*(PX+B1)*(PX*PX+B1*B1)*FA3(A1,A2,A3,A4,B1,B2,B3,B4))
  !-----------------------------------------------------------------------
  !.......................................................................
  !*COMDECK CUCDCD
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
  ! -- THE SIX ARGUMENTS X1,F1,P1,X2,F2,P2 ARE DEFINED AS FOLLOWS:      --
  ! -- F(X1) = F1 , F(X2) = F2 , DF/DX(X1) = P1 , DF/DX(X2) = P2        --
  ! ----------------------------------------------------------------------
  !
  FC3(X1,F1,P1,X2,F2,P2) = &
    &      (ZTWO * (F2 - F1) / (X1 - X2) + (P1 + P2)) / &
    &      ((X1 - X2) * (X1 - X2))
  FC2(X1,F1,P1,X2,F2,P2) = &
    &      (ZTHREE * (X1 + X2) * (F1 - F2) / (X1 - X2) - &
    &       P1 * (X1 + ZTWO * X2) - P2 * (X2 + ZTWO * X1)) / &
    &      ((X1 - X2) * (X1 - X2))
  FC1(X1,F1,P1,X2,F2,P2) = &
    &      (zsix * X1 * X2 * (F2 - F1) / (X1 - X2) + &
    &       X2 * P1 * (2 * X1 + X2) + X1 * P2 * (X1 + ZTWO * X2)) / &
    &      ((X1 - X2) * (X1 - X2))
  FC0(X1,F1,P1,X2,F2,P2) = &
    &      (F1 * X2**2 + F2 * X1**2 - X1 * X2 * (X2 * P1 + X1 * P2) + &
    &       ZTWO * X1 * X2 * (F1 * X2 - F2 * X1) / (X1 - X2)) / &
    &      ((X1 - X2) * (X1 - X2))
  ! ----------------------------------------------------------------------
  ! -- FCDCD0 GIVES THE VALUE OF THE FUNCTION AT POINT PX               --
  ! -- FCDCD0(......,PX) = F(PX)                                        --
  ! ----------------------------------------------------------------------
  FCDCD0(X1,F1,P1,X2,F2,P2,PX) = &
    &              FC0(X1,F1,P1,X2,F2,P2) + &
    &              PX * (FC1(X1,F1,P1,X2,F2,P2) + &
    &                    PX * (FC2(X1,F1,P1,X2,F2,P2) + &
    &                          PX * FC3(X1,F1,P1,X2,F2,P2)))
  ! ----------------------------------------------------------------------
  ! -- FCDCD1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX:    --
  ! -- FCDCD1(......,PX) = DF/DX (PX)                                   --
  ! ----------------------------------------------------------------------
  FCDCD1(X1,F1,P1,X2,F2,P2,PX) = &
    &              FC1(X1,F1,P1,X2,F2,P2) + &
    &              PX * (ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
    &                    ZTHREE * PX * FC3(X1,F1,P1,X2,F2,P2))
  ! ----------------------------------------------------------------------
  ! -- FCDCD2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX:   --
  ! -- FCDCD2(......,PX) = D2F/DX2 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCDCD2(X1,F1,P1,X2,F2,P2,PX) = &
    &             ZTWO * FC2(X1,F1,P1,X2,F2,P2) + &
    &             zsix * FC3(X1,F1,P1,X2,F2,P2) * PX
  ! ----------------------------------------------------------------------
  ! -- FCDCD3 GIVES THE VALUE OF THE THIRD DERIVATIVE OF F(X) AT PX:    --
  ! -- FCDCD3(......,PX) = D3F/DX3 (PX)                                 --
  ! ----------------------------------------------------------------------
  FCDCD3(X1,F1,P1,X2,F2,P2,PX) = &
    &                      zsix * FC3(X1,F1,P1,X2,F2,P2)
  ! ----------------------------------------------------------------------
  ! -- FCDCDM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  FCDCDM1(X1,F1,P1,X2,F2,P2,PX) = &
    &  (PX-X1)*(FC0(X1,F1,P1,X2,F2,P2) + &
    &  0.5_RKIND*(PX+X1)* FC1(X1,F1,P1,X2,F2,P2)+ &
    &  FC2(X1,F1,P1,X2,F2,P2)/3._RKIND*(PX*(PX+X1)+X1*X1) + &
    &  0.25_RKIND*(PX+X1)*(PX*PX+X1*X1)*FC3(X1,F1,P1,X2,F2,P2))
  !
  !.......................................................................
  !*COMDECK QUAQQQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
  ! ----------------------------------------------------------------------
  !
  FB2(A1,A2,A3,B1,B2,B3) = &
    &               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
  FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))- &
    &         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
  FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1 &
    &         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
  ! ----------------------------------------------------------------------
  ! -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
  ! -- FQQQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * (FB1(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * FB2(A1,A2,A3,B1,B2,B3))
  ! ----------------------------------------------------------------------
  ! -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQQQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) + &
    &     ZTWO * PX * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQQQ2(A1,A2,A3,B1,B2,B3) = ZTWO * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM B1 TO PX:
  !
  FQQQM1(A1,A2,A3,B1,B2,B3,PX) = &
    & (PX-B1)*(FB0(A1,A2,A3,B1,B2,B3) + &
    &  0.5_RKIND*(PX+B1)*FB1(A1,A2,A3,B1,B2,B3) + &
    &  FB2(A1,A2,A3,B1,B2,B3)/3._RKIND*(PX*(PX+B1)+B1*B1))
  !.......................................................................
  !*COMDECK QUAQDQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
  ! -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
  ! ----------------------------------------------------------------------
  !
  FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
  FD1(X1,F1,P1,X2,F2) = P1 - ZTWO*X1*FD2(X1,F1,P1,X2,F2)
  FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) + &
    &                                     FD1(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
  ! -- FQDQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) + &
    &                              PX * (FD1(X1,F1,P1,X2,F2) + &
    &                                    PX * FD2(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQDQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) + &
    &                              ZTWO* PX * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQDQ2(X1,F1,P1,X2,F2) = ZTWO * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  !
  FQDQM1(X1,F1,P1,X2,F2,PX) = &
    & (PX-X1)*(FD0(X1,F1,P1,X2,F2) + &
    &  0.5_RKIND*(PX+X1)*FD1(X1,F1,P1,X2,F2) + &
    &  FD2(X1,F1,P1,X2,F2)/3._RKIND*(PX*(PX+X1)+X1*X1))
  !-----------------------------------------------------------------------
  !.......................................................................
  !     LINEAR
  !
  FLINEAR(X1,F1,X2,F2,PX) = F1 + (PX-X1)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  FLINEARM1(X1,F1,X2,F2,PX) = (PX-X1)*(F1+0.5_RKIND*(PX-X1)*(F2-F1)/(X2-X1))
  FLINXYP(X1,F1,P1,PX) = P1*(PX-X1) + F1
  FLINXYPM1(X1,F1,P1,PX) = (PX-X1)*(F1+0.5_RKIND*P1*(PX-X1))

  !-----------------------------------------------------------------------
  ! 0. DEFAULTS
  !
  ICONTDER = 1
  IF (KOPTXPOL .LT. 0) ICONTDER = 0
  IOPTXPOL=ABS(KOPTXPOL)
  !
  ! COMPUTE INT. UP TO EACH INPUT INTERVAL
  IF (KOPTDER .GE. 3) THEN
    ZYINT_XIN(1) = 0._RKIND
    DO I=2,KNIN
      H = PXIN(I) - PXIN(I-1)
      ZYINT_XIN(I) = ZYINT_XIN(I-1) + H/2._RKIND*(PYIN(I) + PYIN(I-1) &
        & - H*H/12._RKIND*(PYINPP(I)+PYINPP(I-1)))
    END DO
  END IF
  !
  ! LOOP OVER XOUT POINTS WHICH CAN BE IN RANDOM ORDER
  DO 100 J=1,KNOUT
    IF ((PXOUT(J) .LT. PXIN(1)) .OR. (PXOUT(J) .GT. PXIN(KNIN))) GO TO 200
    !
    ! 1.1 POINTS INSIDE INTERVAL [XIN(1),XIN(KNIN)]
    ! FIND PXIN INTERVAL BY BI-SECTION
    KLO=1
    KHI=KNIN
10  CONTINUE
    IF (KHI-KLO.GT.1) THEN
      K=(KHI+KLO)/2
      IF(PXIN(K) .GT. PXOUT(J))THEN
        KHI=K
      ELSE
        KLO=K
      ENDIF
      GOTO 10
    ENDIF
    H=PXIN(KHI)-PXIN(KLO)
    H2 = H * H
    SELECT CASE (KOPTDER)
    CASE (0)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
    CASE (1)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
    CASE (2)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=(PXOUT(J)-PXIN(KLO))/H
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
      PYPP(J)=A*PYINPP(KLO)+B*PYINPP(KHI)
    CASE (3)
      A=(PXIN(KHI)-PXOUT(J))/H
      B=ZONE - A
      A2 = A*A
      B2 = B*B
      PY(J)=A*PYIN(KLO)+B*PYIN(KHI)+ ((A**3-A)*PYINPP(KLO)+(B**3-B)*PYINPP(KHI))*H2/6._RKIND
      PYP(J)=(PYIN(KHI)-PYIN(KLO))/H - &
        & ((3._RKIND*A*A-1._RKIND)*PYINPP(KLO)-(3._RKIND*B*B-1._RKIND)*PYINPP(KHI) )*H/6._RKIND
      PYPP(J)=A*PYINPP(KLO)+B*PYINPP(KHI)
      PYINT(J) = ZYINT_XIN(KLO) + H/2._RKIND*(PYIN(KLO)*(1._RKIND-A2) + PYIN(KHI)*B2 &
        &  - (PYINPP(KLO)*(1._RKIND-A2)**2 + PYINPP(KHI)*B2*(2._RKIND-B2)) * H2/12._RKIND)
    END SELECT
    !
    GO TO 100
    !
    !     2 POINT OUTSIDE INTERVAL
    !
200 CONTINUE
    !
    !     2.1 IF KOPTXPOL=0, PRINT WARNING AND RETURN OR STOP
    !
    IF (IOPTXPOL .EQ. 0) THEN
      PRINT *,' POINT PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL [',PXIN(1),',',PXIN(KNIN),']'
      RETURN
      !        STOP 'IOPTXPOL=0'
    ENDIF
    !
    !     2.2 COMPUTE VALUES FOR POINTS ON THE LEFT OF PXIN(1)
    !           EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
    !
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      H = PXIN(2) - PXIN(1)
      ! YP(PXIN(1)) FROM SPLINE
      ! ZY1 REFERS TO FIRST KNOWN POINTS AND ZY2 TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(1) AND PXIN(2) OR PXIN(1)-ALFA*H AND PXIN(1)
      ZY1P = (PYIN(2)-PYIN(1))/H - (2._RKIND*PYINPP(1)+PYINPP(2))*H/6._RKIND
      ZXLFTDEL = PXIN(1) - ALFA*H
      !
      !   2.2.1 SPECIAL PART [XIN(1)-ALFA*H,XIN(1)] IF IOPTXPOL>20
      IF ((PXOUT(J) .GE. ZXLFTDEL) .AND. (IOPTXPOL.GE.21)) THEN
        IF (IOPTXPOL .EQ. 21) THEN
          ! QUADRATIC
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +21
            PY(J) = FQDQ0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2))
            ! INTEGRATES DIRECTLY FROM PXIN(1)
            IF (KOPTDER .GE. 3) PYINT(J) = FQDQM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),PXOUT(J))
          ELSE
            ! KOPTXPOL = -21
            PY(J) = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQQQ2(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3))
            IF (KOPTDER .GE. 3) PYINT(J) = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          END IF
        ELSEIF (IOPTXPOL .GE. 31) THEN
          ! CUBIC PART OF 31, 32
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +31 OR +32
            ZY2P = (PYIN(2)-PYIN(1))/H + (PYINPP(1)+2._RKIND*PYINPP(2))*H/6._RKIND
            PY(J) = FCDCD0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCDCD1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCDCD2(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCDCDM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
          ELSE
            ! KOPTXPOL = -31 OR -32
            PY(J) = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
          END IF
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          PRINT *,'KOPTXPOL 1'
          RETURN ! AVOID STOPS
        END IF
        !
      ELSE
        !
        ! 2.2.2 EXTRAPOLATION FAR LEFT: X<XIN(1)-ALFA*H OR X<XIN(1) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21) .OR. (IOPTXPOL .EQ. 31)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = ZY1P
            ZYINT = 0._RKIND
          CASE(-1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
            ZYINT = 0._RKIND
          CASE (21)
            ZX1 = ZXLFTDEL
            ZY1 = FQDQ0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
            ZYINT = FQDQM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
            ZY1P = FQDQ1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZXLFTDEL)
          CASE (-21)
            ZX1 = ZXLFTDEL
            ZY1 = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZYINT = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
          CASE (31)
            ZY2P = (PYIN(2)-PYIN(1))/H + (PYINPP(1)+2._RKIND*PYINPP(2))*H/6._RKIND
            ZYINT = FCDCDM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
            ZX1 = ZXLFTDEL
            ZY1 = FCDCD0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
            ZY1P = FCDCD1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
          CASE (-31)
            ZX1 = ZXLFTDEL
            ZY1 = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZY1P = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZYINT = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
          END SELECT
          PY(J) = FLINXYP(ZX1,ZY1,ZY1P,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = ZY1P
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FLINXYPM1(ZX1,ZY1,ZY1P,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(1)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PY(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          ! INTEGRATES FROM PXIN(1) DIRECTLY
          IF (KOPTDER .GE. 3) PYINT(J) = (PXOUT(J)-PXIN(1))*ZYIN_XPOL
          !
        ELSE IF ((IOPTXPOL .EQ. 2) .OR. (IOPTXPOL .EQ. 32)) THEN
          ! QUADRATIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (2)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZX2 = PXIN(2)
            ZY2 = PYIN(2)
            ZYINT = 0._RKIND
          CASE (-2)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXIN(1))
            ZX2 = PXIN(2)
            ZY2 = PYIN(2)
            ZYINT = 0._RKIND
          CASE (32)
            ZY2P = (PYIN(2)-PYIN(1))/H + (PYINPP(1)+2._RKIND*PYINPP(2))*H/6._RKIND
            ZYINT = FCDCDM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
            ZX1 = ZXLFTDEL
            ZY1 = FCDCD0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
            ZY1P = FCDCD1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,ZXLFTDEL)
            ZX2 = PXIN(1)
            ZY2 = PYIN(1)
          CASE (-32)
            ZX1 = ZXLFTDEL
            ZY1 = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZY1P = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
            ZX2 = PXIN(1)
            ZY2 = PYIN(1)
            ZYINT = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),ZXLFTDEL)
          END SELECT
          PY(J) = FQDQ0(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(ZX1,ZY1,ZY1P,ZX2,ZY2)
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 3) THEN
          ! CUBIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (3)
            ZY2P = (PYIN(2)-PYIN(1))/H + (PYINPP(1)+2._RKIND*PYINPP(2))*H/6._RKIND
            PY(J) = FCDCD0(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCDCD1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCDCD2(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCDCDM1(PXIN(1),PYIN(1),ZY1P,PXIN(2),PYIN(2),ZY2P,PXOUT(J))
          CASE (-3)
            PY(J) = FCCCC0(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = FCCCCM1(PYIN(1),PYIN(2),PYIN(3),PYIN(4),PXIN(1),PXIN(2),PXIN(3),PXIN(4),PXOUT(J))
          END SELECT
        ELSE
          PRINT *,'OPTION  IOPTXPOL= ',IOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 2'
        END IF
      END IF
      !
    ELSE ! PXOUT(J) .GT. PXIN(KNIN)
      !
      !     2.3 COMPUTE VALUES FOR POINTS ON THE RIGHT OF PXIN(1)
      !         EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
      !
      H = PXIN(KNIN) - PXIN(KNIN-1)
      ! ZYN REFERS TO FIRST KNOWN POINTS AND ZYNN TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(KNIN) AND PXIN(KNIN-1) OR PXIN(KNIN)+ALFA*H AND PXIN(KNIN)
      ! YP(PXIN(KNIN)) FROM SPLINE
      ZYNP = (PYIN(KNIN)-PYIN(KNIN-1))/H + (PYINPP(KNIN-1)+2._RKIND*PYINPP(KNIN))*H/6._RKIND
      ZXRGTDEL = PXIN(KNIN) + ALFA*H
      !
      !   2.3.1 SPECIAL PART ]XIN(KNIN),XIN(KNIN)+ALFA*H] IF IOPTXPOL>20
      IF ((PXOUT(J) .LE. ZXRGTDEL) .AND. (IOPTXPOL.GE.21)) THEN
        ZYINT = ZYINT_XIN(KNIN)
        IF (IOPTXPOL .EQ. 21) THEN
          ! QUADRATIC
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +21
            PY(J) = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1))
            ! INTEGRATES FROM PXIN(KNIN) SO ADD INTEGRAL UP TO PXIN(KNIN)
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          ELSE
            ! KOPTXPOL = -21
            PY(J) = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = &
              & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = &
              & FQQQ2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          END IF
        ELSEIF (IOPTXPOL .GE. 31) THEN
          ! CUBIC PART OF 31, 32
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +31 OR +32
            ZYNNP = (PYIN(KNIN)-PYIN(KNIN-1))/H - (2._RKIND*PYINPP(KNIN-1)+PYINPP(KNIN))*H/6._RKIND
            PY(J) = FCDCD0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = &
              & FCDCD1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = &
              & FCDCD2(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + &
              & FCDCDM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
          ELSE
            ! KOPTXPOL = -31 OR -32
            PY(J) = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
          END IF
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 3'
        END IF
        !
      ELSE
        !
        ! 2.3.2 EXTRAPOLATION FAR RIGHT: X>XIN(KNIN)+ALFA*H OR X>XIN(KNIN) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21) .OR. (IOPTXPOL .EQ. 31)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYINT = ZYINT_XIN(KNIN)
          CASE(-1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FLINEARP(PXIN(KNIN),PYIN(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
            ZYINT = ZYINT_XIN(KNIN)
          CASE (21)
            ZYINT = ZYINT_XIN(KNIN) + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZYNP = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
          CASE (-21)
            ZXN = ZXRGTDEL
            ZYN = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYINT = ZYINT_XIN(KNIN) + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
          CASE (31)
            ZXN = ZXRGTDEL
            ZYNNP = (PYIN(KNIN)-PYIN(KNIN-1))/H - (2._RKIND*PYINPP(KNIN-1)+PYINPP(KNIN))*H/6._RKIND
            ZYINT = ZYINT_XIN(KNIN) + &
              & FCDCDM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
            ZYN = FCDCD0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
            ZYNP = FCDCD1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
          CASE (-31)
            ZXN = ZXRGTDEL
            ZYN = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYNP = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYINT = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
          END SELECT
          PY(J) = FLINXYP(ZXN,ZYN,ZYNP,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = ZYNP
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FLINXYPM1(ZXN,ZYN,ZYNP,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(KNIN)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PY(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + (PXOUT(J)-PXIN(KNIN))*ZYIN_XPOL
          !
        ELSE IF ((IOPTXPOL .EQ. 2) .OR. (IOPTXPOL .EQ. 32)) THEN
          ! QUADRATIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (2)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZXNN = PXIN(KNIN-1)
            ZYNN = PYIN(KNIN-1)
            ZYINT = ZYINT_XIN(KNIN)
          CASE (-2)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN))
            ZXNN = PXIN(KNIN-1)
            ZYNN = PYIN(KNIN-1)
            ZYINT = ZYINT_XIN(KNIN)
          CASE (32)
            ZXN = ZXRGTDEL
            ZYNNP = (PYIN(KNIN)-PYIN(KNIN-1))/H - (2._RKIND*PYINPP(KNIN-1)+PYINPP(KNIN))*H/6._RKIND
            ZYINT = ZYINT_XIN(KNIN) + FCDCDM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
            ZXNN = PXIN(KNIN)
            ZYNN = PYIN(KNIN)
            ZYN = FCDCD0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
            ZYNP = FCDCD1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,ZXRGTDEL)
          CASE (-32)
            ZXN = ZXRGTDEL
            ZYN = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZYNP = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
            ZXNN = PXIN(KNIN)
            ZYNN = PYIN(KNIN)
            ZYINT = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),ZXRGTDEL)
          END SELECT
          PY(J) = FQDQ0(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 1) PYP(J) = FQDQ1(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 2) PYPP(J) = FQDQ2(ZXN,ZYN,ZYNP,ZXNN,ZYNN)
          IF (KOPTDER .GE. 3) PYINT(J) = ZYINT + FQDQM1(ZXN,ZYN,ZYNP,ZXNN,ZYNN,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 3) THEN
          ! CUBIC EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (3)
            ZYNNP = (PYIN(KNIN)-PYIN(KNIN-1))/H - (2._RKIND*PYINPP(KNIN-1)+PYINPP(KNIN))*H/6._RKIND
            PY(J) = FCDCD0(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCDCD1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCDCD2(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + &
              & FCDCDM1(PXIN(KNIN),PYIN(KNIN),ZYNP,PXIN(KNIN-1),PYIN(KNIN-1),ZYNNP,PXOUT(J))
          CASE (-3)
            PY(J) = FCCCC0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 1) PYP(J) = FCCCC1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 2) PYPP(J) = FCCCC2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
            IF (KOPTDER .GE. 3) PYINT(J) = ZYINT_XIN(KNIN) + FCCCCM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PYIN(KNIN-3), &
              & PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN-3),PXOUT(J))
          END SELECT
        ELSE
          PRINT *,'OPTION  KOPTXPOL= ',KOPTXPOL,' NOT YET DEFINED'
          STOP 'KOPTXPOL 4'
        END IF
      END IF
    END IF
    !
100 END DO
  RETURN
END SUBROUTINE SPLIBNDA
!.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!
SUBROUTINE INTLINEAR(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KEXTRAPO)
  !
  !   COMPUTE LINEAR INTERPOLATION OF (PXIN,PYIN) ON (PXOUT,PYOUT).
  !
  !   KOPTDER = 0: COMPUTE ONLY PYOUT (PYOUTP, PYOUTPP NOT USED)
  !   KOPTDER = 1: COMPUTE ALSO 1ST DER. IN PYOUTP (PYOUTPP NOT USED)
  !   KOPTDER = 2: AS 1 AND 2ND DER. IN PYOUTPP
  !   KOPTDER = 3: AS 2 AND INTEGRAL FROM (PXIN(1) IN PYOUTPP
  !
  !   KEXTRAPO = 0: PRINT MESSAGE AND RETURN IF NEED TO EXTRAPOLATE
  !   KEXTRAPO = 1: NORMAL EXTRAPOLATION
  !   KEXTRAPO = 10: Y=Y_EDGE FOR EXTRAPOLATION => CONSTANT EXTRAPOLATION
  !   KEXTRAPO = -10: Y=0 FOR EXTRAPOLATION
  !
  USE PREC_RKIND
  IMPLICIT NONE
  INTEGER :: KNIN, KNOUT, KOPTDER, KEXTRAPO
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN)
  REAL(RKIND) :: PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
  REAL(RKIND) :: ZYINT_XIN(KNIN)
  ! FOR FUNCTIONS
  REAL(RKIND) :: X1, X2, F1, F2, PX, FLINEAR, FLINEARP, FLINEARM1
  INTEGER :: K, KLO, KHI, I, J
  !
  FLINEAR(X1,F1,X2,F2,PX) = F1 + (PX-X1)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  FLINEARM1(X1,F1,X2,F2,PX) = (PX-X1)*(F1+0.5_RKIND*(PX-X1)*(F2-F1)/(X2-X1))
  !.......................................................................
  !
  IF (KOPTDER .GE. 3) THEN
    ZYINT_XIN(1) = 0._RKIND
    DO I=2,KNIN
      ZYINT_XIN(I) = ZYINT_XIN(I-1) + 0.5_RKIND*(PXIN(I)-PXIN(I-1))*(PYIN(I)+PYIN(I-1))
    END DO
  END IF
  !
  DO J=1,KNOUT
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      SELECT CASE (KEXTRAPO)
      CASE (1)
        PYOUT(J) = FLINEAR(PXIN(1),PYIN(1),PXIN(2),PYIN(2),PXOUT(J))
        IF (KOPTDER .GE. 1) PYOUTP(J) = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = FLINEARM1(PXIN(1),PYIN(1),PXIN(2),PYIN(2),PXOUT(J))
      CASE (0)
        PRINT *,' WARNING POINT PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL ON LEFT'
        RETURN
      CASE (10)
        PYOUT(J) = PYIN(1)
        IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = PYIN(1)*(PXOUT(J)-PXIN(1))
      CASE (-10)
        PYOUT(J) = 0._RKIND
        IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = 0._RKIND
      CASE DEFAULT
        PYOUT(J) = FLINEAR(PXIN(1),PYIN(1),PXIN(2),PYIN(2),PXOUT(J))
        IF (KOPTDER .GE. 1) PYOUTP(J) = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = FLINEARM1(PXIN(1),PYIN(1),PXIN(2),PYIN(2),PXOUT(J))
      END SELECT
      !
    ELSE IF (PXOUT(J) .GT. PXIN(KNIN)) THEN
      ! POINT ON RIGHT OF INTERVAL
      SELECT CASE (KEXTRAPO)
      CASE (1)
        PYOUT(J) = FLINEAR(PXIN(KNIN-1),PYIN(KNIN-1),PXIN(KNIN),PYIN(KNIN),PXOUT(J))
        IF (KOPTDER .GE. 1) PYOUTP(J) = FLINEARP(PXIN(KNIN-1),PYIN(KNIN-1),PXIN(KNIN),PYIN(KNIN))
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KNIN) + &
          & FLINEARM1(PXIN(KNIN),PYIN(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
      CASE (0)
        PRINT *,' WARNING POINT PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL ON RIGHT'
        RETURN
      CASE(10)
        PYOUT(J) = PYIN(KNIN)
        IF (KOPTDER .GE. 1) PYOUTP(J) = 0.0_RKIND
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0.0_RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KNIN) + PYIN(KNIN)*(PXOUT(J)-PXIN(KNIN))
      CASE(-10)
        PYOUT(J) = 0.0_RKIND
        IF (KOPTDER .GE. 1) PYOUTP(J) = 0.0_RKIND
        IF (KOPTDER .GE. 2) PYOUTPP(J) = 0.0_RKIND
        IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KNIN)
      END SELECT
      !
    ELSE
      ! FIND PXIN INTERVAL BY BI-SECTION
      KLO=1
      KHI=KNIN
10    CONTINUE
      IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(PXIN(K) .GT. PXOUT(J))THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 10
      ENDIF
      PYOUT(J) = FLINEAR(PXIN(KLO),PYIN(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      IF (KOPTDER .GE. 1) PYOUTP(J) = FLINEARP(PXIN(KLO),PYIN(KLO),PXIN(KHI),PYIN(KHI))
      IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
      IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KLO) + &
        & FLINEARM1(PXIN(KLO),PYIN(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      !
    END IF
  END DO
  !
  RETURN
END SUBROUTINE INTLINEAR
!......................................................................
!     
SUBROUTINE INTQUADRATIC(PXIN,PYIN,KNIN,PXOUT,PYOUT,PYOUTP,PYOUTPP,PYOUTINT,KNOUT,KOPTDER,KOPTXPOL,NBC)
  !
  !   COMPUTE LINEAR INTERPOLATION OF (PXIN,PYIN) ON (PXOUT,PYOUT).
  !
  !   KOPTDER = 0: COMPUTE ONLY PYOUT (PYOUTP, PYOUTPP NOT USED)
  !   KOPTDER = 1: COMPUTE ALSO 1ST DER. IN PYOUTP (PYOUTPP NOT USED)
  !   KOPTDER = 2: AS 1 AND 2ND DER. IN PYOUTPP
  !   KOPTDER = 3: AS 2 AND INTEGRAL FROM (PXIN(1) IN PYOUTPP
  !
  !   KOPTXPOL = 0: SEND MESSAGE IF NEED TO EXTRAPOLATE
  !   KOPTXPOL = 1: LINEAR EXTRAPOLATION WITH CONTINUOUS DERIVATIVE
  !   KOPTXPOL = -1: LINEAR EXTRAPOLATION WITH LAST THREE Y VALUES
  !   KOPTXPOL = 2: QUADRATIC EXTRAPOLATION WITH CONTINUOUS DERIVATIVE
  !   KOPTXPOL = -2: QUADRATIC EXTRAPOLATION WITH LAST THREE Y VALUES
  !   KOPTXPOL = 21: QUADRATIC EXTRAPOLATION WITHIN ONE DELTA_X, THEN LINEAR
  !   KOPTXPOL = -21: SAME AS 21 BUT USES VALUES INSTEAD OF CLOSEST DERIVATIVE FOR EXTRAPOLATION
  !   KOPTXPOL = 10: Y=Y_EDGE FOR EXTRAPOLATION => CONSTANT EXTRAPOLATION
  !   KOPTXPOL = -10: Y=0 FOR EXTRAPOLATION
  !
  ! SOLVES FOR QUADRATICS WITH CONTINUOUS DERIVATIVE AND 
  !   NBC = 0:  2ND DERIVATIVE =0 AT LEFT (DEFAULT)
  !   NBC = 1:  2ND DERIVATIVE =0 AT RIGHT
  !
  USE PREC_RKIND
  IMPLICIT NONE
  REAL(RKIND) :: ZSIX, ZTHREE, ZTWO, ZONE
  PARAMETER(ZSIX=6._RKIND, ZTHREE=3._RKIND, ZTWO=2._RKIND, ZONE=1._RKIND)
  INTEGER KNIN,KNOUT,KOPTDER, KOPTXPOL
  REAL(RKIND) :: PXIN(KNIN), PYIN(KNIN)
  REAL(RKIND) :: PXOUT(KNOUT), PYOUT(KNOUT), PYOUTP(KNOUT), PYOUTPP(KNOUT), PYOUTINT(KNOUT)
  !INTEGER, OPTIONAL ::  NBC
  INTEGER ::  NBC
  !
  REAL(RKIND) :: BCOEF(KNIN-1), ZYINP(KNIN), ZYINT_XIN(KNIN), ZYIN_XPOL, ZYINT, &
    & ZXLFTDEL, ZXRGTDEL, ALFA, ZX1, ZXN, ZY1, ZYN, ZY1P, ZX2, ZXNN, ZYNP, ZY2, ZYNN
  INTEGER I,J,IOPTXPOL, INBC, ICONTDER, K, KLO, KHI
  !
  !
  ! VARIABLES RELATED TO FUNCTIONS:
  REAL(RKIND) :: FQQQ0, FQQQ1, FQQQ2, &
    &  FLINEAR, FLINEARP,  &
    &  FQQQM1, FQDQM1, FLINEARM1, FLINXYP, FLINXYPM1, FPARABOLP
  REAL(RKIND) :: A1, A2, A3, A4, B1, B2, B3, B4, PX
  REAL(RKIND) :: FB0, FB1, FB2, FA0, FA1, FD2, FD1, FD0, FQDQ0, FQDQ1, FQDQ2
  REAL(RKIND) :: X1, F1, P1, X2, F2
  REAL(RKIND) :: X3, F3
  !
  !.......................................................................
  !*COMDECK QUAQQQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE SIX PARAMETERS A1,A2,A3,B1,B2,B3 ARE DEFINED AS FOLLOWS:     --
  ! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3                             --
  ! ----------------------------------------------------------------------
  !
  FB2(A1,A2,A3,B1,B2,B3) = &
    &               ((A1-A2)/(B1-B2)-(A1-A3)/(B1-B3))/(B2-B3)
  FB1(A1,A2,A3,B1,B2,B3) = ((A1-A2)/(B1-B2))- &
    &         FB2(A1,A2,A3,B1,B2,B3)*(B1+B2)
  FB0(A1,A2,A3,B1,B2,B3) = A1-FB1(A1,A2,A3,B1,B2,B3)*B1 &
    &         -FB2(A1,A2,A3,B1,B2,B3)*B1*B1
  ! ----------------------------------------------------------------------
  ! -- FQQQ0 GIVES THE VALUE OF THE FUNCTION AT THE POINT PX            --
  ! -- FQQQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQQQ0(A1,A2,A3,B1,B2,B3,PX) = FB0(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * (FB1(A1,A2,A3,B1,B2,B3) + &
    &                                 PX * FB2(A1,A2,A3,B1,B2,B3))
  ! ----------------------------------------------------------------------
  ! -- FQQQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQQQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQQQ1(A1,A2,A3,B1,B2,B3,PX) = FB1(A1,A2,A3,B1,B2,B3) + &
    &     ZTWO * PX * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQQQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQQQ2(A1,A2,A3,B1,B2,B3) = ZTWO * FB2(A1,A2,A3,B1,B2,B3)
  ! ----------------------------------------------------------------------
  ! -- FQQQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM B1 TO PX:
  !
  FQQQM1(A1,A2,A3,B1,B2,B3,PX) = &
    & (PX-B1)*(FB0(A1,A2,A3,B1,B2,B3) + &
    &  0.5_RKIND*(PX+B1)*FB1(A1,A2,A3,B1,B2,B3) + &
    &  FB2(A1,A2,A3,B1,B2,B3)/3._RKIND*(PX*(PX+B1)+B1*B1))
  !.......................................................................
  !*COMDECK QUAQDQ
  ! ----------------------------------------------------------------------
  ! --     STATEMENT FUNCTION FOR QUADRATIC INTERPOLATION               --
  ! --                         19.01.87            AR        CRPP       --
  ! --                                                                  --
  ! -- QUADRATIC INTERPOLATION OF A FUNCTION F(X)                       --
  ! -- THE FIVE PARAMETERS X1,F1,P1,X2,F2    ARE DEFINED AS FOLLOWS:    --
  ! -- F(X1) = F1 , DF/DX(X1) = P1 , F(X2) = F2                         --
  ! ----------------------------------------------------------------------
  !
  FD2(X1,F1,P1,X2,F2) = ((F2-F1)/(X2-X1) - P1) / (X2-X1)
  FD1(X1,F1,P1,X2,F2) = P1 - ZTWO*X1*FD2(X1,F1,P1,X2,F2)
  FD0(X1,F1,P1,X2,F2) = F1 - X1*(X1*FD2(X1,F1,P1,X2,F2) + &
    &                                     FD1(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ0 GIVES THE VALUE OF THE FUNCTION AT POINT PX                --
  ! -- FQDQ0(......,PX) = F(PX)                                         --
  ! ----------------------------------------------------------------------
  FQDQ0(X1,F1,P1,X2,F2,PX) = FD0(X1,F1,P1,X2,F2) + &
    &                              PX * (FD1(X1,F1,P1,X2,F2) + &
    &                                    PX * FD2(X1,F1,P1,X2,F2))
  ! ----------------------------------------------------------------------
  ! -- FQDQ1 GIVES THE VALUE OF THE FIRST DERIVATIVE OF F(X) AT PX      --
  ! -- FQDQ1(......,PX) = DF/DX (PX)                                    --
  ! ----------------------------------------------------------------------
  FQDQ1(X1,F1,P1,X2,F2,PX) = FD1(X1,F1,P1,X2,F2) + &
    &                              ZTWO* PX * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQ2 GIVES THE VALUE OF THE SECOND DERIVATIVE OF F(X) AT PX     --
  ! -- FQDQ2(......,PX) = D2F/DX2 (PX)                                  --
  ! ----------------------------------------------------------------------
  FQDQ2(X1,F1,P1,X2,F2) = ZTWO * FD2(X1,F1,P1,X2,F2)
  ! ----------------------------------------------------------------------
  ! -- FQDQM1 GIVES THE VALUE OF THE INTEGRAL OF F(X) FROM X1 TO PX:
  !
  FQDQM1(X1,F1,P1,X2,F2,PX) = &
    & (PX-X1)*(FD0(X1,F1,P1,X2,F2) + &
    &  0.5_RKIND*(PX+X1)*FD1(X1,F1,P1,X2,F2) + &
    &  FD2(X1,F1,P1,X2,F2)/3._RKIND*(PX*(PX+X1)+X1*X1))
  !-----------------------------------------------------------------------
  FPARABOLP(X1,X2,X3,F1,F2,F3,PX) = &
    &  ((PX-X1)+(PX-X2))*F3/((X3-X1)*(X3-X2))+ &
    &  ((PX-X1)+(PX-X3))*F2/((X2-X1)*(X2-X3))+ &
    &  ((PX-X2)+(PX-X3))*F1/((X1-X2)*(X1-X3))
  !.......................................................................
  !     LINEAR
  !
  FLINEAR(X1,F1,X2,F2,PX) = F1 + (PX-X1)/(X2-X1) * (F2-F1)
  FLINEARP(X1,F1,X2,F2) = (F2-F1) / (X2-X1)
  FLINEARM1(X1,F1,X2,F2,PX) = (PX-X1)*(F1+0.5_RKIND*(PX-X1)*(F2-F1)/(X2-X1))
  FLINXYP(X1,F1,P1,PX) = P1*(PX-X1) + F1
  FLINXYPM1(X1,F1,P1,PX) = (PX-X1)*(F1+0.5_RKIND*P1*(PX-X1))
  !
  !-----------------------------------------------------------------------
  ! 0. DEFAULTS
  !
  ICONTDER = 1
  IF (KOPTXPOL .LT. 0) ICONTDER = 0
  IOPTXPOL=ABS(KOPTXPOL)
  INBC = 0
  !IF (PRESENT(NBC)) INBC = NBC
  INBC = NBC
  !
  ! DETERMINE DERIVATIVE ZYINP(I) SUCH THAT QUADRATICS ARE THEN DEFINED BY XI, YI, YPI, XI+1, YI+1 WITHIN [XI,XI+1]
  ! VARIOUS OPTIONS ^CHOSEN THROUGH NBC
  !
!!$  ACOEF(1:KNIN-1) = PYIN(1:KNIN-1)
!!$  CCOEF(1:KNIN-1) = PYIN(2:KNIN)
  IF (INBC .EQ. 1) THEN
    ! FIND COEFFICIENTS: Y(X)=AI*(XI+1-X)^2 + BI*(XI+1-X)*(X-XI) + CI*(X-XI)^2 FOR [XI,XI+1], I=1,N-1
    ! AI = YI ; CI=YI+1
    ! BI/HI + BI+1/HI+1 = 2 YI+1 (1/HI+1/HI+1) FOR CONTINUITY OF DERIVATIVES
    ! NEED ONE BOUNDARY CONDITION
    ! d2y/dx2=0 at knin
    BCOEF(KNIN-1) = PYIN(KNIN-1) + PYIN(KNIN)
    DO I=KNIN-2,1,-1
      BCOEF(I) = (PXIN(I+1)-PXIN(I))*(-BCOEF(I+1)/(PXIN(I+2)-PXIN(I+1)) + & 
        & 2._RKIND*PYIN(I+1)*(1._RKIND/(PXIN(I+1)-PXIN(I))+1._RKIND/(PXIN(I+2)-PXIN(I+1))))
    END DO
    ! d2y/dx2=0 at 1
    !equiv    BCOEF(1) = PYIN(1) + PYIN(2)
    !equiv    DO I=2,KNIN-1
    !equiv      BCOEF(I) = (PXIN(I+1)-PXIN(I))*(-BCOEF(I-1)/(PXIN(I)-PXIN(I-1)) + & 
    !equiv        & 2._RKIND*PYIN(I)*(1._RKIND/(PXIN(I)-PXIN(I-1))+1._RKIND/(PXIN(I+1)-PXIN(I))))
    !equiv    END DO
    DO I=1,KNIN-1
      ZYINP(I) = (BCOEF(I)-2._RKIND*PYIN(I))/(PXIN(I+1)-PXIN(I))
    END DO
    ZYINP(KNIN) = (-BCOEF(KNIN-1)+2._RKIND*PYIN(KNIN))/(PXIN(KNIN)-PXIN(KNIN-1))
  ELSE
    ! USE 3 POINTS TO DEFINE QUADRATICS EXCEPT FOR LAST INTERVALS
    DO I=1,KNIN-2
      ZYINP(I) = FPARABOLP(PXIN(I),PXIN(I+1),PXIN(I+2),PYIN(I),PYIN(I+1),PYIN(I+2),PXIN(I))
    END DO
    ZYINP(KNIN-1) = FPARABOLP(PXIN(KNIN-2),PXIN(KNIN-1),PXIN(KNIN),PYIN(KNIN-2),PYIN(KNIN-1),PYIN(KNIN),PXIN(KNIN-1))
    ZYINP(KNIN) = FPARABOLP(PXIN(KNIN-2),PXIN(KNIN-1),PXIN(KNIN),PYIN(KNIN-2),PYIN(KNIN-1),PYIN(KNIN),PXIN(KNIN))
  END IF
  !write(*,'(1p2e14.5)') (pxin(i),zyinp(i),i=1,knin)
  !
  ! COMPUTE INT. UP TO EACH INPUT INTERVAL
  IF (KOPTDER .GE. 3) THEN
    ZYINT_XIN(1) = 0._RKIND
    DO I=1,KNIN-1
      ZYINT_XIN(I+1) = ZYINT_XIN(I) + FQDQM1(PXIN(I),PYIN(I),ZYINP(I),PXIN(I+1),PYIN(I+1),PXIN(I+1))
    END DO
  END IF
  !
  ! LOOP OVER PXOUT POINTS WHICH CAN BE IN RANDOM ORDER
  DO 100 J=1,KNOUT
    IF ((PXOUT(J) .LT. PXIN(1)) .OR. (PXOUT(J) .GT. PXIN(KNIN))) GO TO 200
    !
    ! 1.1 POINTS INSIDE INTERVAL [XIN(1),XIN(KNIN)]
    ! FIND PXIN INTERVAL BY BI-SECTION
    KLO=1
    KHI=KNIN
10  CONTINUE
    IF (KHI-KLO.GT.1) THEN
      K=(KHI+KLO)/2
      IF(PXIN(K) .GT. PXOUT(J))THEN
        KHI=K
      ELSE
        KLO=K
      ENDIF
      GOTO 10
    ENDIF
    SELECT CASE (KOPTDER)
    CASE (0)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    CASE (1)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    CASE (2)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTPP(J)=FQDQ2(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI))
    CASE (3)
      PYOUT(J)=FQDQ0(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTP(J)=FQDQ1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
      PYOUTPP(J)=FQDQ2(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI))
      PYOUTINT(J) = ZYINT_XIN(KLO) + FQDQM1(PXIN(KLO),PYIN(KLO),ZYINP(KLO),PXIN(KHI),PYIN(KHI),PXOUT(J))
    END SELECT
    !
    GO TO 100
    !
    !     2 POINT OUTSIDE INTERVAL
    !
200 CONTINUE
    !
    !     2.1 IF KOPTXPOL=0, PRINT WARNING AND RETURN OR STOP
    !
    IF (IOPTXPOL .EQ. 0) THEN
      PRINT *,' PXOUT(',J,')=',PXOUT(J),' OUTSIDE INTERVAL [',PXIN(1),',',PXIN(KNIN),']'
      RETURN
      !        STOP 'IOPTXPOL=0'
    ENDIF
    !
    !     2.2 COMPUTE VALUES FOR POINTS ON THE LEFT OF PXIN(1)
    !           EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
    !
    IF (PXOUT(J) .LT. PXIN(1)) THEN
      ! ZY1 REFERS TO FIRST KNOWN POINTS AND ZY2 TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(1) AND PXIN(2) OR PXIN(1)-ALFA*H AND PXIN(1)
      ZXLFTDEL = PXIN(1) - ALFA*(PXIN(2) - PXIN(1))
      !
      !   2.2.1 SPECIAL PART [XIN(1)-ALFA*H,XIN(1)] IF IOPTXPOL>20
      IF ((PXOUT(J) .GE. ZXLFTDEL) .AND. (IOPTXPOL.EQ.21)) THEN
        ! QUADRATIC PART
        IF (ICONTDER .EQ. 1) THEN
          PYOUT(J) = FQDQ0(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2))
          ! INTEGRATES DIRECTLY FROM PXIN(1)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = FQDQM1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),PXOUT(J))
        ELSE
          ! KOPTXPOL = -21
          PYOUT(J) = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQQQ2(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3))
          IF (KOPTDER .GE. 3) PYOUTINT(J) = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXOUT(J))
        END IF
        !
      ELSE
        !
        ! 2.2.2 EXTRAPOLATION FAR LEFT: X<XIN(1)-ALFA*H OR X<XIN(1) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = ZYINP(1)
            ZYINT = 0._RKIND
          CASE(-1)
            ZX1 = PXIN(1)
            ZY1 = PYIN(1)
            ZY1P = FLINEARP(PXIN(1),PYIN(1),PXIN(2),PYIN(2))
            ZYINT = 0._RKIND
          CASE (21)
            ZX1 = ZXLFTDEL
            ZY1 = FQDQ0(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
            ZY1P = FQDQ1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
            ZYINT = FQDQM1(PXIN(1),PYIN(1),ZYINP(1),PXIN(2),PYIN(2),ZXLFTDEL)
          CASE (-21)
            ZX1 = ZXLFTDEL
            ZY1 = FQQQ0(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
            ZYINT = FQQQM1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),ZXLFTDEL)
          END SELECT
          PYOUT(J) = FLINXYP(ZX1,ZY1,ZY1P,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = ZY1P
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FLINXYPM1(ZX1,ZY1,ZY1P,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(1)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PYOUT(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          ! INTEGRATES FROM PXIN(1) DIRECTLY
          IF (KOPTDER .GE. 3) PYOUTINT(J) = (PXOUT(J)-PXIN(1))*ZYIN_XPOL
          !
        ELSE
          ! QUADRATIC EXTRAPOLATION
          ZX1 = PXIN(1)
          ZY1 = PYIN(1)
          ZX2 = PXIN(2)
          ZY2 = PYIN(2)
          ZYINT = 0._RKIND
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = 2
            ZY1P = ZYINP(1)
          ELSE
            ! KOPTXPOL = -2
            ZY1P = FQQQ1(PYIN(1),PYIN(2),PYIN(3),PXIN(1),PXIN(2),PXIN(3),PXIN(1))
          ENDIF
          PYOUT(J) = FQDQ0(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(ZX1,ZY1,ZY1P,ZX2,ZY2)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(ZX1,ZY1,ZY1P,ZX2,ZY2,PXOUT(J))
        END IF
      END IF
      !
    ELSE ! PXOUT(J) .GT. PXIN(KNIN)
      !
      !     2.3 COMPUTE VALUES FOR POINTS ON THE RIGHT OF PXIN(1)
      !         EXTRAPOLATION DEPENDS ON VALUE OF KOPTXPOL
      !
      ! ZYN REFERS TO FIRST KNOWN POINTS AND ZYNN TO SECOND KNOWN POINT, 
      ! TYPICALLY AT PXIN(KNIN) AND PXIN(KNIN-1) OR PXIN(KNIN)+ALFA*H AND PXIN(KNIN)
      ! YP(PXIN(KNIN)) FROM SPLINE
      ZXRGTDEL = PXIN(KNIN) + ALFA * (PXIN(KNIN) - PXIN(KNIN-1))
      !
      !   2.3.1 SPECIAL PART ]XIN(KNIN),XIN(KNIN)+ALFA*H] IF IOPTXPOL>20
      IF ((PXOUT(J) .LE. ZXRGTDEL) .AND. (IOPTXPOL.EQ.21)) THEN
        ZYINT = ZYINT_XIN(KNIN)
        ! QUADRATIC
        IF (ICONTDER .EQ. 1) THEN
          ! KOPTXPOL = +21
          PYOUT(J) = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
          ! INTEGRATES FROM PXIN(KNIN) SO ADD INTEGRAL UP TO PXIN(KNIN)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),PXOUT(J))
        ELSE
          ! KOPTXPOL = -21
          PYOUT(J) = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = &
               & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = &
               & FQQQ2(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2))
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + &
               & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXOUT(J))
        END IF
        !
      ELSE
        !
        ! 2.3.2 EXTRAPOLATION FAR RIGHT: X>XIN(KNIN)+ALFA*H OR X>XIN(KNIN) IF NO ALFA PART CONSIDERED
        !
        IF ((IOPTXPOL .EQ. 1) .OR. (IOPTXPOL .EQ. 21)) THEN
          ! LINEAR EXTRAPOLATION
          SELECT CASE (KOPTXPOL)
          CASE (1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = ZYINP(KNIN)
            ZYINT = ZYINT_XIN(KNIN)
          CASE(-1)
            ZXN = PXIN(KNIN)
            ZYN = PYIN(KNIN)
            ZYNP = FLINEARP(PXIN(KNIN),PYIN(KNIN),PXIN(KNIN-1),PYIN(KNIN-1))
            ZYINT = ZYINT_XIN(KNIN)
          CASE (21)
            ZYINT = ZYINT_XIN(KNIN) + FQDQM1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQDQ0(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
            ZYNP = FQDQ1(PXIN(KNIN),PYIN(KNIN),ZYINP(KNIN),PXIN(KNIN-1),PYIN(KNIN-1),ZXRGTDEL)
          CASE (-21)
            ZYINT = ZYINT_XIN(KNIN) + &
              & FQQQM1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZXN = ZXRGTDEL
            ZYN = FQQQ0(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
            ZYNP = FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),ZXRGTDEL)
          END SELECT
          PYOUT(J) = FLINXYP(ZXN,ZYN,ZYNP,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = ZYINP(KNIN)
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FLINXYPM1(ZXN,ZYN,ZYNP,PXOUT(J))
        ELSE IF (IOPTXPOL .EQ. 10) THEN
          ! CONSTANT OUTSIDE PXIN
          IF (ICONTDER .EQ. 1) THEN
            ! KOPTXPOL = +10
            ZYIN_XPOL = PYIN(KNIN)
          ELSE
            ! KOPTXPOL = -10
            ZYIN_XPOL = 0._RKIND
          END IF
          PYOUT(J) = ZYIN_XPOL
          IF (KOPTDER .GE. 1) PYOUTP(J) = 0._RKIND
          IF (KOPTDER .GE. 2) PYOUTPP(J) = 0._RKIND
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT_XIN(KNIN) + (PXOUT(J)-PXIN(KNIN))*ZYIN_XPOL
          !
        ELSE
          ! QUADRATIC EXTRAPOLATION
          ZXN = PXIN(KNIN)
          ZYN = PYIN(KNIN)
          ZXNN = PXIN(KNIN-1)
          ZYNN = PYIN(KNIN-1)
          ZYINT = ZYINT_XIN(KNIN)
          IF (KOPTXPOL .EQ. -2) ZYINP(KNIN) = &
               & FQQQ1(PYIN(KNIN),PYIN(KNIN-1),PYIN(KNIN-2),PXIN(KNIN),PXIN(KNIN-1),PXIN(KNIN-2),PXIN(KNIN))
          PYOUT(J) = FQDQ0(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 1) PYOUTP(J) = FQDQ1(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
          IF (KOPTDER .GE. 2) PYOUTPP(J) = FQDQ2(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN)
          IF (KOPTDER .GE. 3) PYOUTINT(J) = ZYINT + FQDQM1(ZXN,ZYN,ZYINP(KNIN),ZXNN,ZYNN,PXOUT(J))
        END IF
      END IF
    END IF
    !
100 END DO
  RETURN 
END SUBROUTINE INTQUADRATIC
      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(RKIND) ::               AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
     &                   JU, K2, KM, KV, NB, NW
      REAL(RKIND) ::     TEMP
!     ..
!     .. Local Arrays ..
      REAL(RKIND) ::     WORK13( LDWORK, NBMAX ), &
     &                   WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL, &
     &                   SSWAP, STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
     &   RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1 .LT. J+KL ) THEN
!
                        CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
     &                        1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ ) &
     &               CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
     &                          AB( KV, JJ+1 ), LDAB-1, &
     &                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 ) &
     &               INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 ) &
     &            CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use SLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
     &                      IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2 .GT. 0 ) THEN
!
!                 Update A12
!
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2 .GT. 0 ) THEN
!
!                    Update A22
!
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J2, &
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3 .GT. 0 ) THEN
!
!                    Update A32
!
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J2, &
     &                           JB, -ONE, WORK31, LDWORK, &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3 .GT. 0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
     &                        WORK13, LDWORK )
!
                  IF( I2 .GT. 0 ) THEN
!
!                    Update A23
!
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J3, &
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
     &                           LDAB-1 )
                  END IF
!
                  IF( I3 .GT. 0 ) THEN
!
!                    Update A33
!
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J3, &
     &                           JB, -ONE, WORK31, LDWORK, WORK13, &
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1 .LT. J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 ) &
     &            CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of DGBTRF
!
      END

      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
     &                   INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993 
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(RKIND) ::     AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by DGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) REAL array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE
      PARAMETER          ( ONE = 1.0E+0_RKIND )
!     ..
!     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
     &   RETURN
!
      KD = KU + KL + 1
      LNOTI = KL.GT.0
!
      IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) &
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
     &                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
!
         DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
     &                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
!
      ELSE
!
!        Solve A'*X = B.
!
         DO 30 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
     &                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
!
!        Solve L'*X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
     &                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J ) &
     &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of DGBTRS
!
      END

      SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(RKIND) ::     AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
     &   RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 ) &
     &         CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
     &                     AB( KV+1, J ), LDAB-1 )
!
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J ) &
     &            CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
     &                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
     &                       LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 ) &
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of SGBTF2
!
      END

      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(RKIND) ::     A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.EQ.0 ) &
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I ) &
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I ) &
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I ) &
     &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of SLASWP
!
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
     &                 N4 )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
     &       ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
     &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
     &             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
     &             ( IC.GE.162 .AND. IC.LE.169 ) ) &
     &            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
     &   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
     &            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ),RKIND )*1.6E0_RKIND )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
     &      'an illegal value' )
!
!     End of XERBLA
!
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
     &   RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
     &       INTA.GE.145 .AND. INTA.LE.153 .OR. &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR. &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     all_sgbtrf_s.blas.f
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      integer function isamax(n,sx,incx)
!
      USE prec_rkind
      implicit none
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(RKIND) :: sx(*),smax
      integer i,incx,ix,n
!
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
      subroutine scopy(n,sx,incx,sy,incy)
!
      USE prec_rkind
      implicit none
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(RKIND) :: sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
      integer nel
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 nel=7
      m = mod(n,nel)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     &                   BETA, C, LDC )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL(RKIND) ::     ALPHA, BETA
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  SGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(RKIND) ::     TEMP
!     .. Parameters ..
      REAL(RKIND) ::     ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND. &
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND. &
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SGEMM .
!
      END
      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     &                   BETA, Y, INCY )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      REAL(RKIND) ::     ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     .. Local Scalars ..
      REAL(RKIND) ::     TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
     &         .NOT.LSAME( TRANS, 'T' ).AND. &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) &
     &   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SGEMV .
!
      END
      SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      REAL(RKIND) ::     ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  SGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(RKIND) ::     ZERO
      PARAMETER        ( ZERO = 0.0E+0_RKIND )
!     .. Local Scalars ..
      REAL(RKIND) ::     TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of SGER  .
!
      END
      subroutine sscal(n,sa,sx,incx)
!
      USE prec_rkind
      implicit none
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(RKIND) :: sa,sx(*)
      integer i,incx,m,mp1,n,nincx, nel
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 nel=5
      m = mod(n,nel)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
      subroutine sswap (n,sx,incx,sy,incy)
!
      USE prec_rkind
      implicit none
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(RKIND) :: sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n, nel
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 nel = 3
      m = mod(n,nel)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
      SUBROUTINE STBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  STBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(RKIND) ::     ZERO
      PARAMETER        ( ZERO = 0.0E+0_RKIND )
!     .. Local Scalars ..
      REAL(RKIND) ::     TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
     &         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND. &
     &         .NOT.LSAME( TRANS, 'T' ).AND. &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND. &
     &         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STBSV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
     &   RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := inv( A )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT ) &
     &                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT ) &
     &                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT ) &
     &                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT ) &
     &                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := inv( A')*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT ) &
     &               TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT ) &
     &               TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K ) &
     &               KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT ) &
     &               TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT ) &
     &               TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K ) &
     &               KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of STBSV .
!
      END
      SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL(RKIND) ::     ALPHA
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  STRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL(RKIND) ::     TEMP
!     .. Parameters ..
      REAL(RKIND) ::     ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND. &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND. &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT ) &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT ) &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of STRSM .
!
      END


      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993 
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      REAL(RKIND) ::     AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRF computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
      REAL(RKIND) ::     WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMM, SPBTF2, SPOTF2, SSYRK, STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
     &    ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
     &   RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
!
!        Use unblocked code
!
         CALL SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL SPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A12
!
                     CALL STRSM( 'Left', 'Upper', 'Transpose', &
     &                           'Non-unit', IB, I2, ONE, AB( KD+1, I ), &
     &                           LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL SSYRK( 'Upper', 'Transpose', I2, IB, -ONE, &
     &                           AB( KD+1-IB, I+IB ), LDAB-1, ONE, &
     &                           AB( KD+1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
!
!                    Update A13 (in the work array).
!
                     CALL STRSM( 'Left', 'Upper', 'Transpose', &
     &                           'Non-unit', IB, I3, ONE, AB( KD+1, I ), &
     &                           LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     IF( I2.GT.0 ) &
     &                  CALL SGEMM( 'Transpose', 'No Transpose', I2, I3, &
     &                              IB, -ONE, AB( KD+1-IB, I+IB ), &
     &                              LDAB-1, WORK, LDWORK, ONE, &
     &                              AB( 1+IB, I+KD ), LDAB-1 )
!
!                    Update A33
!
                     CALL SSYRK( 'Upper', 'Transpose', I3, IB, -ONE, &
     &                           WORK, LDWORK, ONE, AB( KD+1, I+KD ), &
     &                           LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL SPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A21
!
                     CALL STRSM( 'Right', 'Lower', 'Transpose', &
     &                           'Non-unit', I2, IB, ONE, AB( 1, I ), &
     &                           LDAB-1, AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL SSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
     &                           AB( 1+IB, I ), LDAB-1, ONE, &
     &                           AB( 1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
!
!                    Update A31 (in the work array).
!
                     CALL STRSM( 'Right', 'Lower', 'Transpose', &
     &                           'Non-unit', I3, IB, ONE, AB( 1, I ), &
     &                           LDAB-1, WORK, LDWORK )
!
!                    Update A32
!
                     IF( I2.GT.0 ) &
     &                  CALL SGEMM( 'No transpose', 'Transpose', I3, I2, &
     &                              IB, -ONE, WORK, LDWORK, &
     &                              AB( 1+IB, I ), LDAB-1, ONE, &
     &                              AB( 1+KD-IB, I+IB ), LDAB-1 )
!
!                    Update A33
!
                     CALL SSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
     &                           WORK, LDWORK, ONE, AB( 1, I+KD ), &
     &                           LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
!
  150 CONTINUE
      RETURN
!
!     End of DPBTRF
!
      END

      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL(RKIND) ::     AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite band matrix A using the Cholesky factorization
!  A = U**T*U or A = L*L**T computed by DPBTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangular factor stored in AB;
!          = 'L':  Lower triangular factor stored in AB.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) REAL array, dimension (LDAB,N)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T of the band matrix A, stored in the
!          first KD+1 rows of the array.  The j-th column of U or L is
!          stored in the j-th column of the array AB as follows:
!          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, nel
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           STBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
     &   RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B where A = U'*U.
!
         nel = 1
         DO 10 J = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB, &
     &                  LDAB, B( 1, J ), nel )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, &
     &                  LDAB, B( 1, J ), nel )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO 20 J = 1, NRHS
!
!           Solve L*X = B, overwriting B with X.
!
            CALL STBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, &
     &                  LDAB, B( 1, J ), 1 )
!
!           Solve L'*X = B, overwriting B with X.
!
            CALL STBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB, &
     &                  LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of DPBTRS
!
      END

      SUBROUTINE SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      REAL(RKIND) ::     AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  SPBTF2 computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix, U' is the transpose of U, and
!  L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of super-diagonals of the matrix A if UPLO = 'U',
!          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U'*U or A = L*L' of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, KLD, KN
      REAL(RKIND) ::     AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSCAL, SSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
     &   RETURN
!
      KLD = MAX( 1, LDAB-1 )
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) &
     &         GO TO 30
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
!
!           Compute elements J+1:J+KN of row J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL SSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL SSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD, &
     &                    AB( KD+1, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) &
     &         GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
!
!           Compute elements J+1:J+KN of column J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL SSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL SSYR( 'Lower', KN, -ONE, AB( 2, J ), 1, &
     &                    AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      RETURN
!
   30 CONTINUE
      INFO = J
      RETURN
!
!     End of SPBTF2
!
      END
      SUBROUTINE SPOTF2( UPLO, N, A, LDA, INFO )
!
      USE prec_rkind
      implicit none
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RKIND) ::     ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      REAL(RKIND) ::     AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL(RKIND) ::     SDOT
      EXTERNAL           LSAME, SDOT
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPOTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
     &   RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - SDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
            IF( J.LT.N ) THEN
               CALL SGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
     &                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL SSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - SDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), &
     &            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
            IF( J.LT.N ) THEN
               CALL SGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
     &                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL SSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of SPOTF2
!
      END

      function sdot(n,sx,incx,sy,incy)
!
      USE prec_rkind
      implicit none
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      REAL(RKIND) :: sx(*),sy(*),stemp, sdot
      integer i,incx,incy,ix,iy,m,mp1,n, nel
!
      stemp = 0.0e0_RKIND
      sdot = 0.0e0_RKIND
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 nel=5
      m = mod(n,nel)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
     &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end

      SUBROUTINE SSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      REAL(RKIND) ::     ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  SSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(RKIND) ::     ZERO
      PARAMETER        ( ZERO = 0.0E+0_RKIND )
!     .. Local Scalars ..
      REAL(RKIND) ::     TEMP
      INTEGER            I, INFO, IX, J, JX, KX
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND. &
     &         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYR  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
     &   RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF( LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when A is stored in upper triangle.
!
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
!
!        Form  A  when A is stored in lower triangle.
!
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SSYR  .
!
      END

      SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, &
     &                   BETA, C, LDC )
!
      USE prec_rkind
      implicit none
!     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      REAL(RKIND) ::     ALPHA, BETA
!     .. Array Arguments ..
      REAL(RKIND) ::     A( LDA, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  SSYRK  performs one of the symmetric rank k operations
!
!     C := alpha*A*A' + beta*C,
!
!  or
!
!     C := alpha*A'*A + beta*C,
!
!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
!  in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
!
!              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrix  A.  K must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL(RKIND) ::     TEMP
!     .. Parameters ..
      REAL(RKIND) ::     ONE ,         ZERO
      PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
!
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND. &
     &         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND. &
     &         ( .NOT.LSAME( TRANS, 'T' ) ).AND. &
     &         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYRK ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( N.EQ.0 ).OR. &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  C := alpha*A*A' + beta*C.
!
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
!
!        Form  C := alpha*A'*A + beta*C.
!
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SSYRK .
!
      END
