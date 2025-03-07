%      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
%   the mex function should have the name of the file. In this case interpos
%   Then typically in matlab one calls:
%   >> [a,b,c,..]=interpos(arg1,arg2,...)
%
%   nlhs, nrhs: number of input arguments on left and right-hand side respectively
%   thus one can test and depending on number of arguments choose different options
%
%   plhs, prhs: pointer to the different arguments: plhs(i), i=1,nlhs
%   (plhs and prhs are integer arrays)
%
%   fmex function: mxGetPr gets the pointer value in prhs, plhs
%   .              mxCopyPtrToReal8: copy pointed values to local array
%   .              mxCopyReal8ToPtr: copy local array values to pointed array
%   .              mxGetM: get nb of rows in pointed array
%   .              mxGetN: get nb of columns in pointed array
%   .              mxCreateFull: creates a matrix for the return argument
%
%   In this example, matlab call expected:
%
%   >> [yout{,youtp,youtpp,youtint}] = interpos({kopt,}xin,yin{,xout}{,taus{,nbc,ybc}{,sigma}})
%
%   where {..} means facultative, thus one can have for example the following calls:
%
%   >> [yout]=interpos(xin,yin,yout)
%   >> [yout]=interpos(kopt,xin,yin,xout)
%   >> [yout]=interpos(xin,yin,taus)
%   >> [yout,youtp]=interpos(xin,yin)
%   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,-1.)
%   >> [yout,youtp,youtpp]=interpos(xin,yin,xout,taus,nbc,ybc)
%   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,taus,nbc,ybc)
%   >> [yout,youtp,youtpp,youtint]=interpos(xin,yin,xout,taus,nbc,ybc,sigma)
%
%   NOTE: If xout is a single point, taus should also be given otherwise it assumes it is a taus value:
%          interpos(xin,yin,xxx) means taus=xxx if xxx size of 1
%          interpos(xin,yin,xxx,0.) means xout=xxx, taus=0. if xxx size of 1
%
%   The arguments have the following meaning with respect to the function y(x)
%   that one wants to interpolate/extrapolate or compute its derivatives and integrals:
%
%   kopt: option for interpolation and extrapolation method (default is kopt=13):
%   .     = 1, 11 or 21 : linear interpolation
%   .     = 2, 12 or 22 : quadratic "
%   .     = 3, 13 or 23 : cubic spline interpolation, with tension=taus
%   .     < 10 : send warning if need to extrapolate
%   .     in [11,19] : extrapolate using a lower order than for interpolation
%   .     in [21,29] : extrapolate using same order as interpolation
%
%   xin : array giving the input values of x
%   yin : array giving the input values of y(x=xin(i))
%
%   xout: array of values of x at which the function and/or its derivatives
%   .     are to be computed. If xout is not given, assumes xout=xin
%   yout: interpolated values of y at x=xout(i).
%   youtp: 1st derivative of y at x=xout(i)
%   youtpp: 2nd derivative of y at x=xout(i)
%   youtint: Integral of (y dx) from x=xin(1) to x=xout(i)
%
%   taus  : tension value for cubic spline interpolation. If not given, uses taus=0
%           if taus < 0 (typically -1) uses a default taus value: taus=abs(taus)*default_taus
%           (the default_taus is typically min(delta_x)**3)
%
%   nbc(2): [NBCLFT NBCRGT] (default: [0 0])
%     BOUNDARY CONDITIONS, 4 TYPES DETERMINED BY THE VALUE OF (nbc(1)=NBCLFT
%     and nbc(2) for left and right-hand side BC.):
%
%     0) VALUE OF SECOND DERIVATIVE AT XBCLFT OR RGT IS GIVEN (0 OR 10)
%     1) VALUE OF 1ST        "       "   "     "  "   "   "   (1 OR 11)
%     2) VALUE OF FUNCTION AT XBCLFT OR RGT IS GIVEN          (2 OR 12)
%     if nbc(1)=-1, then assumes periodic boundary condition and uses ybc(1) for the period
%
%     The value of nbc(1 or 2) should be >= 10 if the BC is not at an end point:
%     XBCLFT~=xin(1) or XBCRGT~=xin(end)
%
%     Examples:
%       A good value for radial profiles in rho between [0,1] is to specify
%       the first derivative = 0 at left and second derivative=0 at right:
%       => nbc = [1 0] and ybc=[0. 0.]
%
%       or to obtained the first derivative at right-hand side from a
%       lagrangian interpolation of the last points:
%       => nbc=[1 1] and ybc=[0. 1e32]
%
%       Note: The B.C. can only be given anywhere but within the interval [xin(1),xin(end)]
%       Therefore, say rho=[0.1 ... 1] and one wants to impose zero
%       derivative at rho=0, one should add a point in the input before. One uses sigma to avoid forcing the value
%       calling interpos:
%             rho_eff(1)=0.;
%             rho_eff(2:length(rho)+1) = rho;
%       then  [yout]=interpos([0.;rho],[yin(1);yin],..,taus,[1 0],[0. 0.],[1000;ones(size(rho))]);
%
%   ybc(2 or 4 or 6) (default: [0 0]): [YBCLFT YBCRGT] or [YBCLFT YBCRGT XBCLFT XBCRGT] with
%     THE VALUE IS GIVEN BY YBCLFT OR YBCRGT RESPECTIVELY.
%
%     FOR nbc type 1: IF (YBCLFT OR YBCRGT > 1E31 THEN DER. FROM LAGRANGIAN INTERP.
%     FOR nbc type 1: IF (YBCLFT OR YBCRGT <-1E31 THEN DER. FROM LINEAR     INTERP.
%
%     IF NBCLFT OR NBCRGT IS < 10, PXIN(1) OR PXIN(KNIN) IS USED INSTEAD
%     OF XBCLFT OR XBCRGT, RESPECTIVELY => XBCLFT OR XBCRGT NOT USED
%
%   ybc(6): [YBCLFT YBCRGT XBCLFT XBCRGT PXEXP0 PXEXPDL]
%        Enables one to specify a gaussian wight to taus with PXEXP0 and PXEXPDL:
%              taus_eff =  taus * EXP(-((XIN-PXEXP0)/PXEXPDL)**2)
%        IF PXEXP0 NOT IN [PXIN(1),PXIN(KNIN)], EXP() IGNORED AND taus_eff=cst= taus
%        if ybc(5:6) not given, then  PXEXP0=xin(1)-1. and pxexpdl=1. to get constant taus
%
%   ybc(1)=period for periodic boundary condition, if nbc(1)=-1. Defines the condition that y(x+period)=y(x)
%          In this case xin(end) should not be equal to xin(1)+period, since it is redundant. However interpos just
%          does not use the end point in this case.
%        if period=-1, assumes period=xin(end)-xin(1). Useful if xin goes from 0 to 2pi for example
%
%    sigma: error_bar at each (x,y) point used for the fit. The effective taus value will then
%           be taus(i)=taus .* sigma(i) ./ min(sigma). So uses the relative values of sigma.
%

