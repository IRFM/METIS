function [x,y] = separatrice_mtfile(geom,mode)
%  #ifdef MATLAB73
%  #include "fintrf.h"
%  #endif
%        subroutine mexFunction(nlhs,plhs,nrhs,prhs)
%  C Compilation : mex separatrice.f -lnag
%  C--------------------------------------------------------------
%  
%        IMPLICIT REAL*8 (A-H,O-Z)
%        IMPLICIT integer (I-N)
%  C DECLARATION DES VARIABLES
%  c pointeurs mex   
%        CRONOSINT plhs(*), prhs(*),mxGetPr, mxCreateFull
%        integer nlhs, nrhs, mxgetN, mxgetM
%  c 
%  c Extrait du main de helios.f concernant le trace de la derniere surface
%  c magnetique (il s'agit d'une implementation triviale des formules fournies)
%  c
%  c pointeurs et tailles      
%        CRONOSINT geom_pr,mode_pr,x_pr,y_pr
%        integer*4 ifail 
%        real*8 mode,geom(8),x(125*4),y(125*4)
%        data pi/3.141592653589793/
two = 2.;
one = 1.;
rho = 1.; 
rkappax1  = 1.687;
ptdeltax1 = 0.466;
psip1dg   = 0.;
psim1dg   = 0.;
rkappax2  = 2.001;
ptdeltax2 = 0.568;
psip2dg   = 22.46;
psim2dg   = 67.92;

%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%  
%  C---------------------------------------------------------------------
%  C GESTION DES ENTREES     
%  C entree
%        write(*,*) 'debut sepratrice'
%        geom_pr = mxGetPr(prhs(1))
%        call mxCopyPtrToReal8(geom_pr,geom,8)

rkappax1  = geom(1);
ptdeltax1 = geom(2);
psip1dg   = geom(3);
psim1dg   = geom(4);
rkappax2  = geom(5);
ptdeltax2 = geom(6);
psip2dg   = geom(7);
psim2dg   = geom(8);

%  C mode de calcul
%        mode_pr = mxGetPr(prhs(2))
%        call mxCopyPtrToReal8(mode_pr,mode,1)


%c Initializing geometrical parameters
psimoins1=psim1dg*pi/180.;
psiplus1=psip1dg*pi/180.;
psimoins2=psim2dg*pi/180.;
psiplus2=psip2dg*pi/180.;
tmoins1=tan(psimoins1)*(1.-ptdeltax1)/rkappax1;
tmoins2=tan(psimoins2)*(1.-ptdeltax2)/rkappax2;
%c      write(*,*) 'tmoins1 = ',tmoins1,' tmoins2 =',tmoins2
tplus1=tan(psiplus1)*(1.+ptdeltax1)/rkappax1;
tplus2=tan(psiplus2)*(1.+ptdeltax2)/rkappax2;
%c      write(*,*) 'tplus1 = ',tplus1,' tplus2 =',tplus2
%  c      
%  c Activer manuellement la bonne partie du programme suivant le cas (voir texte
%  c fourni) en remplacant if(one.gt.two) par if(two.gt.one)
%  c Un algorithme automatique peut aussi �re �rit
%  c
%  c plot of the plasma last closed magnetic surface
%  c
%  c case where the upper part is made of 2 portions of ellipse and the
%  c bottom part is made of 2 portions of ellipse
%  c
if (mode  == 1)
    tetaxp2   = 0.5*pi-asin(tplus2/(1.-tplus2));
    tetap2in  = 0.;
    tetap2fi  = tetaxp2;
    nbtetap2  = 125;
    tetap2st  = (tetap2fi-tetap2in)/(nbtetap2-1);
    tetaxm2   = 0.5*pi+asin(tmoins2/(1.-tmoins2));
    tetam2in  = tetaxm2;
    tetam2fi  = pi;
    nbtetam2  = 125;
    tetam2st  = (tetam2fi-tetam2in)/(nbtetam2-1);
    tetaxm1   = 0.5*pi+asin(tmoins1/(1.-tmoins1));
    tetam1in  = pi;
    tetam1fi  = tetaxm1;
    nbtetam1  = 125;
    tetam1st  = (tetam1fi-tetam1in)/(nbtetam1-1);
    tetaxp1   = 0.5*pi-asin(tplus1/(1.-tplus1));
    tetap1in  = -tetaxp1;
    tetap1fi  = 0.;
    nbtetap1  = 125;
    tetap1st  = (tetap1fi-tetap1in)/(nbtetap1-1);
    ii        = 1;
%c bottom plus
    for noteta=1:nbtetap2
      teta     = tetap2in+(noteta-1)*tetap2st;
      alpha0p2 = -(ptdeltax2+(1.-ptdeltax2)*tplus2)/(1.-2.*tplus2);
      alphap2  = (1.+ptdeltax2)*(1.-tplus2)/(1.-2.*tplus2);
      betap2   = rkappax2*(1.-tplus2)/(1.-2.*tplus2)^0.5;
      x(ii)    = rho*(alpha0p2+alphap2*cos(teta));
      y(ii)    = -rho*(betap2*sin(teta));
      ii       = ii+1;
    end   %noteta
%c bottom minus
    for noteta=1:nbtetam2
      teta     = tetam2in+(noteta-1)*tetam2st;
      alpha0m2 = -(ptdeltax2-(1.+ptdeltax2)*tmoins2)/(1.-2.*tmoins2);
      alpham2  = (1.-ptdeltax2)*(1.-tmoins2)/(1.-2.*tmoins2);
      betam2   = rkappax2*(1.-tmoins2)/(1.-2.*tmoins2)^0.5;
      x(ii)    = rho*(alpha0m2+alpham2*cos(teta));
      y(ii)    = -rho*(betam2*sin(teta));
%  c          if (ii .eq. 150) then
%  c            write(*,*) 'y(ii)=',y(ii),betam2,teta,sin(teta)
%  c          endif

      ii       = ii+1;
    end  %noteta
%c upper minus
    for noteta=1:nbtetam1
      teta     = tetam1in+(noteta-1)*tetam1st;
      alpha0m1 = -(ptdeltax1-(1.+ptdeltax1)*tmoins1)/(1.-2.*tmoins1);
      alpham1  = (1.-ptdeltax1)*(1.-tmoins1)/(1.-2.*tmoins1);
      betam1   = rkappax1*(1.-tmoins1)/(1.-2.*tmoins1)^0.5;
      x(ii)    = rho*(alpha0m1+alpham1*cos(teta));
      y(ii)    = rho*(betam1*sin(teta));
      ii       = ii+1;
%c      write(7,*)x,y
    end   %noteta
%c upper plus
    for noteta=1:nbtetap1
      teta     = tetap1in+(noteta-1)*tetap1st;
      alpha0p1 = -(ptdeltax1+(1.-ptdeltax1)*tplus1)/(1.-2.*tplus1);
      alphap1  = (1.+ptdeltax1)*(1.-tplus1)/(1.-2.*tplus1);
      betap1   = rkappax1*(1.-tplus1)/(1.-2.*tplus1)^0.5;
      x(ii)    = rho*(alpha0p1+alphap1*cos(teta));
      y(ii)    = -rho*(betap1*sin(teta));
      ii       = ii+1;
%c      write(7,*)x,y
    end   %noteta
end  %true/false

%  c
%  c case where the upper part is made of 2 portions of ellipse and the
%  c bottom part is made of 1 portion of hyperbola (inside) and
%  c 1 portion of ellipse (outside)
if(mode == 2)
    tetaxp1=0.5*pi-asin(tplus1/(1.-tplus1));
    tetap1in=0.;
    tetap1fi=tetaxp1;
    nbtetap1=125;
    tetap1st=(tetap1fi-tetap1in)/(nbtetap1-1);
    tetaxm1=0.5*pi+asin(tmoins1/(1.-tmoins1));
    tetam1in=tetaxm1;
    tetam1fi=pi;
    nbtetam1=125;
    tetam1st=(tetam1fi-tetam1in)/(nbtetam1-1);
    phixm2 = dasinh((2.*tmoins2-1.)^0.5/(1.-tmoins2));
%c        phixm2=asinh((2.*tmoins2-1.)^0.5/(1.-tmoins2));
    phim2in=0.;
    phim2fi=-phixm2;
    nbphim2=125;
    phim2st=(phim2fi-phim2in)/(nbphim2-1);
    tetaxp2=0.5*pi-asin(tplus2/(1.-tplus2));
    tetap2in=-tetaxp2;
    tetap2fi=0.;
    nbtetap2=125;
    tetap2st=(tetap2fi-tetap2in)/(nbtetap2-1);
%c upper plus
    ii = 1;
    for noteta=1:nbtetap1
      teta=tetap1in+(noteta-1)*tetap1st;
      alpha0p1=-(ptdeltax1+(1.-ptdeltax1)*tplus1)/(1.-2.*tplus1);
      alphap1=(1.+ptdeltax1)*(1.-tplus1)/(1.-2.*tplus1);
      betap1=rkappax1*(1.-tplus1)/(1.-2.*tplus1)^0.5;
      x(ii) =rho*(alpha0p1+alphap1*cos(teta));
      y(ii) =rho*(betap1*sin(teta));
      ii = ii+1
    end  %noteta
%c upper minus
    for noteta=1:nbtetam1
      teta=tetam1in+(noteta-1)*tetam1st;
      alpha0m1=-(ptdeltax1-(1.+ptdeltax1)*tmoins1)/(1.-2.*tmoins1);
      alpham1=(1.-ptdeltax1)*(1.-tmoins1)/(1.-2.*tmoins1);
      betam1=rkappax1*(1.-tmoins1)/(1.-2.*tmoins1)^0.5;
      x(ii)=rho*(alpha0m1+alpham1*cos(teta));
      y(ii)=rho*(betam1*sin(teta));
      ii = ii+1;
    end %noteta
%c bottom minus
    for nophi=1:nbphim2
      phi=phim2in+(nophi-1)*phim2st;
      alpha0m2=-((1.+ptdeltax2)*tmoins2-ptdeltax2)/(2.*tmoins2-1.);
      alpham2=(1.-ptdeltax2)*(1.-tmoins2)/(2.*tmoins2-1.);
      betam2=rkappax2*(1.-tmoins2)/(2.*tmoins2-1.)^0.5;
      x(ii)=rho*(alpha0m2+alpham2*cosh(phi));
      y(ii)=rho*(betam2*sinh(phi));
      ii=ii+1;
    end  %nophi
%c bottom plus
    for noteta=1:nbtetap2
      teta=tetap2in+(noteta-1)*tetap2st;
      alpha0p2=-(ptdeltax2+(1.-ptdeltax2)*tplus2)/(1.-2.*tplus2);
      alphap2=(1.+ptdeltax2)*(1.-tplus2)/(1.-2.*tplus2);
      betap2=rkappax2*(1.-tplus2)/(1.-2.*tplus2)^0.5;
      x(ii)=rho*(alpha0p2+alphap2*cos(teta));
      y(ii)=rho*(betap2*sin(teta));
      ii = ii+1;
    end  %noteta
end   %true/false
%c
%  c case where the upper part is made of 1 portion of hyperbola and
%  c 1 portion of ellipse
%  c and the bottom part of 2 portions of ellipse
%  c
if(mode == 3)
    tetaxp2=0.5*pi-asin(tplus2/(1.-tplus2));
    tetap2in=0.;
    tetap2fi=tetaxp2;
    nbtetap2=125;
    tetap2st=(tetap2fi-tetap2in)/(nbtetap2-1);
    tetaxm2=0.5*pi+asin(tmoins2/(1.-tmoins2));
    tetam2in=tetaxm2;
    tetam2fi=pi;
    nbtetam2=125;
    tetam2st=(tetam2fi-tetam2in)/(nbtetam2-1);
    phixm1=dasinh((2.*tmoins1-1.)^0.5/(1.-tmoins1)) ;
%c         phixm1=asinh((2.*tmoins1-1.)^0.5/(1.-tmoins1))
%c         write(*,*) 'phixm1=',phixm1
    phim1in=0.;
    phim1fi=-phixm1;
    nbphim1=125;
    phim1st=(phim1fi-phim1in)/(nbphim1-1);
    tetaxp1=0.5*pi-asin(tplus1/(1.-tplus1));
    tetap1in=-tetaxp1;
    tetap1fi=0.;
    nbtetap1=125;
    tetap1st=(tetap1fi-tetap1in)/(nbtetap1-1);
    ii = 1;
%c bottom plus
    for noteta=1:nbtetap2
      teta=tetap2in+(noteta-1)*tetap2st;
      alpha0p2=-(ptdeltax2+(1.-ptdeltax2)*tplus2)/(1.-2.*tplus2);
      alphap2=(1.+ptdeltax2)*(1.-tplus2)/(1.-2.*tplus2);
      betap2=rkappax2*(1.-tplus2)/(1.-2.*tplus2)^0.5;
      x(ii)=rho*(alpha0p2+alphap2*cos(teta));
      y(ii)=-rho*(betap2*sin(teta));
      ii = ii+1;
    end %noteta
%c bottom minus
    for noteta=1:nbtetam2
      teta=tetam2in+(noteta-1)*tetam2st;
      alpha0m2=-(ptdeltax2-(1.+ptdeltax2)*tmoins2)/(1.-2.*tmoins2);
      alpham2=(1.-ptdeltax2)*(1.-tmoins2)/(1.-2.*tmoins2);
      betam2=rkappax2*(1.-tmoins2)/(1.-2.*tmoins2)^0.5;
      x(ii)=rho*(alpha0m2+alpham2*cos(teta));
      y(ii)=-rho*(betam2*sin(teta));
      ii = ii+1;
    end  %noteta
% upper minus
    for nophi=1:nbphim1
      phi=phim1in+(nophi-1)*phim1st;
      alpha0m1=-((1.+ptdeltax1)*tmoins1-ptdeltax1)/(2.*tmoins1-1.);
      alpham1=(1.-ptdeltax1)*(1.-tmoins1)/(2.*tmoins1-1.);
      betam1=rkappax1*(1.-tmoins1)/(2.*tmoins1-1.)^0.5;
      x(ii)=rho*(alpha0m1+alpham1*cosh(phi));
      y(ii)=-rho*(betam1*sinh(phi));
      ii = ii+1;
    end %nophi
%c upper plus
    for noteta=1:nbtetap1
      teta=tetap1in+(noteta-1)*tetap1st;
      alpha0p1=-(ptdeltax1+(1.-ptdeltax1)*tplus1)/(1.-2.*tplus1);
      alphap1=(1.+ptdeltax1)*(1.-tplus1)/(1.-2.*tplus1);
      betap1=rkappax1*(1.-tplus1)/(1.-2.*tplus1)^0.5;
      x(ii)=rho*(alpha0p1+alphap1*cos(teta));
      y(ii)=-rho*(betap1*sin(teta));
      ii=ii+1;
    end  %noteta
end %true/false
%  c
%  c case of 2 semi-ellipse with triangularity
%  c 
if (mode == 4)
    rkappax=0.5*(rkappax1+rkappax2);
    ptdeltax=0.5*(ptdeltax1+ptdeltax2);
    dzeta0=0.5*(rkappax1-rkappax2);
%c external part
    tetain=-0.5*pi;
    tetafi=0.5*pi;
    nbteta=250;
    tetast=(tetafi-tetain)/(nbteta-1);
    ii = 1;
    for noteta=1:nbteta
      teta=tetain+(noteta-1)*tetast;
      x(ii)=-ptdeltax+(1.+ptdeltax)*cos(teta);
      y(ii)=dzeta0+rkappax*sin(teta);
      ii = ii+1;
    end   %noteta
%c internal part
    tetain=0.5*pi;
    tetafi=1.5*pi;
    nbteta=250;
    tetast=(tetafi-tetain)/(nbteta-1);
    for noteta=1:nbteta
      teta=tetain+(noteta-1)*tetast;
      x(ii)=-ptdeltax+(1.-ptdeltax)*cos(teta);
      y(ii)=dzeta0+rkappax*sin(teta);
      ii = ii +1; 
    end  %noteta
end     %true/false
%         sizeOUT = 125*4
%  C x
%        write(*,*) 'debut ecriture de x et y'
%        plhs(1) = mxCreateFull(125*4,1,0)
%        x_pr = mxGetPr(plhs(1))
%        call mxCopyReal8ToPtr(x,x_pr,125*4)
%  C x
%        plhs(2) = mxCreateFull(125*4,1,0)
%        y_pr = mxGetPr(plhs(2))
%        call mxCopyReal8ToPtr(y,y_pr,125*4)
%  
%         end
%  c
x=x(:);
y = y(:);
