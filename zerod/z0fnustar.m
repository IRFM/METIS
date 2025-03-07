function [nustare,nustari,tei]=z0fnustar(te,ti,dens,q,r,a,zeff,zi)

% calculation of normalised frequency
% this function is mainly used for the JET and ITPA database 
%  -- INPUTS --
%  te  -- electron temperature [eV]
%  ti  -- ion temperature [eV]
%  dens -- electron density [m-3]
%  q -- q profile
%  r -- major radius (m)
%  a -- minor radius  (m)
%  zeff -- effective charge
%  zi  -- ion charge (Zi=1; deuterium plasma)
% -- OUTPUTS 
%-- nustare and nustari normalised electron/ion collisionality (connection length/trapped particle mean-fre path)  
%- cf ITER definition in ITER physics basis Nuc. Fusion vol.39 page 2203
%	X. Litaudon 03/2003 



%--constant 
qe=1.602e-19;eps0=8.8542e-12;depi=2*pi;me=9.109534e-31;

%--ln coulombien
lnc=31.474 + log(dens.^(-.5).*te);

%--inverse aspect ratio 
eps=a./r;

%thermal speed
vte=((qe*te)/me).^(.5);

%tei  collision  time 
tei=(3*me^(0.5)*eps0*eps0*(depi*qe*te).^1.5)./((qe^4)*zeff.*lnc.*dens);

%nustare 
nustare=r.*q./((eps.^1.5).*vte.*tei);

%-normalised nustar_i
lnci=1.1*lnc; %cf Weysson page 663
vti=((qe*ti)).^(.5);%be careful vti is multiplied by sqrt(m)
%tii=(3*eps0*eps0*(depi)^(3/2)*(mi^(0.5)))/((q^4)*zeff);
tii=(3*sqrt(2)*eps0*eps0*(depi*qe*ti).^1.5)./((qe^4)*zi.*zi.*zeff.*lnci.*dens);% be careful tii is divided by sqrt(m)
nustari=r.*q./((eps.^1.5).*vti.*tii);%vti*tii sqrt(m) does not appear

