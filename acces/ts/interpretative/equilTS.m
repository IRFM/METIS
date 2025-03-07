%  EQUILTS calcul analitique de l'equilibre pour Tore-Supra
%----------------------------------------------------------
% fichier equilTS.m ->  equilTS
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule l'equilibre pour Tore-Supra. Elle
% utilise des formules analytique. Elle suppose que la section du plasma
% est circulaire, que le rapport a/r0 est petit et que le beta est pas
% trop grand.  
%  
% syntaxe du mode 0 :
% 
%    [rhophi,rhog,d,vpr,grho2,r2i,ri,grho2r2,psi,phi,fdia, ...
%     ptot,jmoy,q,bpol,ftrap,b2,b2i,sp,dpsidrho,tau,sqrtg, ...
%     rmoy,r2moy,r2tau2,r3invtau3,r3invtau] = equilTS(r0,a,ip,beli,wdia,b0);
% 
% entrees du mode 0 :  
%      
%     r0   = grand rayon (m)   [nbt,1]
%     a    = petit rayon (m)   [nbt,1]
%     ip   = courant plasma (MA)    [nbt,1]
%     beli = beta_dia + li/2        [nbt,1]
%     wdia = contenu enregetique total du plama (MJ)    [nbt,1]
%     b0   = champ toroidal (T)  [nbt,1]
%         
% syntaxe du mode 1 :
% 
%    function [rho,rhog,d,vpr,grho2,r2i,ri,grho2r2,psi,phi,fdia, ...
%              ptot,jmoy,q,bpol,ftrap,b2,b2i,sp,dpsidrho,tau,sqrtg, ...
%              rmoy,r2moy,r2tau2,r3invtau3,r3invtau] =  ...
%              equilTS(r0,a,d0,b0,rhoe,jmoy,ptot,{mode_rho});
%         
% entrees du mode 1 :
% 
%     r0   = grand rayon (m)    [nbt,1]
%     a    = petit rayon (m)    [nbt,1]
%     d0   = decentrement de Shafranov au centre (m)    [nbt,1]
%            ou profil de decentrement (m) [nbt,nbrho]
%     b0   = champ toroidal (T)      [nbt,1]      
%     rhoe = coordonnees normalisee    [nbt,nbrho_in] :
%               soit r/a  (mode_rho = 0)
%               soit sqrt(Phi)/sqrt(Phi(a))  (mode_rho = 1)
%     jmoy = profil de densite de courant (A*m^-2)   [nbt,nbrho_in] 
%     ptot = profil de pression totale (keV*m^-3)  [nbt,nbrho_in]
%            ou wdia [nbt,1] (MJ)
%
% sorties :
% 
%      rhophi     =  sqrt(Phi / R B0)  coordonnee de flux toroidal (m)   [nbt,nbrho]
%      rhog       =  r   coordonnee geometrique (m) [nbt,nbrho]
%      d          =  decentrement de Shafranov (m) [nbt,nbrho]
%      vpr        =  covolume dV/drho  (m^2) [nbt,nbrho] 
%      grho2      =  <|gradient(rho)|^2> [nbt,nbrho] 
%      r2i        =  <1/R^2> [nbt,nbrho] 
%      ri         =  <1/R> [nbt,nbrho] 
%      grho2r2    =  <|gradient(rho)|^2/R^2> [nbt,nbrho] 
%      psi        =  flux poloidal [nbt,nbrho] 
%      phi        =  flux toroidal (approximation Fdia = R0 B0) [nbt,nbrho] 
%      fdia       =  fonction diamagnetique [nbt,nbrho] 
%      ptot       =  profil de pression totale (keV*m^-3) [nbt,nbrho] 
%      jmoy       =  profil de densite de courant <Jphi/R>/<1/R> (A*m^-2) [nbt,nbrho] 
%      q          =  profil de facteur de securite  [nbt,nbrho] 
%      bpol       =  <Bpoloidal> (T)  [nbt,nbrho] 
%      ftrap      =  fraction de particules piegees  [nbt,nbrho] 
%      b2         =  <B^2> (T^2)  [nbt,nbrho] 
%      b2i        =  <1/B^2> (T^-2)  [nbt,nbrho] 
%      sp         =  perimetre de la surface magnetique (2*pi*a*rho, m)  [nbt,nbrho] 
%      dpsidrho   =  derivee du flux poloidal [nbt,nbrho] 
%      tau        =  <tau>
%      sqrtg      =  <sqrt(g)>
%      rmoy       =  <R>
%      r2moy      =  <R^2>
%      r2tau2     =  <R^2*tau^2>
%      r3invtau3  =  <R^3?tau^3>
%      r3invtau   =  <R^3/tau>
%      
% remarque : si b0 <0 alors les calcul des moyennes de B qui sont long ne sont pas fait et b0 =abs(b0)
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 24/02/2000.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [rhophi,rhog,d,vpr,grho2,r2i,ri,grho2r2,psi,phi,fdia, ...
          ptot,jmoy,q,bpol,ftrap,b2,b2i,sp,dpsidrho,tau,sqrtg, ...
          rmoy,r2moy,r2tau2,r3invtau3,r3invtau]= equilTS(varargin);

         
if length(varargin) <6 
	error('nombre d''arguments incorrect !')
end
if length(varargin) < 7
	mode    = 0;
	r0      = varargin{1};
	a       = varargin{2};
	ip      = varargin{3};
	beli    = varargin{4};
	wdia    = varargin{5};
	b0      = varargin{6};
else
	mode    = 1;
	r0      = varargin{1};
	a       = varargin{2};
	d0      = varargin{3};
	b0      = varargin{4};
	rhoe    = varargin{5};
	jmoy    = varargin{6};
	ptot    = varargin{7};
	if length(varargin) ==8
		mode_rho = varargin{8};
	else
		mode_rho = 1;
	end
	
end

if any(b0<0)
	b0=abs(b0);
	modeb2=0;
else
	modeb2=1;
end

warning off        
         
% les constante
ee          =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)

% changement de dimension
nbrho = 101;
vt=ones(size(a));
ve=ones(1,nbrho);
r0=r0*ve;
a=a*ve;
if mode ==0
	beli=beli*ve;
	wdia=wdia*ve .* 1e6;
	ip=ip*ve .* 1e6;
elseif size(ptot,2)==1
        wdia=ptot*ve .* 1e6;	
end
if length(b0)==1
	b0=b0.* vt;
end
b0=b0*ve;

% 1 - les coordonnees geometriques
rho   = linspace(0,1,nbrho);
r2m1  = 1 - rho .^2;
rho   = vt * rho;
r2m1  = vt * r2m1;

% 2 - separation des parametres et calcul de d0
if mode == 0
	beta     = wdia ./ ((3/8) .*mu0 .* r0 .* ip .^2);
	li       = 2 .* (beli - beta);
	ra    = a ./ r0;
	d0 = a .^ 2 ./ r0 ./ 2 .* beli .* (1 - ra);
   d  = d0 .* r2m1;
elseif size(d0,2) ==1
	d0 = d0 *ve ;
	d  = d0 .* r2m1;
else
	if mode_rho == 0
		d = tsplinet(rhoe,d0,rho);
	else
		warning('Ce mode ne marche pas -> valeur approchee de d');
		d = tsplinet(rhoe,d0,rho);
	end
end

% 3 - calcul de la derivee de d 
dp = pdederive(a .* rho,d,0,2,2,1);

% 4 - coefficient geometrique 
vpr = (4*pi^2) .* a .* rho .* (r0 + d + a.* rho .* dp ./ 2);

sp  = 2 .* pi .* a .* rho;

r2i = (4*pi^2) ./ vpr .* (dp + (a .* rho - (r0 + d ) .* dp) ./ ...
       sqrt((r0 +d) .^ 2 - (a .* rho) .^2 ));
r2i(:,1) = 2 .* r2i(:,2) - r2i(:,3);
       
ri = (4*pi^2) ./ vpr .* a .* rho;
ri(:,1) = 2 .* ri(:,2) - ri(:,3);

grho2r2 = (4*pi^2) ./ vpr .* (a .* rho ./ sqrt((r0 +d) .^ 2 - (a .* rho) .^2 ) - ...
          dp ./ sqrt(1 - dp .^ 2)) .* a .* rho ./ (a .* rho - (r0 + d ) .* dp);    
grho2r2(:,1) = 2 .* grho2r2(:,2) - grho2r2(:,3);
       
grho2 = (4*pi^2) ./ (vpr .* dp) .* a .* rho .* (a .* rho +  ...
        ((r0 + d) .* dp - a .* rho) ./ sqrt(1 - dp .^ 2));  
grho2(:,1) = 2 .* grho2(:,2) - grho2(:,3);


arho =a .*rho;
qq =1;
tau       = 2.*pi.^2.*arho.^2.*(r0.*dp.^2+2.*d+2.*r0+d.*dp.^2+2.*arho.*dp)./vpr;
tau(:,1) = 2 .* tau(:,2) - tau(:,3);             
sqrtg     = 1./2.*pi.^2.*arho.^2.*(8.*r0.^2+4.*arho.^2+8.*d.^2+3.*arho.^2.*dp.^2+ ...
            16.*d.*arho.*dp+16.*r0.*d+16.*arho.*r0.*dp+8.*r0.*d.*dp.^2+4.*d.^2.*dp.^2+4.*r0.^2.*dp.^2)./vpr;
sqrtg(:,1) = 2 .* sqrtg(:,2) - sqrtg(:,3);             
rmoy      = 2.*pi.^2.*arho.*(arho.^2+2.*r0.^2+2.*d.^2+4.*r0.*d+2.*arho.*r0.*dp+2.*d.*arho.*dp)./vpr;
rmoy(:,1) = 2 .* rmoy(:,2) - rmoy(:,3);             
r2moy     = 1./2.*pi.^2.*arho.*(12.*d.*arho.^2+24.*r0.^2.*d+12.*r0.^2.*arho.*dp+8.*r0.^3+ ...
            24.*r0.*d.^2+3.*arho.^3.*dp+12.*r0.*arho.^2+8.*d.^3+12.*d.^2.*arho.*dp+24.*r0.*d.*arho.*dp)./vpr;
r2moy(:,1) = 2 .* r2moy(:,2) - r2moy(:,3);             
r2tau2    = 1./4.*pi.^2.*arho.^3.*(48.*r0.*d.^2+24.*r0.^3.*dp.^2+72.*r0.*d.^2.*dp.^2+ ...
            36.*d.*arho.*r0.*dp.^3+54.*r0.*arho.^2.*dp.^2+48.*r0.^2.*d+54.*d.*arho.^2.*dp.^2+ ...
            18.*arho.*r0.^2.*dp.^3+24.*d.^3.*dp.^2+5.*arho.^3.*dp.^3+18.*arho.^3.*dp+ ...
            144.*r0.*d.*arho.*dp+72.*r0.^2.*d.*dp.^2+18.*d.^2.*arho.*dp.^3+ ...
            72.*r0.^2.*arho.*dp+24.*r0.*arho.^2+24.*d.*arho.^2+16.*r0.^3+16.*d.^3+ ...
            72.*d.^2.*arho.*dp)./vpr;
r2tau2(:,1) = 2 .* r2tau2(:,2) - r2tau2(:,3);             
r3invtau3 = 2.*pi.^2./vpr.*(8.*arho.^4.*dp.^2+2.*d.^4.*dp.^4-8.*r0.^3.*dp.^5.*arho- ...
            24.*arho.^3.*r0.*dp.^3-24.*d.*r0.^2.*dp.^5.*arho+16.*arho.^3.*r0.*dp- ...
            8.*d.^3.*dp.^5.*arho+8.*d.*r0.^3.*dp.^4-24.*d.^2.*r0.*dp.^5.*arho- ...
            6.*arho.^4+48.*d.*r0.*dp.^4.*arho.^2+2.*r0.^4.*dp.^4-24.*arho.^3.*d.*dp.^3+ ...
            16.*arho.^3.*d.*dp+12.*d.^2.*r0.^2.*dp.^4-12.*d.^2.*dp.^2.*arho.^2+ ...
            24.*r0.^2.*dp.^4.*arho.^2-24.*d.*r0.*dp.^2.*arho.^2-12.*r0.^2.*dp.^2.*arho.^2+ ...
            24.*d.^2.*dp.^4.*arho.^2+8.*d.^3.*r0.*dp.^4)./dp.^4./arho.^2./(dp-1)./ ...
            (dp+1)./(dp.^2-1).^(1./2).*qq+2.*pi.^2./vpr.*(-12.*d.^2.*dp.^2.*arho.^2.* ...
            (dp.^2-1).^(1./2)-12.*r0.^2.*dp.^2.*arho.^2.*(dp.^2-1).^(1./2)+16.*arho.^3.*d.* ...
            dp.*(dp.^2-1).^(1./2)+16.*arho.^3.*r0.*dp.*(dp.^2-1).^(1./2)-16.*arho.^3.*r0.*dp.^3.* ...
            (dp.^2-1).^(1./2)+12.*r0.^2.*dp.^4.*arho.^2.*(dp.^2-1).^(1./2)+12.*d.^2.*dp.^4.* ...
            arho.^2.*(dp.^2-1).^(1./2)-24.*d.*r0.*dp.^2.*arho.^2.*(dp.^2-1).^(1./2)+24.*d.*r0.* ...
            dp.^4.*arho.^2.*(dp.^2-1).^(1./2)-16.*arho.^3.*d.*dp.^3.*(dp.^2-1).^(1./2)+5.*arho.^4.* ...
            dp.^2.*(dp.^2-1).^(1./2)+arho.^4.*dp.^4.*(dp.^2-1).^(1./2)-6.*arho.^4.* ...
            (dp.^2-1).^(1./2))./dp.^4./arho.^2./(dp-1)./(dp+1)./(dp.^2-1).^(1./2);
r3invtau3(:,1) = 2 .* r3invtau3(:,2) - r3invtau3(:,3);             
r3invtau = 2.*pi./vpr.*(8.*d.*r0.^3.*pi+6.*r0.^2.*arho.^2.*pi+12.*d.^2.*r0.^2.*pi+ ...
            2.*d.^4.*pi+8.*d.^3.*r0.*pi+6.*d.^2.*arho.^2.*pi+3./4.*arho.^4.*pi+12.*r0.*d.*arho.^2.*pi+2.*r0.^4.*pi);
r3invtau(:,1) = 2 .* r3invtau(:,2) - r3invtau(:,3);             

        
% 10 -calcul de phi 1ere version
%  la version a partir de fdia n'est pas assez precise !
phi = 2 .* pi .* r0 .* b0 .* (r0 + d - sqrt((r0 + d) .^ 2 - a .^2 .* rho .^ 2));
rhophi  = sqrt( phi ./ pi ./ b0);

% 6 - calcul de <Jphi/R> (relier a Ip et li)
if mode == 0
	piqj = (exp(li)-1.65)./0.89;
	ind  = find(piqj <1);
	piqj(ind) = 1 * ones(1,length(ind));
	ind  = find(piqj >10);
	piqj(ind) = 10 * ones(1,length(ind));
	jmoy = r2m1 .^piqj + 1e-3 .* r2m1;
	% normalisation de j
	ipj  = (2*pi) .* a(:,1) .^ 2 .* trapz(rho(1,:),rho .*jmoy,2);
	j0   = ((ip(:,1) ./ipj) * ve); 
	jmoy = jmoy .* j0 ;
else
	if mode_rho == 0
		jmoy = tsplinet(rhoe,jmoy,rho);
	else
		jmoy = tsplinet(rhoe,jmoy,(rhophi./(rhophi(:,end)*ve)));
	end
end
jr   = jmoy .* ri;

% 7 - calcul de psi
dpsidrho      = 2 .* pi .* mu0 .* a .* cumtrapz(rho(1,:),vpr .* jr,2) ./ vpr ./ grho2r2;
dpsidrho(:,1) = zeros(size(a,1),1);
psi           = a .* cumtrapz(rho(1,:),dpsidrho,2);

% 8 - calcul du profil de Ptot
if (mode == 0)|(size(ptot,2)==1)
        % calcul de la forme
        ptot    = - a .* cumtrapz(rho(1,:),jr .* dpsidrho ./ (2*pi),2);
        ptot    = ptot - ptot(:,end) * ve;
        % normalisation
	wtot    = (3/2) .* a(:,1) .* trapz(rho(1,:),ptot.*vpr,2);
	wm      = wdia(:,1);
	ind     = find(wm<0);
	wm(ind) = zeros(1,length(ind));
	p0      = wm ./ wtot;
	ptot    = ptot .* (p0*ve);
else
	ptot    = ptot .* ee .* 1e3;    % keV*m^-3 -> Pa
	if mode_rho == 0
		ptot = tsplinet(rhoe,ptot,rho);
	else
		ptot = tsplinet(rhoe,ptot,(rhophi./(rhophi(:,end)*ve)));
	end
end

% 9 - calcul de la fontion diamagnetique 
df2drho  = mu0 .* ( - jr .* dpsidrho ./ (2*pi) - pdederive(a .* rho,ptot,0,2,2,1)) ./ r2i;
fdia2    = 2 .* a .* cumtrapz(rho(1,:),df2drho,2);
% normalisation
fdia     = sqrt(fdia2 - fdia2(:,end) * ve + r0.^2 .* b0.^2);
df2drho  = pdederive(a .* rho,fdia,0,2,2,1);

% 8 - bis :  si mode = 0 petite convergence sur la pression
if (mode == 0)  |(size(ptot,2)==1)
        % calcul de la forme
        ptot    = a .* cumtrapz(rho(1,:),- jr .* dpsidrho ./ (2*pi) -  ...
                    df2drho ./ 2. / mu0 .* r2i,2);
        ptot    = ptot - ptot(:,end) * ve;
        % normalisation
	wtot    = (3/2) .* a(:,1) .* trapz(rho(1,:),ptot.*vpr,2);
	wm      = wdia(:,1);
	ind     = find(wm<0);
	wm(ind) = zeros(1,length(ind));
	p0      = wm ./ wtot;
	ptot    = ptot .* (p0*ve);
end

% 11 - facteur de securite et bpol 
bpol   = dpsidrho ./ (2*pi) ./ (r0 + d + a .* rho .* dp ./ 2);
q      = fdia .* a .* rho ./ bpol .* r2i;
q(:,1) = 2 .* q(:,2) - q(:,3);

% 12 - fraction de pieges
epsi  = a .* rho ./ (r0 +d);
ftrap = 1 - (1 - epsi) .^ 2 ./ sqrt(1 - epsi.^2) ./ (1 + 1.46 .* sqrt(epsi)); 

if modeb2 == 1
    % 13 - <B^2>
    b2   = bpol.^2 + fdia .^2 .* r2i;

    % 14 - <1/B^2>
    theta = linspace(0,2*pi,73)';
    va    = ones(size(theta));
    b2i   = zeros(size(fdia));
    for k = 1:size(a,1)
            inter = (va * r0(k,:) + va * d(k,:) + cos(theta) * (a(k,:) .* rho(k,:))) .^ 3 ./ ...
                    ((va * ((dpsidrho(k,:) ./ (2*pi)) .^ 2)) ./ ( 1 + cos(theta) * dp(k,:)) .^ 2 + ...
                    va * fdia(k,:) .^ 2) ;
            b2i(k,:) = (2*pi) .* trapz(theta,inter .* (va * (a(k,:) .* rho(k,:))) .* ( 1 + cos(theta) * dp(k,:)),1) ./ vpr(k,:);
            b2i(k,1) = 2 * b2i(k,2) - b2i(k,3);        
    end
else
	b2=[];
	b2i=[];
end

% 15 - recalcul de phi :
phi   = a .* cumtrapz(rho(1,:),fdia .* vpr .* r2i,2) ./ (2*pi);
phi   = phi - phi(:,1) * ve;

% 16 - coordonnees de sortie et d'unites
rhog = a .* rho;
ptot = ptot ./ (ee * 1e3);

warning on
