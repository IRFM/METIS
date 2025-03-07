% na  = densite au temps t
% bord = condition au limite au bord 0 -> densite donnee, 1 -> flux donne
% cbord = valeur au bord de la densite ou du flux
% nap = densite au temps t+dt
% rhomax   = equi.rhomax
% drhomaxdt   = equi.drhomaxdt
% grho2    = equi.grho2
% vpr      = equi.vpr
% dvprdt   = equi.dvprdt
% dn       = coeficient de diffusion au temps t
% dnp      = coefficient de diffusion au temps t+dt
% vn       = vitesse de convection au temps t
% vnp      = vitesse de convection au temps t+dt
% sn       = source de matiere au temps t
% snp      = source de matiere au temps t+dt
% x        = coordonnees normalisee
% dt       = intervalle de temps d'integration
% dtmax    = intervalle de temps maximum pour la convergence
% cn       = mode implicite ou C-N
function nap=zdiffusedensite(na,bord,cbord,rhomax,drhomaxdt,grho2,vpr,dvprdt,dn,dnp,vn,vnp,sn,snp,x,dt,dtmax,cn)


if nargin < 18
	cn = 0.5; % implicite pure
end

% constantes
nbrho = length(na);
dx    = mean(diff(x));

% calcul des derivees spatiales premiere
% la derivee premiere est nulle en rho = 0 (=> en x=0)
grho2d1          = pdederive(x,grho2,0,2,2,1);
vprd1            = pdederive(x,vpr,0,2,2,1);
vnd1             = pdederive(x,vn,0,2,2,1);
vnpd1            = pdederive(x,vnp,0,2,2,1);

% coef de l'equation sur Ne, la reliant a Ne (derivee 2)
warning off
Ann  =  grho2.*dn./rhomax.^2;
Annp =  grho2.*dnp./rhomax.^2;

% coef de l'equation sur Ne, la reliant a Ne (derivee 1)
Bnn  =  1./rhomax.*grho2.*vn+1./vpr./rhomax.^2.*dn.*vprd1.*grho2 ...
       +1./rhomax.^2.*dn.*grho2d1 ...
       +1./rhomax.*x.*drhomaxdt;
Bnnp =  1./rhomax.*grho2.*vnp+1./vpr./rhomax.^2.*dnp.*vprd1.*grho2 ...
       +1./rhomax.^2.*dnp.*grho2d1 ...
       +1./rhomax.*x.*drhomaxdt;

% coef de l'equation sur Ne, la reliant a Ne
Cnn  =  (1./vpr./rhomax.*x.*drhomaxdt+1./vpr./rhomax.*grho2.*vn).*vprd1-1./vpr.*dvprdt ...
        +1./rhomax.*grho2d1.*vn+1./rhomax.*grho2.*vnd1;
Cnnp  =  (1./vpr./rhomax.*x.*drhomaxdt+1./vpr./rhomax.*grho2.*vnp).*vprd1-1./vpr.*dvprdt ...
        +1./rhomax.*grho2d1.*vnp+1./rhomax.*grho2.*vnpd1;


% source pour Ni
Dnn     =  sn ;
Dnnp    =  snp ;
warning on 

% remplissage des matrices A, B et C
A      = zreshape(Ann,nbrho,1,1);
B      = zreshape(Bnn,nbrho,1,1);
C      = zreshape(Cnn,nbrho,1,1);
D      = zreshape(Dnn,nbrho,1,1);
AP     = zreshape(Annp,nbrho,1,1);
BP     = zreshape(Bnnp,nbrho,1,1);
CP     = zreshape(Cnnp,nbrho,1,1);
DP     = zreshape(Dnnp,nbrho,1,1);

% condition au centre
V0 = zeros(2,1);
T0 = zeros(2,1,3);
% deriveree nulle au centre
T0(:,1,2) = ones(2,1,1);

% condition au bord
V1 = zeros(2,1);
T1 = zeros(2,1,3);
% au temps t ->  densite donnee
V1(1,1)    = na(end);
T1(1,1,1)  = 1;

% au temps t+ dt -> selon mode 
if bord == 0
	V1(2,1)    = cbord;
	T1(2,1,1)  = 1;
else
	V1(2,1)    = cbord;
	T1(2,1,1) = - vnp(end) .* equi.grho2(end);
	t1(2,1,2) = - dnp(end) .* equi.grho2(end) ./ equi.rhomax;
end


% diffusion pour un pas de temps
if dt <= dtmax
   nap=pde1dsolver(A,B,C,D,AP,BP,CP,DP,na',0 .* na',V0,T0,V1,T1,0,cn,dx,dt)';
else
   dtc = 0;
   nap = na;
   while ((dt -dtc) > 0)
      dti = min(dtmax,dt - dtc + eps);
      dtc = dtc + dti;
      nap = pde1dsolver(A,B,C,D,AP,BP,CP,DP,nap',0 .* na',V0,T0,V1,T1,0,cn,dx,dti)';
      V1(1,1)    = nap(end);
   end
end

function s=zreshape(e,nx,ny,nz)

% continuite au centre pour les NaN
if ~isfinite(e(1))
	e(1) = (61/46) .* e(2) - (9/23) .* e(3) + (3/46) .* e(4);
end

% suppression des NaN
ind = find(~ isfinite(e));
if ~isempty(ind)
	e(ind) = zeros(1,length(ind));
end

% mise en forme du vecteur
s = reshape(e,nx,ny,nz);
