% function [] = build_local_neo_input_from_METIS_default(flnm,rova,timec)
function [neoin_mod ynorm yspec] = build_local_neo_inputs(option,z0dinput,profil,neoin_mod,rova,timec,write,calc_DV)

if ~exist('write')
  write = 0;
end

signJ = 1; % Counterclockwise for WEST when viewed from above -> +1
signB = 1; % Same direction in WEST

[aa itime] = min(abs(profil.temps-timec));
[aa itimez] = min(abs(z0dinput.cons.temps-timec));

Raxe = profil.Raxe(itime,1);
BT = z0dinput.geo.b0(itime);

% Phys cst
eps0 = 8.8542e-12;
mp = 1.6726e-27;
me = 9.1096e-31;
qe = 1.602e-19;

% Fraction of hydrogen minority
nch = option.cmin;
  
switch option.gaz
case 1
   zi = 1;
   mi = 1*mp;
case 2
   zi = 1;
   mi = 2*mp;
case 3
   zi = 1;
   mi = (2 .* (1 - z0dinput.cons.iso(itimez)) +  3 .* z0dinput.cons.iso(itimez))*mp;
case 4
   % If helium be careful, completely blind on Ti on WEST not even neutrons to measure
   zi = 2;
   mi = 4*mp;
end


rmina = (profil.epsi(itime,end).*profil.Raxe(itime,end));
rmin = profil.epsi(itime,:).*profil.Raxe(itime,:);
rmaj = profil.Raxe(itime,:);
rmin_ov_a = rmin./rmina;

kk = profil.kx(itime,:);
dd = profil.dx(itime,:);

te = profil.tep(itime,:);
ti = profil.tip(itime,:);
ne = profil.nep(itime,:);

n1 = profil.nzp(itime,:);
n2 = profil.nzp(itime,:)*option.rimp;
n3 = profil.nwp(itime,:);
  
z1 = option.zimp;
z2 = option.zmax;
if ~isfield(option,'Sn_fraction')
   Sn_fraction = 0;
 else
   Sn_fraction = option.Sn_fraction;
 end

z3 = round((1 - Sn_fraction) .* z0wavez(te) + Sn_fraction .* z0snavez(te));
ni = (ne-z1.*n1-z2.*n2-z3.*n3)./(1+nch)./zi;
n4 = nch.*ni;



qmag = profil.qjli(itime,:);

% Loads neo inputs (change template to account for x species etc)

neoin_mod.RMIN_OVER_A=rova;
rok = rova.*rmina;
neoin_mod.RMAJ_OVER_A=interpos(rmin,rmaj,rok,1e-6,[1 0], [0 0]);
neoin_mod.RMAJ_OVER_A=neoin_mod.RMAJ_OVER_A./rmina;

% Interp profiles on rmin at rova

teneo= interpos(rmin,te,rok,1e-6,[1 0], [0 0]);
alte = interpos(rmin,-gradient(te)./gradient(rmin)./te.*rmin,rok,1e-6,[1 0], [0 0]);

tineo= interpos(rmin,ti,rok,1e-6,[1 0], [0 0]);
alti = interpos(rmin,-gradient(ti)./gradient(rmin)./ti.*rmin,rok,1e-6,[1 0], [0 0]);

neneo= interpos(rmin,ne,rok,1e-6,[1 0], [0 0]);
alne = interpos(rmin,-gradient(ne)./gradient(rmin)./ne.*rmin,rok,1e-6,[1 0], [0 0]);

nineo= interpos(rmin,ni,rok,1e-6,[1 0], [0 0]);
%alni = -dy./nineo.*rmina;

n1neo= interpos(rmin,n1,rok,1e-6,[1 0], [0 0]);
aln1 = interpos(rmin,-gradient(n1)./gradient(rmin)./n1.*rmin,rok,1e-6,[1 0], [0 0]);

n2neo= interpos(rmin,n2,rok,1e-6,[1 0], [0 0]);
aln2 = interpos(rmin,-gradient(n2)./gradient(rmin)./n2.*rmin,rok,1e-6,[1 0], [0 0]);

n3neo= interpos(rmin,n3,rok,1e-6,[1 0], [0 0]);
aln3 = interpos(rmin,-gradient(n3)./gradient(rmin)./n3.*rmin,rok,1e-6,[1 0], [0 0]);

% If flat profile used, force a nonzero gradient to compute D and V
if max(aln3)<1e-3
  aln3 = teneo./teneo.*0.1;
end

n4neo= interpos(rmin,n4,rok,1e-6,[1 0], [0 0]);
aln4 = interpos(rmin,-gradient(n4)./gradient(rmin)./n4.*rmin,rok,1e-6,[1 0], [0 0]);

[z3neo]= interpos(rmin,z3,rok,1e-6,[1 0], [0 0]);

if calc_DV >= 1 && calc_DV <= 5
   switch calc_DV
       case 1
           aln1 = aln1 * 1.05;
       case 2
           aln2 = aln2 * 1.05;
       case 3
           aln3 = aln3 * 1.05;
       case 4
           aln4 = aln4 * 1.05;
       case 5
           alne = alne * 1.05;
   end
elseif calc_DV == 6
   aln1 = aln1 * 1.05;
   aln2 = aln2 * 1.05;
   aln3 = aln3 * 1.05;
   aln4 = aln4 * 1.05;
   alne = alne * 1.05;
end
alni = (neneo.*alne-aln1.*n1neo.*z1-aln2.*n2neo.*z2-aln3.*n3neo.*z3neo)./nineo-aln4.*nch;

kneo= interpos(rmin,kk,rok,1e-6,[1 0], [0 0]);
dkneo = interpos(rmin,gradient(kk)./gradient(rmin),rok,1e-6,[1 0], [0 0]);

tmp= interpos(rmin,rmaj,rok,1e-6,[1 0], [0 0]);
shift = interpos(rmin,gradient(rmaj)./gradient(rmin),rok,1e-6,[1 0], [0 0]);

ddneo= interpos(rmin,dd,rok,1e-6,[1 0], [0 0]);
dddneo = interpos(rmin,gradient(dd)./gradient(rmin),rok,1e-6,[1 0], [0 0]);

coulog = 23-log(zi.^2*sqrt(2.*(ni./1e6)*zi.^2))+1.5*log(ti);
tauii = (3*eps0^2*mi^2*(2*pi*qe/mi.*ti).^1.5)./(qe^4*zi^4.*ni.*coulog);

qneo= interpos(rmin,qmag,rok,1e-6,[1 0], [0 0]);
smag = interpos(rmin,gradient(qmag)./gradient(rmin)./qmag.*rmin,rok,1e-6,[1 0], [0 0]);

tauneo = interpos(rmin,tauii,rok,1e-6,[0 0], [0 0]);

tnorm = tineo;
vnorm = sqrt(tnorm*qe./mp/2);
nnorm = nineo;
mnorm = mi;
anorm = rmina;
rnorm = neoin_mod.RMAJ_OVER_A*rmina;

%Charges
neoin_mod.Z_1 = 1;
neoin_mod.Z_2 = z1;
neoin_mod.Z_3 = z2;
neoin_mod.Z_4 = z3neo;
neoin_mod.Z_5 = 1;
neoin_mod.Z_6 = -1;

%Masses
neoin_mod.MASS_1 = 1;
neoin_mod.MASS_2 = z1 * 2 * mp / mi;
neoin_mod.MASS_3 = z2 * 2 * mp / mi;
neoin_mod.MASS_4 = 92 * 2 * mp / mi;
neoin_mod.MASS_5 = 1 * mp / mi;
neoin_mod.MASS_6 = 0.00027245;

%Densities
neoin_mod.DENS_1 = 1.;
neoin_mod.DENS_2 = max(n1neo ./ nnorm,1e-8);
neoin_mod.DENS_3 = max(n2neo ./ nnorm,1e-8);
neoin_mod.DENS_4 = max(n3neo ./ nnorm,1e-8);
neoin_mod.DENS_5 = max(n4neo ./ nnorm,1e-8);
neoin_mod.DENS_6 = max(neneo ./ nnorm,1e-8);

%Density gradients
neoin_mod.DLNNDR_1 = alni;
neoin_mod.DLNNDR_2 = aln1;
neoin_mod.DLNNDR_3 = aln2;
neoin_mod.DLNNDR_4 = aln3;
neoin_mod.DLNNDR_5 = aln4;
neoin_mod.DLNNDR_6 = alne;

%Temperature gradients
neoin_mod.DLNTDR_1 = alti;
neoin_mod.DLNTDR_2 = alti;
neoin_mod.DLNTDR_3 = alti;
neoin_mod.DLNTDR_4 = alti;
neoin_mod.DLNTDR_5 = alti;
neoin_mod.DLNTDR_6 = alte;

%Temperature
neoin_mod.TEMP_1 = 1;
neoin_mod.TEMP_2 = 1;
neoin_mod.TEMP_3 = 1;
neoin_mod.TEMP_4 = 1;
neoin_mod.TEMP_5 = 1;
neoin_mod.TEMP_6 = teneo ./ tineo;


% Diffusion de l'ion
neoin_mod.NU_1 = rmina ./ tauneo ./ vnorm; 


% Geometry
neoin_mod.KAPPA=  kneo;
neoin_mod.S_KAPPA=  dkneo./kneo.*rok;
neoin_mod.DELTA=  ddneo;
neoin_mod.S_DELTA=  dddneo.*rok;
neoin_mod.SHIFT= shift;
neoin_mod.Q=  qneo;
neoin_mod.SHEAR=  smag;

% General
neoin_mod.OMEGA_ROT=  0.;
neoin_mod.OMEGA_ROT_DERIV= 0.;

% Here we use Bunit = BT for simple circular geometry etc
% To be changed consistently with the toroidal magnetic flux

rhotmp=sqrt(abs(profil.phi(itime,:))/pi./BT);  
yy= interpos(rmin,rhotmp,rok,1e-6,[0 0], [0 0]);
drhodr = interpos(rmin,gradient(rhotmp)./gradient(rmin),rok,1e-6,[1 0], [0 0]);
rhott = yy;
bunit = 1./rok.*BT.*rhott.*drhodr;

neoin_mod.RHO_STAR=  vnorm*mnorm./zi./qe./bunit./rmina;

ynorm.rhostar = vnorm*mnorm./zi./qe./bunit./rmina;
ynorm.bunit = bunit;
ynorm.tnorm = tnorm;
ynorm.vnorm = vnorm;
ynorm.nnorm = nnorm;
ynorm.mnorm = mnorm;
ynorm.anorm = anorm;
ynorm.rnorm = rnorm;
ynorm.rhotn = interpos(rmin,sqrt(profil.phi(itime,:)./profil.phi(itime,end)),rok,1e-6,[0 0], [0 0]);

yspec.n1ovni = n1neo./nnorm;
yspec.n2ovni = n2neo./nnorm;
yspec.n3ovni = n3neo./nnorm;
yspec.n4ovni = n4neo./nnorm;

yspec.alnz1 = aln1;
yspec.alnz2 = aln2;
yspec.alnz3 = aln3;
yspec.alnz4 = aln4;
yspec.alti = alti;

end
