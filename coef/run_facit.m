function [yfacit_out] = run_facit(facitparams,option, z0dinput, profil,time_in,iplot,Aimp,Zimp,Nimp,fac);

if ~exist('iplot')
  iplot = 0;
end

path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')));
path_to_METIS_FACIT = [path_to_METIS_coef '/FACIT_TOOLS'];
addpath(path_to_METIS_FACIT);

eps0=8.8542e-12;
mu0=4*pi*1e-7;
me=9.1096e-31;
mp=1.6726e-27;
qe=1.602e-19;

ionregime = facitparams.ION_REGIME;
optasym = facitparams.OPT_ASYM;
labgeo = facitparams.LAB_GEO;

switch option.gaz
case 1
   zi = 1;
   mi = 1*mp;
case 2   %for WEST (D)
   zi = 1;
   mi = 2*mp;
case 3
   zi = 1;
   mi = (2 .* (1 - interp1(z0dinput.cons.temps,z0dinput.cons.iso,time_in)) +  3 .* interp1(z0dinput.cons.temps,z0dinput.cons.iso,time_in))*mp;
case 4
   zi = 2;
   mi = 4*mp;
end

irok = find(profil.xli>0);
irok = 1:length(profil.xli);
RHON = profil.xli(irok);
nx = length(RHON);
nxli = length(profil.xli);

% 0D quantities
amin = interp1(z0dinput.cons.temps,z0dinput.geo.a,time_in);
R0 =  interp1(z0dinput.cons.temps,z0dinput.geo.R,time_in);
B0 =  interp1(z0dinput.cons.temps,z0dinput.geo.b0,time_in);
ASPCT = amin./R0;

% 1D quantities
TE = interp1(profil.temps,profil.tep(:,irok),time_in);
TI = interp1(profil.temps,profil.tip(:,irok),time_in);
NE = interp1(profil.temps,profil.nep(:,irok),time_in);


if ~isfield(option,'Sn_fraction')
  Sn_fraction = 0;
else
  Sn_fraction = option.Sn_fraction;
end

if facitparams.MULTIS == 2 %z1
  % for z1
  NI = interp1(profil.temps,profil.nzp(:,irok),time_in).*fac;
  MI = 2*option.zimp.*mp;
  ZI = option.zimp;

 elseif  facitparams.MULTIS == 3 %z2
  % for z2
  NI = interp1(profil.temps,profil.nzp(:,irok).*option.rimp,time_in).*fac;
  MI = 2*option.zmax.*mp;
  ZI = option.zmax;

%elseif  facitparams.MULTIS == 4 %W
% for z3
%  NI = interp1(profil.temps,profil.nwp(:,irok),time_in).*fac;
%  MI = 184.*mp;
%  ZI = round((1 - Sn_fraction) .* z0wavez(TE) + Sn_fraction .* z0snavez(TE));
% zi=30;

elseif  facitparams.MULTIS == 4   %H
  % for H
  NI = interp1(profil.temps,profil.n1p(:,irok).*option.cmin,time_in).*fac;
  MI = mp;
  ZI = 1;

else   %main ions
  % Does it work if Helium chosen (n1p that is)?  
  NI = interp1(profil.temps,profil.n1p(:,irok),time_in);
  MI = mi;
  ZI = zi;
end

Ai = MI./mp;
Zi = ZI;

QPSI = interp1(profil.temps,profil.qjli(:,irok),time_in);
EPSILON = ASPCT*RHON;
psi = interp1(profil.temps,profil.psi(:,irok),time_in);
dpsidx = gradient(psi)./gradient(RHON);
F = R0*B0*ones(size(RHON));
% Zimp = round(z0wavez(TE));



%dpsidx = amin^2*B0*RHON./QPSI;

% Tungsten density
% Nimp = interp1(profil.temps,profil.nwp(:,irok),time_in);
% Aimp=184;

Machi = 0.*RHON;

% Asymmetries
AsymPhi = zeros(nx,2); 
AsymN = zeros(nx,2); 

% Fill in structure
in.AsymN = AsymN;
in.AsymPhi= AsymPhi;
in.x = RHON;
in.dpsidx = dpsidx;
in.FV = F;
in.Te  = TE;
in.Ti =	TI;
in.gradTi = gradient(TI,RHON);
in.Ne = NE;
in.qmag = QPSI;
in.Zi = Zi;
in.Ai = Ai;
in.Machi = Machi;
in.invaspct = amin/R0 ;
in.B0 = B0;
in.R0 = R0;
in.ionregime = ionregime ;
in.details = 0;
in.geom = labgeo;
in.pol_asym = optasym;
in.Ni = NI;
% in.jacob = jacob;
in.gradNi = gradient(NI,RHON);


in.Aimp = Aimp;
in.Zimp = Zimp;
in.Nimp = Nimp;
gradNimp_in = gradient(Nimp,RHON);
in.gradNimp = gradNimp_in;
in.grad_ln_Nimp = gradNimp_in./Nimp;

% runs FACIT in the current path and creates inputs/outputs

currentpath = pwd;
if isempty(getappdata(0,'local'))
  FACITruns = [path_to_METIS_FACIT '/runs_facit'];
else
  FACITruns = [getappdata(0,'local') '/runs_facit'];
end  
unix(sprintf(['mkdir -p ' FACITruns]));
cd(FACITruns);
filename = [pwd '/FACIT_IO.m'];
if exist(filename, 'file') ~= 2
	eval(['!cp ',path_to_METIS_FACIT,'/runs/FACIT_IO.m .'])
end

out = FACIT_IO(in,path_to_METIS_FACIT);
cd(currentpath);

% Compute tungsten flux (not correct if poloidal assymetries but too lazy to do it now)
out.fluxImp = - out.Da_neo.*in.gradNimp'./amin+in.Nimp'.*out.Va_neo;
yfacit_out = out;

if ~iplot
  return
end

figure(5001)
nlig = 3;
ncol = 3;

ip=0;

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,TE,RHON,TI,'--')
grid
title('Te & Ti')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,NE,RHON,NI)
grid
title('Ne, Ni')
legend('Ne','Ni')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,out.Va_neo)
grid
title('Vra neo')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,out.Va_cl)
grid
title('Vra cl')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
semilogy(RHON,out.Da_neo)
grid
title('Da neo')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
semilogy(RHON,out.Da_cl)
grid
title('Da cl')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,out.dmin)
grid
title('dmin')
xlabel('r/a')

ip=ip+1;
subplot(nlig,ncol,ip)
plot(RHON,out.dmaj)
grid
title('dmaj')
xlabel('r/a')



