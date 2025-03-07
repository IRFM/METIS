% interface to METIS LH model
% input:
%	  cons.temps    = time slices vector [n_time * 1]
%     cons.ip       = plasma current (A) [n_time * 1]
%     cons.plh      = LH input power (W) [n_time * 1]
%
%  	  profil.xli    = Lao coordinate (r/a) [1 * 21]  
%  	  profil.Raxe   = magnetic axis of each flux surface (m) [n_time * 21] 
%  	  profil.epsi   = inverse aspect ratio (a(x) / Raxe(x)) of each flux surface (m) [n_time * 21] 
%  	  profil.fdia   = diamagnetic function (R*B_T in T.m)   [n_time * 21] 
%  	  profil.qjli   = safety factor  [n_time * 21]
%  	  profil.nep    = electron density (m^-3) [n_time * 21]
%  	  profil.tep    = electron temperature (eV) [n_time * 21]
%  	  profil.rmx    = toroidal flux coordinate (sqrt(phi_tor/pi/B0) in m) [n_time * 21]
%  	  profil.spr    = dS/dxli surface element (m^2) [n_time * 21]
%  	  profil.vpr    = dV/dxli surface element (m^3) [n_time * 21]
%  	  profil.zeff   = effective charge [n_time * 21]
%     profil.epar   = parallel electrique field (V/m) [n_time * 21]
%  
%  	  option.gaz         = 1 -> H, 2 -> D , 3 -> DT ,4 -> He, 5 -> D+He3 and 11 -> p-B11
%  	  option.freqlh      = Lower Hybrid frequency (GHz)
%     option.etalh       = launcher directivity defined as the fraction of total LH power in the co-current peak
%  	  option.npar0       = launched parallel refractive index of LH at antenna
%  	  option.wlh         = wlh is the width of LH antenna active part (m)
%  	  option.npar_neg    = parallel refractive index of negative peak in the spectrum at the launcher; if = 0, used npar_neg = -npar0
%     option.fupshift    = factor applied to kinetic resonance position: n_par_Landau = fupshift * 6.5 / sqrt(Te); default = 1
%
% remarks: 
%  1/ parameters internally fixed
%       option.xlh     = 0.2 -> backup value for maximum of deposition if the  model do not converge
%       option.dlh     = 0.3 -> backup value for the deposition  width if the  model do not converge
%       option.lhmode  = 0   -> eta_lh is the directivity
%
%  2/ force the model to newmodel
%       option.upshiftmode = 'newmodel';
%
% output:
%
%        time          = copy of cons.temps [n_time * 1]
%        plh_tot       = LH absorbed power (W)  [n_time * 1]
%        ilh           = LH current drive (A)  [n_time * 1]
%        x             = copy of profil.xli    [1 * 21]
%        plh           = LH absorbed power density (W/m^3) [n_time * 21]
%        jlh           = LH current drive density (A/m^2) [n_time * 21]
%        efficiency    = nomalized current drive efficiency (A/W/m^2) [n_time * 1]
%
% test from metis :
%
% [time,plh_tot,ilh,x,plh,jlh,efficiency] = external_call_metis_lh_model(post.z0dinput.cons,post.profil0d,post.z0dinput.option);
% figure;
% subplot(2,2,1)
% plot(time,plh_tot,post.zerod.temps,post.zerod.plh)
% xlabel('time')
% ylabel('P_{LH} (W)');
% legend('EXT','METIS');
% subplot(2,2,2)
% plot(time,ilh,post.zerod.temps,post.zerod.ilh)
% xlabel('time')
% ylabel('I_{LH} (A)');
% legend('EXT','METIS');
% subplot(2,2,3)
% plot(time,efficiency,post.zerod.temps,post.zerod.efficiency)
% xlabel('time')
% ylabel('efficiency');
% legend('EXT','METIS');
% figure;
% subplot(1,2,1)
% zplotprof(gca,time,x,plh,'color','b');
% hold on
% zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.plh,'color','r');
% xlabel('r/a');
% ylabel('P_{LH} (W/m^3)');
% legend('EXT','METIS');
% subplot(1,2,2)
% zplotprof(gca,time,x,jlh,'color','b');
% hold on
% zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jlh,'color','r');
% xlabel('r/a');
% ylabel('J_{LH} (A/m^2)');
% legend('EXT','METIS');



function [time,plh_tot,ilh,x,plh,jlh,efficiency] = external_call_metis_lh_model(cons,profil,option)


% internal METIS parameters 
transitoire = 1;
% backup value
option.xlh  = 0.2;
option.dlh  = 0.3;
option.lhmode  = 0;
% force the model
option.upshiftmode = 'newmodel';
% no ripple 
friplh = 0;
vt = ones(size(cons.temps));
ve = ones(size(profil.xli));
x    = profil.xli;
zeff = trapz(x, profil.zeff,2);
% directivite
switch option.lhmode
case 0
   directivity = abs(option.etalh);
case 3
   directivity = abs(option.etalh);
case 4
   directivity = abs(option.etalh);
otherwise
   directivity = 0.75;
end
directivity = max(0,min(1,directivity));
xlhin  =  option.xlh * vt;
dlhin  =  option.dlh * vt;
geo.a  = profil.Raxe(:,end) .* profil.epsi(:,end);
geo.R  = profil.Raxe(:,end);
geo.b0 = profil.fdia(:,end) ./ profil.Raxe(:,end);
%

switch gaz
    case 1
        agaz    = 1 .* vt;
        zgaz    = 1 .*vt;
    case 2
        agaz    = 2 .* vt;
        zgaz    = 1 .* vt;
    case 3
        agaz   = 2.5 .* vt;
        zgaz    = 1 .* vt;
    case 5
        agaz   = 2.5 .* vt;
        zgaz    = 1.5 .* vt;
    case 11
        agaz   = 1 .* vt;
        zgaz    = 1 .* vt;
    otherwise
        agaz    = 4 .* vt;
        zgaz    = 2 .* vt;
end
agaz = agaz * vt;
zgaz = zgaz * vt;
% facteur de propagation LHCD
%qcyl          =  5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ R .* ( 1 + (a ./ R) .^ 2);
qcyl          =  5 .* geo.a .^ 2 .* geo.b0 ./ (cons.ip./1e6) ./ geo.R;
qbord = profil.qjli(:,end);
switch option.upshiftmode
case {'newmodel','newmodel + tail'}
  upshift       = option.fupshift .* vt;
otherwise
  upshift       = option.fupshift .* max(eps,qbord ./ qcyl - 1);
end
%
[xvoid,fplh,xlh,dlh,efficiency,rapnegpos] = z0lhacc2lobes(option.freqlh.*1e9,option.npar0,option.wlh,agaz,zgaz,cons.temps,profil.xli,profil.nep,profil.tep,...
				profil.qjli,profil.Raxe,profil.rmx,profil.spr,profil.vpr,geo.b0,max(1,cons.plh),xlhin,dlhin,transitoire,directivity, ...
				friplh,upshift,0,option.upshiftmode,option.npar_neg);
if isfield(profil,'zeff')
	fjlh = imag(fplh) ./ (5 + profil.zeff) .* (1 - profil.epsi .^ ((5 + profil.zeff) ./ 2 ./ (1 + profil.zeff)));
else
	fjlh = imag(fplh);
end
fplh = abs(real(fplh));
normplh = max(1,cons.plh) ./ trapz(x,profil.vpr .* fplh,2);
fplh = fplh .* (normplh*ve);
%
%  profil de densite de courant
%  cette facon de le decrire le rend compatible avec le module sans lobe negatif
%  
jlh    = fplh .* ((fplh > 0) + (rapnegpos * ve) .* (fplh <=0));
indjlh = find(any(fjlh~=0,2));
if ~isempty(indjlh)
	jlh(indjlh,:) = fjlh(indjlh,:);
end
% normalisation sur le courant LH
indbadlh = find(trapz(x,profil.spr .* jlh,2) < eps);
if ~isempty(indbadlh)		
	jlh(indbadlh,:) = ones(length(indbadlh),1) *mean(fplh .* (fplh > 0),1);
end	
% ajout de l'effet du Zeff
efficiency = efficiency ./ (5 + zeff) .* (1 - (xlh .* geo.a ./ geo.R) .^ ((5 + zeff) ./ 2 ./ (1 + zeff))); 
if isfield(profil,'zeff') && ~isempty(indjlh)
	efficiency_full    =  trapz(x,profil.spr .* fjlh ,2) ./ max(1,trapz(x,profil.spr .* fplh,2)) .* trapz(x, profil.nep,2);
	efficiency(indjlh) =  efficiency_full(indjlh);
end
vloop_lh    = 2 .* pi .*  trapz(profil.xli,profil.Raxe .* profil.epar .* fplh .* profil.spr,2) ./ ...
	                          max(1,trapz(profil.xli,fplh .* profil.spr,2));
nbar   =  trapz(x, profil.nep,2);
xlh    =  max(1,cons.plh) ./ nbar ./ geo.R;
rhot   =  8 .* nbar .* geo.R ./ xlh ./ efficiency .^ 2 .* (3 + zeff) ./ (5 + zeff) .^ 2; % ok
ilh    = (cons.plh > 1) .* (efficiency .* xlh + vloop_lh ./ rhot);
jlh    = jlh  .* ((ilh./  max(eps,trapz(x,profil.spr .* jlh,2))) * ve);
%
% profil de source de chaleur eletronique LH
%
plh        = abs(fplh);
plh_tot    = trapz(x,profil.vpr .* fplh,2);
time       = cons.temps;
	
	