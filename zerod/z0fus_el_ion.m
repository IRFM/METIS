function [ecrit_ave,pfus_el,pfus_ion,pfus_shape_el,pfus_shape_ion] = z0fus_el_ion(zerod,profil0d,option,zu1,zu2)

% resampling (for testing)
if isfield(profil0d,'temps') && (length(zerod.temps) ~= length(profil0d.temps))
    pfus_th = interp1(zerod.temps,zerod.pfus_th,profil0d.temps,'nearest','extrap');
    nDi = interp1(zerod.temps,zerod.nDm,profil0d.temps,'nearest','extrap');
    nTi = interp1(zerod.temps,zerod.nTm,profil0d.temps,'nearest','extrap');
    n1i = interp1(zerod.temps,zerod.n1m,profil0d.temps,'nearest','extrap');
else
    pfus_th = zerod.pfus_th;
    nDi = zerod.nDm;
    nTi = zerod.nTm;
    n1i = zerod.n1m;
end

% vectors
ve = ones(size(profil0d.xli)); 
% for testing
if nargin < 5
  % valeurs utiles pour les impuretees
  zu1 = (option.zimp + option.rimp .* option.zmax);
  zu2 = (option.zimp .^ 2 + option.rimp .* option.zmax .^ 2);
elseif length(zu1) > 1
  zu1 = zu1 * ve;
  zu2 = zu2 * ve;
end

% improve precision
[A_el,Z_el] = chargemasse;
dd   = abs(Z_el - option.zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* option.zimp;
end
dd   = abs(Z_el - option.zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* option.zmax;
end



% mass factor
n1  = profil0d.n1p .* ((n1i./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .*  ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
nD  = profil0d.n1p .* ((nDi./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .*  ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
nT  = profil0d.n1p .* ((nTi./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
switch option.gaz
  case 11
      nB = nT;
      nT = zeros(size(nB));
      % probleme d'arrondis numeriques
      n1  = max(n1,nD);
      nH  = n1 - nD;
  otherwise
      nB = zeros(size(nT));
      % probleme d'arrondis numeriques
      n1  = max(n1,nD+nT);
      nH  = n1 - nD -nT;
end
%
[zw,zw2] = z0wavez(profil0d.tep);
[zsn,zsn2] = z0snavez(profil0d.tep);

if isfield(profil0d,'nwp')
    if option.Sn_fraction > 0
        fact = nD ./2 + nT ./ 3 + nH  + option.zimp .^ 2 .* profil0d.nzp ./ aimp + ...
               option.rimp .* option.zmax .^ 2 .* profil0d.nzp ./ amax +  ...
               25  .* nB  ./ 11.0093054 + (1 - option.Sn_fraction) .* zw2 .* profil0d.nwp ./ 183.84 + ...
               option.Sn_fraction  .* zsn2 .* profil0d.nwp ./ 118.71;
    else        
        fact = nD ./2 + nT ./ 3 + nH +  option.zimp .^ 2 .* profil0d.nzp ./ aimp + ...
               option.rimp .* option.zmax .^ 2 .* profil0d.nzp ./ amax +  ...
               25  .* nB  ./ 11.0093054 + zw2 .* profil0d.nwp ./ 183.84;
    end
    switch option.gaz
        case 5
             fact = fact + (4/3) .* profil0d.nhep + option.frhe0 .* profil0d.nep;
       otherwise
            fact = fact + profil0d.nhep;
    end
else
    fact = nD ./2 + nT ./ 3 + (n1 - nT - nD) + profil0d.nhep + zu2 .* profil0d.nzp ./ aimp;   
end
%fact = nD ./2 + nT ./ 3 + (n1 - nT - nD) + profil0d.nhep + zu2 .* profil0d.nzp ./ aimp;
fact = (fact ./ profil0d.nep);
% temps de thermalisation (plutot central, calibree avec spot)
%lnldei     = 15.2 - 0.5 .* log(profil0d.nep./1e20) + log(profil0d.tep ./1e3);
%taus       = 6.27e8 .* profil0d.tep  .^ (3/2) ./ (profil0d.nep ./ 1e6) ./ lnldei;
%ecrit = 14.8 .* profil0d.tep  .* (4 .* profil0d.zeff) .^ (2/3);
% 4 is for He4 nucleus
ecrit  = 14.8 .* profil0d.tep .* (4 .^ (3/2)  .* fact) .^ (2/3);

% fusion
fpfus    = profil0d.salf;
fpfus    = fpfus .* ((pfus_th./  max(1,trapz(profil0d.xli,profil0d.vpr .* fpfus,2))) * ve);
% deposition
switch option.gaz
    case 11
        % from p-B11 reactions
        % the spectrum of alpha particles is broad
        % references: Plasma 2023, 6(3), 379-392; https://doi.org/10.3390/plasma6030026
        % V. F. Dmitriev, Physics of Atomic
        % Nuclei, 2009, Vol. 72, No. 7, pp. 11651167, DOI: 10.1134/S1063778809070084
        % approximated with one alpha at 4 MeV and two at 2.3Mev
        pfus_shape_ion = min(1 - eps,(option.alpha_channeling + (2./3) .* zfract0(ecrit,2.3e6) + (1/3) .* zfract0(ecrit,4e6))) .* fpfus;
   otherwise
        % from DT reactions
        pfus_shape_ion = min(1 - eps,(option.alpha_channeling + zfract0(ecrit,3.56e6))) .* fpfus;
end
pfus_shape_el  = max(0,fpfus - pfus_shape_ion);
% zerod data 
pfus_el  = trapz(profil0d.xli,profil0d.vpr .* pfus_shape_el,2);
pfus_ion = trapz(profil0d.xli,profil0d.vpr .* pfus_shape_ion,2);
ecrit_ave = trapz(profil0d.xli,profil0d.vpr .* (pfus_shape_el + pfus_shape_ion) .* ecrit,2) ./ ...
            trapz(profil0d.xli,profil0d.vpr .* max(1,pfus_shape_el + pfus_shape_ion) ,2);

% if not for testing
if nargin >= 5
  return
end
% control figure

% fraction de la puissance deposee sur les ions
% reference : D. Sigmar and J. Joyce, Nuclear fusion, 11(3), 1971.
% formule etablie pour DT, (fit a partir d'une figure)
%  0.2 < te < 80 keV 
t=[0,10,20,30,40,50,60,70,80,100].*1e3;
f=[0.01,0.146,0.269,0.377,0.473,0.55,0.615,0.665,0.698,0.72];
fion = interp1(t,f,profil0d.tep,'pchip');
pfus_shape_ion_ctrl = min(1 - eps,(option.alpha_channeling + fion)) .* fpfus;
pfus_ion_ctrl = trapz(profil0d.xli,profil0d.vpr .* pfus_shape_ion,2);

if ~isfield(profil0d,'temps')
  if isfield(zerod,'temps')
	profil0d.temps = zerod.temps;
  else
	profil0d.temps = 1:length(pfus_el);
	zerod.temps    = linpace(min(profil0d.temps),max(profil0d.temps),length(zerod.pion_fus));
  end
end

figure(21)
plot(profil0d.temps,pfus_el ./ 1e6,'r',profil0d.temps,pfus_ion ./ 1e6,'b', ...
     zerod.temps,zerod.pel_fus ./ 1e6,'m',zerod.temps,zerod.pion_fus ./ 1e6,'c', ...
     zerod.temps,zerod.pfus_th ./ 1e6,'g',profil0d.temps,(pfus_el + pfus_ion)  ./ 1e6,'--k', ...
     profil0d.temps,pfus_ion_ctrl./1e6,'.k');
xlabel('time (s)');
ylabel('MW')
legend('P_{el} (profile)','P_{ion} (profile)','P_{el} (0d)','P_{ion} (0d)')
edition2
drawnow
