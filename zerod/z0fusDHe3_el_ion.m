function [ecrit_alpha,pfus_el,pfus_ion,pfus_shape_el,pfus_shape_ion] = z0fusDHe3_el_ion(zerod,profil0d,option,zu1,zu2)

% resampling (for testing)
if isfield(profil0d,'temps') && (length(zerod.temps) ~= length(profil0d.temps))
    pfus_th = interp1(zerod.temps,zerod.pfus_th,profil0d.temps,'nearest','extrap');
    n1i = interp1(zerod.temps,zerod.n1m,profil0d.temps,'nearest','extrap');
    nDi = interp1(zerod.temps,zerod.nDm,profil0d.temps,'nearest','extrap');
    nTi = interp1(zerod.temps,zerod.nTm,profil0d.temps,'nearest','extrap');
    nhe4i = interp1(zerod.temps,zerod.nhe4m,profil0d.temps,'nearest','extrap');
    % He4 simply factor on nep
else
    pfus_th = zerod.pfus_th;
    nDi = zerod.nDm;
    nTi = zerod.nTm;
    n1i = zerod.n1m;
    nhe4i = zerod.nhe4m;
    
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

% mass factor
n1  = profil0d.n1p .* ((n1i./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .*  ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
nD  = profil0d.n1p .* ((nDi./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .*  ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
nT  = profil0d.n1p .* ((nTi./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
nhe4  = profil0d.n1p .* ((nhe4i./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* ...
      trapz(profil0d.xli,profil0d.vpr,2)) * ve);
% probl�me d'arrondis numeriques
n1  = max(n1,nD+nT);
aimp = ceil(zu1.* (7/3));
fact = nD .*2 + nT .* 3 + (n1 - nT - nD) + profil0d.nhep.*3 + nhe4.*4 + profil0d.nzp .* aimp;
fact = (fact ./ profil0d.nip);
% temps de thermalisation (plutot central, calibree avec spot)
%lnldei     = 15.2 - 0.5 .* log(profil0d.nep./1e20) + log(profil0d.tep ./1e3);
%taus       = 6.27e8 .* profil0d.tep  .^ (3/2) ./ (profil0d.nep ./ 1e6) ./ lnldei;
%ecrit = 14.8 .* profil0d.tep  .* (4 .* profil0d.zeff) .^ (2/3);
ecrit.he4_DHe3 = 14.8 .* profil0d.tep .* (4 .^ (3/2)  .* fact) .^ (2/3);
ecrit.p_DHe3   = 14.8 .* profil0d.tep .* (1 .^ (3/2)  .* fact) .^ (2/3);
ecrit.he3_DDn  = 14.8 .* profil0d.tep .* (3 .^ (3/2)  .* fact) .^ (2/3);
ecrit.t_DDp    = 14.8 .* profil0d.tep .* (3 .^ (3/2)  .* fact) .^ (2/3);
ecrit.p_DDp    = 14.8 .* profil0d.tep .* (1 .^ (3/2)  .* fact) .^ (2/3);
ecrit.he4_DT   = 14.8 .* profil0d.tep .* (4 .^ (3/2)  .* fact) .^ (2/3);

% fusion
% fpfus    = profil0d.salf;
% fpfus    = fpfus .* ((pfus_th./  max(1,trapz(profil0d.xli,profil0d.vpr .* fpfus,2))) * ve);

fpfus  =  profil0d.palf0_tot;
fpfus   = fpfus .* ((pfus_th./  max(1,trapz(profil0d.xli,profil0d.vpr .* fpfus,2))) * ve);

% energy of fast ion produced by each reaction
dhe3.he4   = 3.71e6; % (eV) 
dhe3.p     = 14.64e6; % (eV)
dd_n.he3   = 0.82e6; % (eV)
dd_p.t     = 1.01e6; % (eV)
dd_p.p     = 3.02e6; % (ev)
dt.he4     = 3.56e6; % (eV)

% deposition
% pfus_shape_ion = min(1 - eps,(option.alpha_channeling + zfract0(ecrit,3.56e6))) .* fpfus;
% pfus_shape_el  = max(0,fpfus - pfus_shape_ion);

% deposition of each fast ion should be calculated separately because the shape
% of each fast ion distribution is different due to their cross-section is different at
% different temperature
pfus_shape_ion  =   min(1 - eps,(option.alpha_channeling + zfract0(ecrit.he4_DHe3, dhe3.he4))) .* profil0d.palf0_he4_DHe3 ...
                  + min(1 - eps,(option.alpha_channeling + zfract0(ecrit.p_DHe3  , dhe3.p  ))) .* profil0d.palf0_p_DHe3   ...
                  + min(1 - eps,(option.alpha_channeling + zfract0(ecrit.he3_DDn , dd_n.he3))) .* profil0d.palf0_he3_DDn  ...
                  + min(1 - eps,(option.alpha_channeling + zfract0(ecrit.t_DDp   , dd_p.t  ))) .* profil0d.palf0_t_DDp    ...
                  + min(1 - eps,(option.alpha_channeling + zfract0(ecrit.p_DDp   , dd_p.p  ))) .* profil0d.palf0_p_DDp    ... 
                  + min(1 - eps,(option.alpha_channeling + zfract0(ecrit.he4_DT  , dt.he4  ))) .* profil0d.palf0_he4_DT   ;

pfus_shape_el  = max(0,fpfus - pfus_shape_ion);

% zerod data 
pfus_el  = trapz(profil0d.xli,profil0d.vpr .* pfus_shape_el,2);
pfus_ion = trapz(profil0d.xli,profil0d.vpr .* pfus_shape_ion,2);
ecrit_alpha.he4_DHe3 = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_he4_DHe3 .* ecrit.he4_DHe3,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_he4_DHe3) ,2);
ecrit_alpha.p_DHe3   = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_p_DHe3 .* ecrit.p_DHe3,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_p_DHe3) ,2);
ecrit_alpha.he3_DDn  = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_he3_DDn .* ecrit.he3_DDn,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_he3_DDn) ,2);
ecrit_alpha.t_DDp    = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_t_DDp .* ecrit.t_DDp,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_t_DDp) ,2);
ecrit_alpha.p_DDp    = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_p_DDp .* ecrit.p_DDp,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_p_DDp) ,2);
ecrit_alpha.he4_DT   = trapz(profil0d.xli,profil0d.vpr .* profil0d.palf0_he4_DT .* ecrit.he4_DT,2) ./ ...
                       trapz(profil0d.xli,profil0d.vpr .* max(1,profil0d.palf0_he4_DT) ,2);

% ecrit_alpha.alpha   = ecrit_alpha.he4_DHe3 + ecrit_alpha.he4_DT;
% ecrit_alpha.proton  = ecrit_alpha.p_DHe3   + ecrit_alpha.p_DDp;
% ecrit_alpha.helium3 = ecrit_alpha.he3_DDn;
% ecrit_alpha.tritium = ecrit_alpha.t_DDp;


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