function xiioxie = z0xiioxie(A,profli,tep,tip,stiff)


% partie gradient critique
lmin   = profli.rmx(:,end) ./ 21;
lmax   = 2 .* pi .* profli.qjli(:,end) .* profli.Raxe(:,end);
rolne  = min(1./lmin,max(1./lmax,2 .* abs(profli.nep(:,5) - profli.nep(:,end-5)) ./ max(profli.nep(:,11),1e13) ./ profli.rmx(:,11))) .* profli.Raxe(:,11);
rolni  = min(1./lmin,max(1./lmax,2 .* abs(profli.nip(:,5) - profli.nip(:,end-5)) ./ max(profli.nip(:,11),1e13) ./ profli.rmx(:,11))) .* profli.Raxe(:,11);
rolti  = min(1./lmin,max(1./lmax,2 .* abs(tip(:,5) - tip(:,end-5)) ./ max(tip(:,11),30) ./  profli.rmx(:,11))) .*  profli.Raxe(:,11);
rolte  = min(1./lmin,max(1./lmax,2 .* abs(tep(:,5) - tep(:,end-5)) ./ max(tep(:,11),30) ./  profli.rmx(:,11))) .*  profli.Raxe(:,11);
Kt     = profli.ftrap(:,11) ./ (1 - profli.ftrap(:,11));

% ref : Casaiti 2003
% separation faible gradient de densite et fort gradient
% utilisation de la limite donnees par P. Mantica
% dependnace en shear données par C. Bourdelle
% pas de prise en compte de ETG (faible transport induit)
shear          = (profli.qjli(:,5) - profli.qjli(:,end - 5)) ./ (profli.xli(:,5) - profli.xli(:,end - 5))  .* ...
                  profli.xli(:,11) ./ profli.qjli(:,11);
b2             = profli.bpol(:,11) .^ 2 + (profli.fdia(:,11) .* profli.ri(:,11)) .^ 2;
rhos           = 4.57e-3 .* sqrt(A .* tep(:,11) ./ 1e3 ./ b2);
% en presence de gradient de densite
rolte_tem_gn   = 20 ./ 9 ./ Kt + (2/3) .* rolne + Kt ./ 2 .* (1 -  rolne) .^ 2; 
%rolte_etg_cr   = (1.33 + 1.91 .* abs(shear) ./ profli.qjli(:,11)) .* (1 + profli.zeff(:,11) .* tep(:,11) ./ max(13.6,tip(:,11)));
rolti_itg_gn       = (2/3) .* rolni + (20/9) .* tip(:,11) ./ max(13.6,profli.zeff(:,11) .* tep(:,11)) +  ...
                      0.5 .* (1  - rolni + 0.25 .* rolni .^ 2) .* profli.zeff(:,11) .* tep(:,11) ./ max(13.6,tip(:,11));
%sans gradient de densite
% Article de Casati 
lambda_e        = 0.25 +  2/3 .* abs(shear);
rolte_tem_0gn   =  lambda_e ./ profli.ftrap(:,11) .*  (4.4 + 1.1 .*  ...
                   min(1.5,profli.zeff(:,11) .* tep(:,11) ./ max(13.6,tip(:,11)))); 
% article Fourment (sans le 1/2 pour coller aux experiences avec les ITGs)
sf2             = (sign(shear) >= 0) .* (1.1 + 1.4 .* abs(shear) + 1.9 .* abs(shear) ./ profli.qjli(:,11)) + ...
                  (sign(shear) < 0) .* (0.9 + 1.6 .* abs(shear)  + 9.9 .* abs(shear) ./ profli.qjli(:,11)); 
rolti_itg_0gn  =  (4/3) .* (1 + tip(:,11) ./ max(13.6,profli.zeff(:,11) .* tep(:,11))) .* sf2;

% transition (tend vers 0 quand le gradient de densite est negligeable)
trans_e = (1 + tanh(rolne - 2 .* (1 + profli.zeff(:,11) .* tep(:,11) ./ max(13.6,tip(:,11))))) ./ 2;
rolte_cr = rolte_tem_gn .* trans_e + (1 - trans_e) .* rolte_tem_0gn;
trans_i = (1 + tanh(rolni - 2 .* (1 + tip(:,11) ./ max(13.6,profli.zeff(:,11) .* tep(:,11))))) ./ 2;
rolti_cr = rolti_itg_gn .* trans_i + (1 - trans_i) .* rolti_itg_0gn;

% les coefficients de transport
xie  = (stiff .* profli.ftrap(:,11) .* sqrt(Kt) .* max(0, sqrt(rolte - rolte_cr) .* (rolte > rolte_cr)) + 1);
xii  = (stiff  .* max(0, sqrt(rolti - rolti_cr) .* (rolti > rolti_cr)) + 1);
%  terme croise de ITG dans Xie
% ref  : F. Ryter et al NF 2011 113016
% attention la sortie  est Kii /Kie !
xiioxie = xii ./ (xie +  0.5 .* tep(:,11) ./ max(13.6,tip(:,11)) .* (xii - 1)) .* profli.nip(:,11) ./ profli.nep(:,11);
