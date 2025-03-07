function [tcvtemp,tcvprof,tcveq,tcvpsi] = zlecttcv(choc,racine,relect)
%
% -- lecture des donnees ftu --
% [tcvtemp,tcvprof,tcveq,tcvpsi] = zlecttcv(choc)
%
% input
% choc   : numero de choc TCV
%
% Auteur V. Basiuk
% 3 decembre 2004
% version 3.0
%
%

directory = [racine,int2str(choc)];
if ~exist(racine,'dir')
  eval(['!mkdir ',racine])
end

if ~exist(directory,'dir')
  eval(['!mkdir ',directory])
end

filetemp = [directory,'/tcvtemp.mat'];
fileprof = [directory,'/tcvprof.mat'];
fileeq = [directory,'/tcveq.mat'];
%
% relect = 1 : relecture des fichiers
%
relect = 1;
if relect == 1
%
% donnees plasmas
%
  liste                    = [];
  liste{1}               = 'results::I_p';
  liste{end+1}      = 'magnetics::r_axis';
  liste{end+1}      = 'magnetics::z_axis';
  liste{end+1}      = 'magnetics::rbphi';
  liste{end+1}      = 'results::l_i';
  liste{end+1}      = 'magnetics::vloop[*,"001"]';
  liste{end+1}      = 'results::beta_pol';
  liste{end+1}      = 'results::j_tor';
  liste{end+1}      = 'results::q_psi';
  liste{end+1}      = 'results::r_contour';
  liste{end+1}      = 'results::z_contour';
%
% lecture de Te, rho associe
%
%
% donn�s brutes temps, te ne et psi (liuqe 1)
%
  liste{end+1}      = 'results::thomson:times';
  liste{end+1}      = 'results::te_thomson:foo';
  liste{end+1}      = 'results::ne_thomson:foo';
  liste{end+1}      = 'results::thomson:psiscatvol';
%
% donn�s fitter de Te et ne
%
  liste{end+1}      = 'results::proffit.avg_time:teft';
  liste{end+1}      = 'results::proffit.avg_time:time';
  liste{end+1}      = 'results::proffit.avg_time:neft';
%
% lecture puissance injectee
%
  liste{end+1}      = 'results::toray.input:freq';
  liste{end+1}      = 'results::toray.input:phi_toray';
  liste{end+1}      = 'results::toray.input:theta_toray';
  liste{end+1}      = 'results::toray.input:P_GYRO';
%
% lecture de psi
%
  liste{end+1}      = 'results::psi';
%
% donnees geometriques
%
  liste{end+1} = 'results::kappa_95';
  liste{end+1} = 'results::kappa_edge';
  liste{end+1} = 'results::delta_95';
  liste{end+1} = 'results::delta_edge';
%
% lecture impuretes (<zeff>profil de zeff)
%
  liste{end+1} = 'results::zx:foo';
  liste{end+1} = 'results::zx:error_bar';
  liste{end+1} = 'results::x_tomo_zx:zx';
  liste{end+1} = 'results::x_tomo_zx:position';
%
% psi toolbox equilibre, rs et zs coordonnees radial et vertical des surfaces de flux en
% fonction du temps
%
  liste{end+1} = 'results::psitbx:rs';
  liste{end+1} = 'results::psitbx:zs';
%
% donn�s manquantes
%
  liste{end+1} = 'results::psitbx:vol';
  liste{end+1} = 'results::total_energy';
  liste{end+1} = 'results::surface_flux';

   [tcvtemp,tcvprof,tcveq]   = cgcgettcv(choc,liste);
%
% sauvegarde donn�s temporelles
%
  eval(['save ',directory,'/tcvtemp tcvtemp'])
%
% sauvegarde donn�s profils
%
  eval(['save ',directory,'/tcvprof tcvprof'])
%
% sauvegarde donn�s equilibres
%
  eval(['save ',directory,'/tcveq tcveq'])
else
  load(filetemp)
  load(fileprof);
  load(fileeq);
end


tcvpsi = [];



