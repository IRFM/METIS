% cette fonction tabule <sigma*v> pour l'interraction faisceau (D) - plasma (T)
% ref :  D. R. Mikkelsen, NF vol 29 p 1113, 1989
function sigmavnbidplasmahe3

% plot de verification 
ecm  = linspace(0,4.7e6,1001);
sddt  =  z0csdhe3(ecm./1e3);
figure
loglog(ecm,sddt,'r');
xlabel('Ecm (eV)');
ylabel('Sigma (mb)');
drawnow




% constantes utiles
mp          =   1.6726485e-27;            % masse au repos du proton (kg)
e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)


% le premier indice est l'energie du faisceau (eV)
enbi = logspace(3,6.5,101)';
% le deuxieme est la temperature du plasma cible (eV)
eti  = logspace(3,5,51);
% resevation de la sortie
svnbidhe3 = NaN .* ones(101,51) ;


la            = sqrt(2 .*e .*  1000./ mp ./ 3) % borne inferieur 1keV
lb            = sqrt(2 .*e .*  4.7e6./ mp ./ 2) % borne superieur 4700 keV

% boucle sur les energie de faisceau
for k=1:length(enbi)

% boucle sur les temperature plasma
   for l=1:length(eti)
      enbic         = enbi(k);
      etic          = eti(l);
      vb            = sqrt(2 .*e .*  enbic./ mp ./ 2);
      vth           = sqrt(2 .*e .*  etic ./ mp ./ 3);
      la            = sqrt(2 .*e .*  1000./ mp ./ 3); % borne inferieur 1keV
      lb            = sqrt(2 .*e .*  4.7e6./ mp ./ 2); % borne superieur 4700 keV
      v             = logspace(log(la)/log(10),log(lb)/log(10),100001);
      svcb          = trapz(v,svb2int(v,vth,vb));
      %warning off
      %svc           = quad('svb2int',la,lb,1e-28,0,vth,vb);
      %warning on
      %svc           = svc ./ (vb .* vth .* sqrt(pi)) ;  
      svcb           = svcb ./ (vb .* vth .* sqrt(pi)) ;  
      svnbidhe3(k,l)  = svcb;
      fprintf(' svnbidhe3(%g,%g) = %g  @ %g   (m^3*s)\n',vb,vth,svcb,svcb);
   end
end
save sigmav_nbi_d_plasma_he3 enbi eti svnbidhe3

