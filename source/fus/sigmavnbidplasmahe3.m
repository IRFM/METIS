% cette fonction tabule <sigma*v> pour l'interraction faisceau (D) - plasma (He3)
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
phys = cphys;

% le premier indice est l'energie du faisceau (eV)
enbi = logspace(3,6.5,101)';
% le deuxieme est la temperature du plasma cible (eV)
eti  = logspace(3,5,51);
% resevation de la sortie
svnbidhe3 = NaN .* ones(101,51) ;

% reduce masse in eV
% for pure D and He3
mu = 3.016029322 * 2.0141017778 / (3.016029322 + 2.0141017778) * phys.ua;
la = sqrt(2 .* phys.e .* 1000 ./ mu);
lb = sqrt(2 .* phys.e .* 4.7e6 ./ mu); % to prevent NaN litle smaller
v  = logspace(log10(la),log10(lb),100001);

% boucle sur les energie de faisceau
for k=1:length(enbi)

% boucle sur les temperature plasma
   for l=1:length(eti)
      enbic         = enbi(k);
      etic          = eti(l);
      vb            = sqrt(2 .* phys.e .*  enbic./ phys.ua ./ 2.0141017778);
      vth           = sqrt(2 .* phys.e .*  etic ./ phys.ua ./ 3.016029322);
      svcb          = trapz(v,svb2int_dhe3(v,vth,vb));
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

