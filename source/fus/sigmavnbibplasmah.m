% cette fonction tabule <sigma*v> pour l'interraction faisceau (B) - plasma (p)
% ref :  D. R. Mikkelsen, NF vol 29 p 1113, 1989
function sigmavnbibplasmah

% plot de verification 
ecm  = linspace(1e3,1e7,100001);
sdpb  =  pb11_cross_section_tentori(ecm);
figure
loglog(ecm,sdpb,'r');
xlabel('Ecm (eV)');
ylabel('Sigma (m^2)');
drawnow




% constantes utiles
phys = cphys;


% le premier indice est l'energie du faisceau (eV)
%enbi = logspace(3,6.5,101)';
% take into account cross section ranges of energy
enbi = cat(2,linspace(1e3,399e3,301),linspace(400e3,799e3,101),logspace(log10(800e3),6.5,101))';
% le deuxieme est la temperature du plasma cible (eV)
%eti  = logspace(3,6,101);
eti  = cat(2,linspace(1e3,399e3,301),linspace(400e3,799e3,101),logspace(log10(670e3),6,71))';
% resevation de la sortie
svnbibp = NaN .* ones(length(enbi),length(eti));
%
% reduce masse in eV
% for pure B11 and proton
mu = 11.0093054 * 1.00782503207 / (11.0093054 + 1.00782503207) * phys.ua; % kg
la = sqrt(2 .* phys.e .* 1000 ./ mu);
lb = sqrt(2 .* phys.e .* 9.75e6 ./ mu); % to prevent NaN litle smaller
v  = logspace(log10(la),log10(lb),200001);


% boucle sur les energie de faisceau
for k=1:length(enbi)

% boucle sur les temperature plasma
   for l=1:length(eti)
      enbic         = enbi(k);
      etic          = eti(l);
      vb            = sqrt(2 .* phys.e .*  enbic./ phys.ua ./ 11.0093054);
      vth           = sqrt(2 .* phys.e .*  etic ./ phys.mp);
      svcb          = trapz(v,svb2int_pb(v,vth,vb));
      svcb           = svcb ./ (vb .* vth .* sqrt(pi)) ;  
      svnbibp(k,l)  = svcb;
      fprintf(' svnbibp(%g,%g) = %g  (m^3*s)\n',vb,vth,svcb);
   end
end
save sigmav_nbi_b_plasma_p enbi eti svnbibp

