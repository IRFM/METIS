% reference : P.C. Stangeby NF 51 (2011) 063001
% with fmom = 1 and fcond = 1;
% Test selon Stangeby problem 5.2
% qpar = 1e8 * ones(1,101)'; lc = 100* ones(1,101)'; zeff =ones(1,101)'; meff = 2.5 .* ones(1,101)'; fpe =1.* ones(1,101)'; fei = 2.* ones(1,101)'; neu = logspace(18,20,101)';tite =  ones(1,101)';
% [teu,tet,net,indbad] = z0twopoints(qpar,neu,lc,zeff,meff,fpe,fei,tite);
% figure;
% subplot(3,1,1)
% plot(neu,teu)
% hold on
% plot(3.68e19,120,'or');
% title('test of 2 points model')
% ylabel('teu (eV)')
% subplot(3,1,2)
% semilogy(neu,tet)
% hold on
% plot(6e19,10,'or',3.6e19,20,'or');
% ylabel('tet (eV)')
% subplot(3,1,3)
% semilogy(neu,net)
% ylabel('net (eV)')
%% Stangeby Nucl. Fusion 51 (2011) 063001 (16pp)
%% other test
% Pheat = [4.7 9.2 12] .* 1e6;
% nt    = 10^20 .* [ 2.6 3.0 3.3];
% qtot  = [23.5 36 56] .* 1e6;  % take into account fpower correction
% gamma_eff = [5.16 6.85 9.69];
% neu       = 10^20 .* [ 0.21 0.25 0.35];
% Teu   = [40 50 75];
% Pouter_leg  = [1.27 1.9 2.75] .* 1e6;
% fpower      = [0.43 0.63 0.56];
% fmom         = [1 1 1];
% fcond         = [1 1 1];
% f_ei       = [6.19 4.8 2.51];
% TDu        = [200 190  113 ]; % ?
% Tiu        = [160 160 130];
% fpe       = [0.16 0.2 0.82];
% fie       = 1 + Tiu ./ Teu;
% Tet   = [5 5 5];
% meff = [2 2 2];
% zeff = [2 2 2];
% R = 1.65; q= 3.9;
%  lc   = pi .* R .* q; 
% les formules sont differentes pour la vitesse du son
%[teu_out,tet_out,net_out,indbad] = z0twopoints(qtot',neu',lc',zeff',meff',fpe',f_ei',1,gamma_eff');
%nt,net_out
%Teu,teu_out
%Tet,tet_out


function [teu,tet,net,indbad,noconv] = z0twopoints(qpar,neu,lc,zeff,meff,fpe,fei,tite,gamma)

% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)

% securite
neu  = max(1e13,neu);
qpar = max(1,qpar);
noconv = 0;

% thermal // Spitzer heat conductivity coefficient
koe = 2000 ./ zeff;
% sheath factor (include ion and electron)
if nargin < 9
	gamma = 7;
end

% coefficients
A = 4 ./ fei ./ neu;
B = qpar ./ gamma ./ phys.e ./ sqrt((zeff + tite) .* phys.e ./ (phys.ua .* meff));
%B = qpar ./ gamma ./ phys.e ./ sqrt(2 .* phys.e ./ (phys.ua .* meff));
C = (7/2) .* fpe .* qpar .* lc ./ koe;
D = (A .* B) .^ (7/2);

% resolution
c0 = -D;
c1 = C;
c2 = zeros(size(qpar));
c3 = ones(size(qpar));
[x1,x2,x3]= zpolyroot(c0,c1,c2,c3);

% selection
xs = cat(2,x1,x2,x3);
mask = (abs(imag(xs)) < sqrt(eps)) & ( real(xs) >= 0) & (isfinite(xs));
data_ok = double(sum(mask,2) > 0);
%tau  = sum(real(xs).* mask,2) ./ max(1,sum(mask,2));
tau  = sum(abs(xs).* mask,2) ./ max(1,sum(mask,2));

% (x1 - tau) .* (x2 - tau) .* (x3 - tau) 
% calcul des sorties
tet = tau .^ (4/7);
net = B ./ max(eps,tet .^ (3/2));
teu = A .* net .* tet;

% points non valide (changement de formule)
% cas tet << teu
indeps = find(tet < eps);
if ~isempty(indeps) 
    %fprintf('TeT < eps: %d / %d \n',length(indeps),length(tet));
    %noconv = 1;
    %fprintf('&');
    teu(indeps) = (C(indeps) + 1e-1 .^ (7/2)) .^ (2/7);
    %tet(indeps) = (qpar(indeps) ./ (gamma .* phys.e * (teu(indeps) ./ A(indeps)) .* sqrt(2 .* phys.e ./ (phys.ua .* meff(indeps)))) ) .^ 2;     
    tet(indeps) = 0.1;     
    net(indeps) = teu(indeps) ./ tet(indeps) ./ A(indeps);
end

%  % cas tet ~ teu
%  indone = find(tet > teu);
%  if ~isempty(indone)
%      %fprintf('Teu>=TeT : %d / %d \n',length(indone),length(tet));
%      %fprintf(':');
%      noconv = 1;
%      net(indone) = 1 ./ A(indone);
%      tet(indone) = (B(indone) ./ net(indone)) .^ (2/3); 
%      %teu(indone) = tet(indone);
%      teu1 = tet(indone) .*  (1 + (2/7) .* C(indone) ./ (tet(indone) .^ (7/2)));
%      teu2 = (C(indone) + tet(indone) .^ (7/2)) .^ (2/7);
%      err1 = abs(teu1 .^ (7/2) - tet(indone) .^ (7/2) - C(indone));
%      err2 = abs(teu2 .^ (7/2) - tet(indone) .^ (7/2) - C(indone));
%      teu(indone) = teu1 .* (err1 <= err2) + teu2 .* (err1 > err2);
%      net(indone) = teu(indone) ./ tet(indone) ./ A(indone);
%      tet(indone) = (B(indone) ./ net(indone)) .^ (2/3); 
%      teu1 = tet(indone) .*  (1 + (2/7) .* C(indone) ./ (tet(indone) .^ (7/2)));
%      teu2 = (C(indone) + tet(indone) .^ (7/2)) .^ (2/7);
%      err1 = abs(teu1 .^ (7/2) - tet(indone) .^ (7/2) - C(indone));
%      err2 = abs(teu2 .^ (7/2) - tet(indone) .^ (7/2) - C(indone));
%      teu(indone) = teu1 .* (err1 <= err2) + teu2 .* (err1 > err2);
%      tet(indone) = (B(indone) ./ net(indone)) .^ (2/3); 
%      net(indone) = (teu(indone) ./ tet(indone)) ./ A(indone);
%  end

% verification
% equation 1
err1 = eps + abs(4 .* net .* tet - fei .* neu .* teu) ./ (4 .* net .* max(eps,tet) + fei .* neu .* teu) .* 2;
% equation 2
err2 = eps + abs(qpar - gamma .* phys.e .* net .* tet .* sqrt((zeff + tite) .* phys.e .* tet ./ (phys.ua .* meff))) ./ qpar;
% equation 3
err3 = eps + abs(teu .^ (7/2) - tet .^ (7/2) - (7/2) .* fpe .* qpar .* lc ./ koe) ./ ((7/2) .* fpe .* qpar .* lc ./ koe);
% backup for invalid
indbad = find((data_ok == 0) | (err1 > 1e-3) | (err2 > 1e-3) | (err3 > 1) | ...
              (tau == 0) | ~isfinite(tau));
          
% add security on data type
teu = max(eps/10,real(teu));
tet = max(eps/10,real(tet));
net = max(eps/10,real(net));

% pas de plot
return

%
% valeur max de teu si tet ~ 0
temax0 = ((7/2) .* fpe .* qpar .* lc ./ koe) .^ (2/7);

ind=1:length(qpar);
data_wrong = double((data_ok == 0) | (err1 > 1e-3) | (err2 > 1e-3) | (err3 > 1));
data_wrong(data_wrong == 0) = NaN;
data_ok(data_ok == 0) = NaN;
figure(21); clf;
subplot(4,1,1)
semilogy(ind,err1,'b',ind,err2,'r',ind,err3,'g',ind,data_ok,'o',ind,data_wrong,'*');
subplot(4,1,2)
plot(ind,teu,ind,tet,ind,temax0);
subplot(4,1,3)
semilogy(ind,neu,ind,net);
subplot(4,1,4)
plot(ind,qpar);

drawnow
%keyboard
