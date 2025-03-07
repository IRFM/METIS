% nbar   = nl /2/a (10^20 m^-3)
% ip     = courant  plasma (MA)
% Bt     = champ toroidal / R*Bt soit la rigidite magnetique (T)
% ploss  = puissance perdue a l'etat stationnaire (MW)
% pin    = puissance de chauffage y compris ohmique (MW)
% R      = grand rayon plasma (m)
% K      = elongation plasma (b/a)
% Ka     = elongation plasma (S/pi/a^2)
% ep     = a/R (su)
% meff   = masse effective (su)
% zeff   = charge effective (su)
% loi de confinement et de transition L-H
function [wrlw,plossl2h,tauthl,tauh,tauhe_l,tauhe_h,tau_na,tau_nas,tau_ped_sc] =  ...
         zscale0(nbar,ip,Bt,ploss,pin,R,K,a,meff,zeff,Vp,Sp,q95,sext,modescale,l2hsc,pohm,nem,eta,ae,tau_bohm,prad_bol,delta,Ka95,d95,time)
     
% compatibilite ascendante
if nargin < 23    
    delta = zeros(size(nbar));
end

% selon la fonction 
ne   = 10 .* nbar;
ep   = a ./ R;
b    = a .* K;
Ka   = Vp ./ (2*pi^2.*R.*a.^2);
qcyl  = 5 .* Ka .* a .^ 2 .* Bt ./ R ./ ip;
qcyl_star  = 5 .* (1 + Ka .^ 2) ./ 2 .* a .^ 2 .* Bt ./ R ./ ip;
Fq    = q95 ./ qcyl;
Btout = Bt .* sqrt((R ./ (R + a)) .^ 2 + (a ./ R ./ q95) .^ 2);
fa    = 1 - (2 ./ ( 1 + max(1.1,meff))) .^ 0.5;
fa(fa<=0) = 0.1 .* meff(fa<0);
Fa    = 0.1 .* meff ./ fa;
Hg    = nbar .* pi .* K .* a .^ 2 ./ ip;
nsat  = 0.06 .* ip .* R .* sqrt(meff) ./ K ./ a .^ 2.5;
nstar = min( nsat,nbar); 
F_Gr    = max(eps,nbar .* pi .*  a .^ 2 ./ ip);
frad    = min(1,prad_bol ./ max(eps,ploss .* (ploss > prad_bol) + pin .* (ploss <= prad_bol)));

%
% transition d'apres J.G. Cordey rapport JET-P(85)28
rap   = min(1,max(0,erfc(max(pin,pohm) ./ max(ip/1e3,pohm) - 1)));


% tau neo alcator; Wesson p 181
tau_na  = 0.07 .* nbar  .* a .* R .^ 2 .* qcyl_star .* sqrt(meff);
tau_nas = 0.07 .* nstar .* a .* R .^ 2 .* qcyl_star .* sqrt(meff);

% rlw
wrlw  = 2.6e-2 .* ne .^ (3/4) .* zeff .^ (1/4)  .*  Bt .^ (1/2) .* ip  .^ (1/2) .* (R .* a .* b) .^ (11/12) + ...
    1.2e-2 .* ip .* (R .* a .* b) .^ (1/2) .* ploss .* zeff .^ (-1/2); %MJ
wrlw  = wrlw .* 1e6;
switch l2hsc
    case 1
        % loi ITER LH02noZeff (Ryter,PPCF,2002,A415)
        plossl2h   = 0.042 .* nbar .^ 0.73 .* Btout .^0.74 .* sext .^0.98;
    case 2
        % loi ITER LH02Zeff (Takizuka,PPCF, 2004)
        plossl2h   = 0.072 .* nbar .^ 0.7 .* Btout .^0.7 .* sext .^0.9 .* (zeff ./ 2) .^0.7 .* Fa .^ 0.5;
    case 3
        % loi Y R Martin, Journal of Physics 2008 Conference series 123 p 012033
        plossl2h   = 0.0488 .* nbar .^ 0.717 .* Bt .^ 0.803 .* sext .^ 0.941;
    case 4
        % loi NLM-7 ref : A. Murari et al, NF 52 (2012) p 0630016-
        plossl2h   = 0.074 .* nbar .^ 0.776 .* Bt .^ 0.728 .* sext .^ 1.02 .* Ka .^ -1.266;
    case 5
        % loi NLM-11 ref : A. Murari et al, NF 52 (2012) p 0630016-
        plossl2h   = 0.068 .* nbar .^ 0.765 .* Bt .^ 0.764 .* sext .^ 1.035 .* Ka .^ -1.152 .* q95 .^ -0.051 .* ip .^ -0.034;
    case 10
        % multi machine scaling from reference : A. Murari et al, NF 53 (2013) p 043001-043013, formula 6
        % the parametr K is not clearly defined, it seem that is b/a from ITER parameter numerical application
        plossl2h   = 5.119 .* nbar + 1.204 .* Bt + 0.026 .* sext + 1.069 .* Bt .* nbar - 1.651 - 1.106 .* K - 0.039 .* Bt .^ 2 - 0.166 .* nbar .^ 2;
    case 28
        % Low density case - Ryter et al, Nucl. Fusion 54 (2014) 083003 (9pp), equation 4
        plossl2h   = 0.36 .* ip .^ 0.27 .* Bt.^ 1.25 .* R .^ 1.23 .* (R./a) .^ 0.08;
        
    case 30
        % Metal DB: Status of TC 26: L H/H L scaling in the presence of
        % Metallic walls  JET ILW: E . Delabie, E.R. Solano, 
        % C . Maggi , J. Hillesheim ; AUG: F . Ryter , M. Cavedon  , 
        % G . Birkenmeier , U. Plank , C mod : J . Hughes 
        % case for horizontal target slide 23
        plossl2h   = 0.045 .* nbar .^ 1.08 .* Bt .^ 0.56 .* sext .^ 1.0 .* (2 ./ meff) .^ 0.96;
        
    case 31
        % Metal DB: Status of TC 26: L H/H L scaling in the presence of
        % Metallic walls  JET ILW: E . Delabie, E.R. Solano, 
        % C . Maggi , J. Hillesheim ; AUG: F . Ryter , M. Cavedon  , 
        % G . Birkenmeier , U. Plank , C mod : J . Hughes 
        % case for vertical/corner target slide 23
        plossl2h   = 1.93 .* 0.045 .* nbar .^ 1.08 .* Bt .^ 0.56 .* sext .^ 1.0 .* (2 ./ meff) .^ 0.96;
        
    otherwise
        % seuil mode H loi LH99(1)
        plossl2h   = 2.84 .* nbar .^ 0.58 .* Bt .^ 0.82 .* R .* a .^ 0.81 ./ meff; % MW
end
plossl2h   = plossl2h .* 1e6;

switch modescale
    case 0
        % loi standard ITER
        % ITERH-96P(th)
        tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        
    case 1
        % version optimisee  pour le startup
        % ITERH-96P(th)
        taul  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        
        % OH debut de choc
        % V.Lukash et all , Plamsa devices and operations, vol 13 , 2005, p 143-156
        % G. Bracco and K. Thomsen, NF, 37 (1997)
        %
        tauoh = 0.14 .* a .* R .* Hg .* Bt .^ 1.1 ./ (1 + 0.63 .* Bt .* Hg .^ 2);
        tauthl = rap .* tauoh + (1-rap) .* taul;
        
    case 2
        % scaling 2 termes recommande par ITPA
        % McDonald et al,
        % NF 47 (2007) p147
        tauthl   =  6.93e-2 .* ip .^ 0.62 .* ploss .^ (-0.53) .* ne .^ 0.64 .* R .^ 2.12 .* ep .^ 1.15;
        taupied  =  2.4e-2  .* ip .^ 1.64 .* ploss .^ (-0.44) .* ne .^ (-0.18) .* R .^ 1.03 .* ep .^ (-0.39);
        tauh     =  taupied + tauthl;
        
    case 3
        % cas sans dependance en beta (IAEA/ITPA 2004)
        % cordey core  ref : NF 43 (2003), 670-674 core fit 2 equation 8 (sans dependance en beta)
        % version etude Hoang@ITPA
        %tauthl  = 0.151 .* ip .^ 0.68 .* R .^ 2.32 .* ne .^ 0.59 .* Bt .^ 0.13 .*  ...
        %          Ka .^ -0.34 .* ep .^ 1.96 .* meff .^ 0.34  .* ploss .^ (0.42 -1);
        taupied  = exp(-3.87) .* ip .^ 1.6 .* R .^ 1.03 .* ne .^ (-0.16) .* ep .^ (-0.26) .* ploss .^ (0.6-1);
        
        % ITPA @AIEA 2004 from McDonald  PPCF 46 (2004),A215- , #11
        tauh   =  0.028 .* ip .^ 0.83 .* Bt .^ 0.07 .* R .^ 1.81 .* ne .^ 0.49 .* ...
            a .^ 0.3 .* K .^ 0.75 .* meff .^ 0.14 .* ploss .^  (-0.55);    % s
        tauthl = tauh - taupied;
    case 4
        % fit de Wdia
        
        % loi standard ITER
        % ITERH-96P(th)
        tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        
        %     % ITPA @AIEA 2004 from McDonald  PPCF 46 (2004),A215- , #11
        %     tauh   =  0.028 .* ip .^ 0.83 .* Bt .^ 0.07 .* R .^ 1.81 .* ne .^ 0.49 .* ...
        %                 a .^ 0.3 .* K .^ 0.75 .* meff .^ 0.14 .* ploss .^  (-0.55);    % s
        %     % Didier Elbeze
        %     tauthl = 0.016 .* ip .^ 0.84 .* Bt .^ 0.12 .* R .^ 1.55 .* ne .^  0.46 .* ...
        %                 a .^ - 0.83 .* (pi .* K .* a .^ 2) .^ 0.63 .* meff .^ 0.21 .* ploss .^ (-0.76);
    case 5
        % scaling EIV ITPA ref. MCDonald, NF 47 (2007) 147-174
        % ITERH-96P(th)
        tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        
        % ITERH-EIV(y,2)
        tauh   = 5.55e-2 .* ip .^ 0.75 .* Bt .^ 0.32 .* ne .^ 0.35 .* ploss .^ -0.62 .* ...
            R .^ 2 .* Ka .^ 1.14 .* ep .^ 0.76 .* meff .^ 0.06;    % s
        
    case 6
        % OH Wesson p 181
        % version optimisee  pour le startup
        % ITERH-96P(th)
        taul  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        
        % OH debut de choc
        tauoh = 0.07 .* nstar .* a .* R .^ 2 .* qcyl;
        % transition d'apres J.G. Cordey rapport JET-P(85)28
        tauthl = rap .* tauoh + (1-rap) .* taul;
        
    case 7
        % loi standard ITER
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        % de la modelisation ITM
        tauthl = tauh ./ 2;
        
    case 8
        
        % remplissage par defaut
        tauthl = mono_scaling(R,ep,Ka,Bt,nem,ip,meff,ploss,1.4215e-12);
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s
        
        tauthl = min(tauthl,tauh);
        % appel de la loi definie par l'utilisateur
        if ~isempty(which('metis_user_scaling'))
            try
                metis_user_scaling;
            catch
                error(sprintf('error in metis_user_scaling: %s',lasterr));
            end
        end
        
        %figure(51);clf;semilogy(tauoh);drawnow
        
    case 9
        
        % robust fit on L-mode data base (rms = 0.19 on whole valid points of the data base)
        tauthl  = 0.0796  .* ip .^ 0.77 .* Bt .^ 0.00 .* nbar .^ 0.46 .* ploss .^ -0.68 .* ...
            R .^ 2.20 .* Ka .^ 1.00 .* ep .^ 0.37 .* meff .^ 0.26; % s
        
        
        %     % loi standard ITER
        %     % ITERH-96P(th)
        %     tauthl_ref  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
        %           R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s
        %
        %     figure(21);clf;plot(tauthl,'r');hold on;plot(tauthl_ref,'b');drawnow
        %
        
        % pour les points ou le calcul altenatif echoue
        % ITPA @AIEA 2004 from McDonald  PPCF 46 (2004),A215- , #11
        tauh   =  0.028 .* ip .^ 0.83 .* Bt .^ 0.07 .* R .^ 1.81 .* ne .^ 0.49 .* ...
            a .^ 0.3 .* K .^ 0.75 .* meff .^ 0.14 .* ploss .^  (-0.55);    % s
    case 10
        
        % remplissage par defaut
        tauthl = mono_scaling(R,ep,Ka,Bt,nem,ip,meff,ploss,1.4215e-12);
        
        % fit robuste de O. Sauter and Y. martin in Nuclear Fusion Vol 40 (2000) p 955-964
        % qeng = qcyl , formule 9 , rms 16% sur la base complete, gyro-bohm
        tauh = 0.0307 .* a .* Ka .* Bt .* sqrt(ne) ./ qcyl .* (ploss ./ Vp) .^ -(2/3);
        
        % securite
        tauthl = min(tauthl,tauh);
        
    case 11
        
        % Elbeze 32nd EPS 2005 P1.033, Tarragona (EIV) Bohm
        tauthl = 0.0041  .* ip .^ 0.73 .* Bt .^ 0.16 .* ne .^ 0.43 .* ploss .^ -0.66 .* ...
            R .^ 1.97 .* (pi .* a .^ 2 .* Ka) .^ 1.14 .* a .^ -2.18 .* meff .^ 0.39; % s
        
        
        %  figure(21);clf;plot(1:length(ploss),23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
        %           R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2,'b',1:length(ploss),tauthl,'r');drawnow
        
        % fit robuste de O. Sauter and Y. martin in Nuclear Fusion Vol 40 (2000) p 955-964
        % qeng = qcyl , formule 9 , rms 16% sur la base complete, gyro-bohm
        tauh = 0.0307 .* a .* Ka .* Bt .* sqrt(ne) ./ qcyl .* (ploss ./ Vp) .^ -(2/3);
        
        % securite
        tauthl = min(tauthl,tauh);
        
    case {12,13}
        
        % avec limitation deu contenu en energie pour "advanced inductive scenario"
        
        % deriver de la definition de betan et du scaling en P^-1 de la reference :
        % T. Luce et al, NF 54 (2014) 013015
        %taumax = 0.33 .* ip .* a .* Ka  .* F_Gr .^(1/3) .* R .* Bt ./ max(1,ploss);
        % a dependance en F_Gr est mieux capturee par le scaling de Lang de la reference :
        %  Lang P.T. et al 2012 24th IAEA Fusion Energy Conf.,
        %  (San Diego, CA, 2012) EX/P4-01 http://www-naweb.iaea.org/napc/physics/FEC/FEC2012/papers/33 EXP401.pdf
        taumax = 0.33 .* ip .* a .* Ka  .* F_Gr .^ (-0.22 .* log(F_Gr)) .* R .* Bt ./ max(1,ploss);
        
        % Elbeze 32nd EPS 2005 P1.033, Tarragona (EIV) Bohm
        tauthl = 0.0041  .* ip .^ 0.73 .* Bt .^ 0.16 .* ne .^ 0.43 .* ploss .^ -0.66 .* ...
            R .^ 1.97 .* (pi .* a .^ 2 .* Ka) .^ 1.14 .* a .^ -2.18 .* meff .^ 0.39; % s
        
        
        %  figure(21);clf;plot(1:length(ploss),23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
        %           R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2,'b',1:length(ploss),tauthl,'r');drawnow
        
        % fit robuste de O. Sauter and Y. martin in Nuclear Fusion Vol 40 (2000) p 955-964
        % qeng = qcyl , formule 9 , rms 16% sur la base complete, gyro-bohm
        tauh = 0.0307 .* a .* Ka .* Bt .* sqrt(ne) ./ qcyl .* (ploss ./ Vp) .^ -(2/3);
        %t=1:length(a);
        %figure(21);clf;plot(t,taumax,'r',t,tauh,'b');drawnow
        % limitation
        tauh = min(taumax,tauh);
        % securite
        tauthl = min(tauthl,tauh);
        
    case {14,15}
        % ref : A. Murari et al, Nuclear Fusion 57 (2017) 126017-
        % Tau_H1: equation 16 and Tau_L1: equation 21
        tauh   = 2.73e-2 .* ip .^ 0.575 .* Ka .^ 1.348 .* R .^ 2.006 .* ploss .^ -0.441 .* ne .^ 0.423 ./ ...
            (1 + exp(- 0.989 .* (frad .* ne ./ Bt) .^ -0.440));
        
        %
        tauthl = 7.4e-2 .* ip .^ 0.465 .* R .^ 2.012 .* ep .^ 0.946 .* Ka .^ 1.193 .* ploss .^ -0.353 .*  ...
            exp(-0.235 .* frad .^ 0.660 ./ meff .^ 1.279 ./ ip .^ 0.646);
        
        
        if modescale == 15
            % deriver de la definition de betan et du scaling en P^-1 de la reference :
            % T. Luce et al, NF 54 (2014) 013015
            %taumax = 0.33 .* ip .* a .* Ka  .* F_Gr .^(1/3) .* R .* Bt ./ max(1,ploss);
            % a dependance en F_Gr est mieux capturee par le scaling de Lang de la reference :
            %  Lang P.T. et al 2012 24th IAEA Fusion Energy Conf.,
            %  (San Diego, CA, 2012) EX/P4-01 http://www-naweb.iaea.org/napc/physics/FEC/FEC2012/papers/33 EXP401.pdf
            taumax = 0.33 .* ip .* a .* Ka  .* F_Gr .^ (-0.22 .* log(F_Gr)) .* R .* Bt ./ max(1,ploss);
            % limitation
            tauh = min(taumax,tauh);
            % securite
            tauthl = min(tauthl,tauh);
        end
        
        %figure(21);clf;plot(ploss,tauh,'r',ploss,tauthl,'b');drawnow
        
    case 16 
        % Scalings for tokamak energy confinement, 
        % P.N. Yushmanov et al 1990 Nucl. Fusion 30 1999  (ITER-89P)
        % equation 3 in L mode and  ITERH-98P(y,2) in H-mode
        % both using Ploss for high radiative scenario (with option for
        % Ploss in METIS.
        %
        % ITER-89P
        tauthl  = 0.048 .* meff .^ 0.5 .* ip .^ 0.85 .* R .^ 1.2 .* a .^ 0.3 .* ...
                  K .^ 0.5 .* nbar .^ 0.1 .* Bt .^ 0.2 .* ploss .^ -0.5;
        
        % ITERH-98P(y,2)
        tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s

    case 17
        % Scalings for tokamak energy confinement, 
        % P.N. Yushmanov et al 1990 Nucl. Fusion 30 1999  (ITER-89P)
        % equation 3 in L mode and  add pedestal in H-mode
        % both using Ploss for high radiative scenario (with option for
        % Ploss in METIS  + triangularity (Ka)
        % ITER-89P
        tauthl  = 0.048 .* meff .^ 0.5 .* ip .^ 0.85 .* R .^ 1.2 .* a .^ 0.3 .* ...
                  K .^ 0.5 .* nbar .^ 0.1 .* Bt .^ 0.2 .* ploss .^ -0.5;
        %
        % cordey core  ref : NF 43 (2003), 670-674 core fit 2 equation 8 (sans dependance en beta)
        % version etude Hoang@ITPA
        taupied  = exp(-3.87) .* ip .^ 1.6 .* R .^ 1.03 .* ne .^ (-0.16) .* ep .^ (-0.26) .* ploss .^ (0.6-1);
        tauh = tauthl + taupied;
        
     case 18
        % Scalings for tokamak energy confinement, 
        % P.N. Yushmanov et al 1990 Nucl. Fusion 30 1999  (ITER-89P)
        % equation 3 in L mode and  ITERH-98P(y,2) in H-mode
        % both using Ploss for high radiative scenario (with option for
        % Ploss in METIS.
        %
        % ITER-89P
        tauthl  = 0.048 .* meff .^ 0.5 .* ip .^ 0.85 .* R .^ 1.2 .* a .^ 0.3 .* ...
                  K .^ 0.5 .* nbar .^ 0.1 .* Bt .^ 0.2 .* ploss .^ -0.5;
        
        % Sizing up plasmas using dimensionless parametersa…
        % C. C. Petty
        % Phys. Plasmas 15, 080501 (2008)
        % https://doi.org/10.1063/1.2961043
        % equation 36
        tauh   = 0.052 .* ip .^ 0.75 .* Bt .^ 0.30 .* ne .^ 0.32 .* ...
                 ploss .^ -0.47 .* R .^ 2.09 .* Ka  .^ 0.88 .* ep .^ 0.84;    % s (meff exponent is 0)

     case 19
        % Scalings for tokamak energy confinement, 
        % P.N. Yushmanov et al 1990 Nucl. Fusion 30 1999  (ITER-89P)
        % equation 3 in L mode and  ITERH-98P(y,2) in H-mode
        % both using Ploss for high radiative scenario (with option for
        % Ploss in METIS.
        %
        % ITER-89P
        tauthl  = 0.048 .* meff .^ 0.5 .* ip .^ 0.85 .* R .^ 1.2 .* a .^ 0.3 .* ...
                  K .^ 0.5 .* nbar .^ 0.1 .* Bt .^ 0.2 .* ploss .^ -0.5;
        
        % pedestal will add later
        tauh   = tauthl;

   case 20
        
        
        % Elbeze 32nd EPS 2005 P1.033, Tarragona (EIV) Bohm
        tauthl = 0.0041  .* ip .^ 0.73 .* Bt .^ 0.16 .* ne .^ 0.43 .* ploss .^ -0.66 .* ...
            R .^ 1.97 .* (pi .* a .^ 2 .* Ka) .^ 1.14 .* a .^ -2.18 .* meff .^ 0.39; % s
        
        % new  ITPA scaling 2020: ITPA20
        % The updated ITPA global H-mode confinement database,Geert Verdoolaege et al 2021 Nucl. Fusion, https://doi.org/10.1088/1741-4326/abdb91
        tauh   =  0.053 .* ip .^ 0.98 .* Bt .^ 0.22 .* R .^ (1.71-0.35) .* ne .^ 0.24 .* ...
            a .^ 0.35 .* Ka95 .^ 0.8 .* meff .^ 0.2 .* ploss .^  (-0.669) .* (1 + abs(d95)) .^ 0.36;    % s
        % attention il s'agit de Ka_95 et delta_95
        warning('Using formulation of ITPA20 from preprint using Ka95 and d95: consider to switch to option 21 to used published value of thi scaling !')
        % securite
        tauthl = min(tauthl,tauh);
        
    case 21
        
        
        % Elbeze 32nd EPS 2005 P1.033, Tarragona (EIV) Bohm
        tauthl = 0.0041  .* ip .^ 0.73 .* Bt .^ 0.16 .* ne .^ 0.43 .* ploss .^ -0.66 .* ...
            R .^ 1.97 .* (pi .* a .^ 2 .* Ka) .^ 1.14 .* a .^ -2.18 .* meff .^ 0.39; % s
        
        % new  ITPA scaling 2020: ITPA20
        % The updated ITPA global H-mode confinement database,Geert Verdoolaege et al 2021 Nucl. Fusion, https://doi.org/10.1088/1741-4326/abdb91
        tauh   =  0.053 .* ip .^ 0.98 .* Bt .^ 0.22 .* R .^ (1.71-0.35) .* ne .^ 0.24 .* ...
            a .^ 0.35 .* Ka .^ 0.8 .* meff .^ 0.2 .* ploss .^  (-0.669) .* (1 + abs(delta)) .^ 0.36;    % s
        % securite
        tauthl = min(tauthl,tauh);
        
        
end

% volume du plasma (rapport groupe dimensionnement)
% de JT60
%raphese    = 3 + 7 ./ (ne/3);
% iter physics basis p 2226 (effectif avec recyclage)
tauhe  = 0.74 .* Vp .^ 0.7 .* ip  .^ 0.31 .* (pin ./ nbar).^ -0.57;
tauhe_l = tauhe;
tauhe_h = tauhe ./ max(1e-3,tauthl) .* tauh;
% volume du plasma (rapport groupe dimensionnement)
% de JT60
%raphese    = 3 + 7 ./ (ne/3);
% iter physics basis p 2226 (effectif avec recyclage)
%tauhe  = 0.74 .* Vp .^ 0.7 .* ip  .^ 0.31 .* (pin ./ nbar).^ -0.57;
%tauhe_l = (tauhe + raphese .* tauthl) ./ 2 ;
%tauhe_h = (tauhe ./ max(1e-3,tauthl) .* tauh  + raphese .* tauh) ./ 2;

% from scaling 2 termes recommande par ITPA
% McDonald et al, 
% NF 47 (2007) p147
tau_ped_sc = 2.4e-2  .* ip .^ 1.64 .* ploss .^ (-0.44) .* ne .^ (-0.18) .* R .^ 1.03 .* ep .^ (-0.39);

% use excternal data if available
if isappdata(0,'TAUE_EXP')
    tau_ext  = getappdata(0,'TAUE_EXP');
    if isfield(tau_ext,'temps') && ~isempty(tau_ext.temps)
        if isfield(tau_ext,'tau_bohm') && ~isempty(tau_ext.tau_bohm)
            tau_bohm = import_external_tau(tau_ext.temps,tau_ext.tau_bohm,tau_bohm,time);
        end
        if isfield(tau_ext,'tauthl') && ~isempty(tau_ext.tauthl)
            tauthl = import_external_tau(tau_ext.temps,tau_ext.tauthl,tauthl,time);            
        end
        if isfield(tau_ext,'tauh') && ~isempty(tau_ext.tauh)
            tauh = import_external_tau(tau_ext.temps,tau_ext.tauh,tauh,time);            
        end
        if isfield(tau_ext,'tauhe_l') && ~isempty(tau_ext.tauhe_l)
            tauhe_l = import_external_tau(tau_ext.temps,tau_ext.tauhe_l,tauhe_l,time);                        
        end
        if isfield(tau_ext,'tauhe_h') && ~isempty(tau_ext.tauhe_h)
            tauhe_h = import_external_tau(tau_ext.temps,tau_ext.tauhe_h,tauhe_h,time);                        
        end
    end
end

% securite
% utilisation de la regle ref : V. A. Belyakov et al , Plasma Devices and operations vol 11, 2003, p 193-201
taum    = R + a .* K;
tauthl  = min(max(tau_bohm,tauthl),taum);
tauh    = min(max(tau_bohm,tauh),taum);
tauhe_l = min(max(tau_bohm,tauhe_l), 10 .* taum);
tauhe_h = min(max(tau_bohm,tauhe_h), 10 .* taum);


%disp(' in scale0')
%keyboard
%function taug = mono_scaling(R,ep,Ka,Bt,nem,q95,meff,ploss,c0,m,nu)

%ee           =   1.602176462e-19;   

%  taug =  c0.^(-2./(4+m)).* 3 .* 4.^(-1./(4+m)) .* ee.^((4+3.*m)./(4+m)) .* pi.^ (2.*(2+m)./(4+m)) .*  ...
%          ep.^(4.*(2+m)./(4+m)) .* R.^(5.*(2+m) ./(4+m)) .* ...
%          ((4.*(1+Ka.^2).^2)./(4 + 2.*Ka.^2)).^(1./(4+m)).* Ka.^((2+m)./(4+m)) .*  ...
%      	Bt.^(2.*(1+m)./(4+m)) .* q95.^(-2.*nu./(4+m)) .* ...
%  	(1+ae.*eta.*f).^(-2./(4+m)) .* (1+ae.*eta) .* ...
%          nem.^((2+m)./(4+m)) .* meff.^(-m./(4+m)) .* ploss.^(-(2+m)./(4+m));

  


%  taug =  c0 .*  ...
%          ep.^(4.*(2+m)./(4+m)) .* R.^(5.*(2+m) ./(4+m)) .* ...
%          ((4.*(1+Ka.^2).^2)./(4 + 2.*Ka.^2)).^(1./(4+m)).* Ka.^((2+m)./(4+m)) .*  ...
%      	Bt.^(2.*(1+m)./(4+m)) .* q95.^ nu .* ...
%  	(1+ae.*eta.*f).^(-2./(4+m)) .* (1+ae.*eta) .* ...
%          nem.^((2+m)./(4+m)) .* meff.^(-m./(4+m)) .* ploss.^(-(2+m)./(4+m));

%   taug =  c0 .*  ...
%          ep.^(4.*(2+m)./(4+m)) .* R.^(5.*(2+m) ./(4+m)) .* ...
%          ((4.*(1+Ka.^2).^2)./(4 + 2.*Ka.^2)).^(1./(4+m)).* Ka.^((2+m)./(4+m)) .*  ...
%      	Bt.^(2.*(1+m)./(4+m)) .* q95.^ nu  .* ...
%          nem.^((2+m)./(4+m)) .* meff.^(-m./(4+m)) .* ploss.^(-(2+m)./(4+m));


function taug = mono_scaling(R,ep,Ka,Bt,nem,ip,meff,ploss,c0)

% The Quasi-Isotropic Stochastic Magnetic Field and Transport 
% p 121, The Quasi-Isotropic Stochastic Magnetic Field and Transport,O.G. Bakunin,
% Review of Plasma Physics vol 24
 taug =  c0 .^ (-2/5) .* 5.73e-18 .* ep .^(9/5) .* R .^ (13/5) .* nem .^ (3/5) .* ip .^ (2/5) .* Bt .^ (2/5) .* ...
        ((1+Ka.^2)./(2+Ka.^2)).^ (1/5) .* Ka .^ (3/5)  .* ploss .^ (-3/5);


function val_out = import_external_tau(time_in,val_in,val_origin,time_out)

val_out          = interp1_ex(time_in,val_in,time_out,'linear','extrap');
indnok           = find(~isfinite(val_out));
val_out(indnok)  = val_origin(indnok);
val_out          = max(0,val_out);


        
