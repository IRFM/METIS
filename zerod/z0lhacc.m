% diagramme LH
function [x,fpout,xlh,dlh,lc,hc,acc,landau,efficiency,effacc,fjlh] = z0lhacc(flh,npar0,width,agaz,zgaz,temps,x,nep,tep,qp,Raxe,rmx,spr,vpr,Bt,plh,xlhin,dlhin,transitoire,upshift,plotonoff,upshiftmode)


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


% test des entrees
if nargin < 21
	plotonoff = 0;
end 
if nargin < 22
	upshiftmode = 'linear';
end
if npar0 < 0
	lobe = -1;
	npar0 = abs(npar0);
else
	lobe = 1;
end
% securite
npar0 = max(1,npar0);

% vecture utils
ve = ones(size(x));
vt = ones(size(tep,1),1);

% champ total
btor      = Bt * ve;
bpol      = sign(vt * x) .* rmx .* btor ./ Raxe ./ qp;
btot      = sqrt( btor .^2  + bpol .^2);
qp        = max(1,min(qp(:,end) *ve,qp));
sp        = pdederive(x,qp,2,2,2,1) ./ qp .* (vt*x); 
% calcul du up-shift
kinetic_factor = ones(size(qp));
switch upshiftmode
case '1/q'
	upshift   = (upshift * ve)  .* ((1 ./ qp - 1 ./ (max(qp,[],2) *ve)) ./ (1 ./ (min(qp,[],2) *ve) - 1 ./ (max(qp,[],2) *ve))) .^ 2;
case 'Bpol'
	upbp   =  - cumtrapz(x(end:-1:1),bpol(:,end:-1:1) ./ btor,2);
	upbp   =   upbp(:,end:-1:1) ./ max(eps,max(upbp,[],2) * ve);
	upshift   = min(10,max(-0.9,(upshift * ve)  .* upbp));
	
case 'x^2'
	upshift   = (upshift * ve)  .* (1 - (rmx ./ (rmx(:,end) * ve)) .^ 2);
case 'sqrt(x)'
	upshift   = (upshift * ve)  .* (1 - sqrt(rmx ./ (rmx(:,end) * ve)));
case 'null'
	upshift   = zeros(size(qp));
case {'newmodel','newmodel + tail','newmodel retune'}
        kinetic_factor = upshift * ve;
        kinetic_factor(kinetic_factor == 0) = 1;
        upshift   = zeros(size(qp));
case 'step@edge'
	upshift   = (upshift * ve); 
otherwise
	upshift   = (upshift * ve)  .* (1 - rmx ./ (rmx(:,end) * ve));
end

% log lambda 
lnl = 14.9 - 0.5 .* log(nep./1e20) + log(tep./1e3);

%  if lobe >0
%  	figure(21);clf;plot(x,npar0 .* (upshift + 1));drawnow
%  end

% frequence
w    = 2 .* pi .* flh;
wpe  = sqrt(nep .* phys.e .^ 2 ./ phys.me ./ phys.epsi0);
wpi  = sqrt(nep ./ (zgaz *ve) .* phys.e .^ 2 ./ phys.mp ./ (agaz * ve) ./ phys.epsi0);
wce  = phys.e ./ phys.me .* btot;

% SDP
S = 1 + wpe .^2 ./ wce .^2 - wpi .^ 2 ./ w .^ 2;
D = wpe .^ 2 ./ w ./ wce;
P = 1 - wpe .^2 ./ w .^2 - wpi .^2 ./ w .^ 2;


% Landau
%  fnbl   =  max(1,6.5 - 1.56 .* sqrt( 2 .* rmx ./ (Raxe + rmx)));
%  landau = fnbl ./ sqrt(max(30,tep) ./1e3);
landau = 6.5 ./ sqrt(max(30,tep) ./1e3);
% tail model
switch upshiftmode
    case {'newmodel','newmodel retune'}
        upshift   = zeros(size(qp));
    case 'newmodel + tail'
        upshift = max(0,landau(:,1) - npar0) * (1 - x);
end

% courbes
lc     = (npar0 + upshift) ./ (1 + rmx ./ qp ./ Raxe .* sqrt(max(eps,- P ./ S)));
hc     = (npar0 + upshift) ./ (1 - rmx ./ qp ./ Raxe .* sqrt(max(eps,- P ./ S)));
acc    = wpe ./ wce + sqrt(1 + wpe .^ 2 ./ wce .^ 2 - wpi .^2 ./ w .^2);
acc(imag(acc) ~= 0)    =  1000;
indbad  = find((hc <lc) | (hc > (6.5./sqrt(0.03))));
hc(indbad) = (6.5./sqrt(0.03));
hc(find(imag(hc))) = npar0;
hc(find(~isfinite(hc))) = npar0;

% largeur naturelle

dn0  = phys.c ./flh ./ width; % depend du grill
switch upshiftmode
    case {'newmodel','newmodel + tail','newmodel retune'}
        % d'apres les simulations completes ALOHA/C3PO/LUKE, cet effet n'existe pas
        % largeur naturelle
        dn0  = phys.c ./flh ./ width; % depend du grill
        % accessibility
        effacc  = min(1,exp((npar0 -  acc(:,end))./dn0));
        % limite de propagation (decale d'une largeur pour entrer en service quand l'onde approche la coupure)
        % si l'onde est absorbee il faut pendre landau
        factacc = min(1,exp((landau -  (acc + dn0))./dn0));
        %largeur au lieu d'absorbtion
        % natural + upshift + poloidal extention effect (assume = toroidal extention, small in the simulation)
        % meilleur choix:
        switch upshiftmode
            case 'newmodel'
                dlhabs = (dn0  + npar0 .* width ./ Raxe) .* (1 + landau ./ npar0);
            otherwise
                dlhabs = (dn0  + npar0 .* width ./ Raxe);
        end
        % caustique/refraction/ ...
        switch upshiftmode
            case {'newmodel retune'}
                factlc  = min(1,exp((landau - lc) ./ dlhabs));
                facthc  = min(1,exp((hc - landau) ./ dlhabs));
                % probability of absorbtion choix 1
                pabs = exp(-(kinetic_factor .* landau - (npar0 + dlhabs ./ 4 + upshift) ) .^ 2  ./ dlhabs .^ 2) .*  factacc .* factlc .* facthc;
            otherwise
                factlc  = min(1,exp((landau - (lc + dlhabs ./ 2)) ./ dlhabs));
                facthc  = min(1,exp(((hc - dlhabs ./ 2) - landau) ./ dlhabs));
                % probability of absorbtion choix 1
                pabs = exp(-(kinetic_factor .* landau - (npar0 + dlhabs ./ 2 + upshift) ) .^ 2  ./ dlhabs .^ 2) .*  factacc .* factlc .* facthc;
        end
        pabs = pabs ./ max(eps,(max(pabs,[],2) * ve));
        
        % 	figure(21);clf
        % 	subplot(2,2,1)
        % 	plot(x,factacc)
        % 	subplot(2,2,2)
        % 	plot(x,dlhabs)
        % 	subplot(2,2,3)
        % 	plot(x,facthc .* factlc)
        % 	subplot(2,2,4)
        % 	plot(x,exp(-(kinetic_factor .* landau - (npar0 + dlhabs ./ 2 + upshift) ) .^ 2  ./ dlhabs .^ 2))
        %     drawnow
    otherwise
        
        if lobe < 0
            dn0 = pi .* dn0;
        end
        effacc = min(1,exp((npar0 -  acc(:,end))./dn0));
        
        % par default l'absorption a lieu au minimum de la courbe de landau
        %dhc  = abs(landau - hc);
        dhc  = landau;
        %mask = (landau == (min(landau,[],2) * ve));
        mask = (dhc == (min(dhc,[],2) * ve));
        indmat = ones(size(temps)) * (1:length(x));
        indlh = max(1,round(sum(mask .* indmat,2) ./ sum(mask,2)));
        vali  =  zeros(size(temps));
        
        % condition de landau + effet de largeur de spectre
        % l'onde ne peut pas se coupler sous acc mais se propage
        % recherche du n// de premiere absorbtion
        landau2 = landau;
        landau2((landau < acc) | (landau < lc) | (landau > hc)) = Inf;
        dd  = abs(landau2 - npar0 - upshift);
        [ddmin,indlhs] = min(dd,[],2);
        % indice  ok
        indok = find(isfinite(ddmin));
        indlh(indok) = indlhs(indok);
        if length(temps) > 5
            try
                indlh = round(medfilt1(indlh,3));
            catch
                fprintf('T');
                indlh  = round((cat(1,indlh(1),indlh(1:end-1)) + indlh + cat(1,indlh(2:end),indlh(end))) ./ 3);
            end
        end
        vali(indok) = 1;
        valimat = vali * ve;
        % nouvelle valeur en sortie
        xlh  = x(indlh)';
        
        nparabs  = diag(landau(:,indlh));
        nparwave = npar0 + upshift;
        nparu  = diag(nparwave(:,indlh));
        accabs  = diag(acc(:,indlh));
        lcabs   = diag(lc(:,indlh));
        indnonok = find(~isfinite(nparabs));
        indok = find(isfinite(nparabs));
        if ~isempty(indnonok)
            nparabs(indnonok) = mean(nparabs(indnonok));
        end
        nparabs = nparabs * ve;
        nparu   = nparu   * ve;
        
        % correction imprecision en 21 points
        nparabs = max(nparu,nparabs);
        
        % largeur spatiale ?
        % gaussienne d'absorbsion
        % gaussienne ecart a l'optimum . la largeur est change par la propagation de l'onde (dispersion)
        % l'onde est evanescente en dehors du domaine de propagation, elle est de preference absorbee dans le domaine de propagation
        % dn0 represente le decalage du a la non linearite de l'absorption
        fact   = 1 + 0.5 .* max(0,(nparabs - nparu) ./ dn0);
        % wave + space;
        dlhabs = max(dn0 .* fact + width ./ Raxe ,((~vali) .* npar0) * ve);
        maskxmax  = (xlh * ve) > (vt * x);
        pabs        = exp(-0.5 .* (nparabs - nparu + npar0 + upshift - landau) .^ 2  ./ dlhabs .^ 2);
        pabs_u      = exp(-0.5 .* (nparabs - nparu + npar0 + upshift - min(hc,landau)) .^ 2  ./ dlhabs .^ 2);
        % la puissance decroit rapidement apres le maximum si la caustique haute est en dessous de l'absorbtion landau
        indi = find(vali);
        pabs(indi,:)   = pabs(indi,:) .* (maskxmax(indi,:) == 0) + pabs_u(indi,:) .* (maskxmax(indi,:) ~= 0);
        pabs   = pabs ./ max(eps,(max(pabs,[],2) * ve));
        
end

% forme par defaut
fplh    = exp(-(vt*x -(xlhin * ve)).^ 2 ./ ((max(dlhin,0.05) * ve) .^ 2) ./ 2);
%disp([mean(xlhin),mean(dlhin)])


% repli si pas valide sur la forme par defaut
%indni = find(~vali);
indni = find(any(~isfinite(pabs),2));
pabs(indni,:) = fplh(indni,:);

% normalisation
pabs   = pabs ./ max(eps,(max(pabs,[],2) * ve));


%%%%%%%%%%%%%%%%%%
% debut efficacite
%%%%%%%%%%%%%%%%%%%
switch upshiftmode
    case  {'newmodel','newmodel + tail','newmodel retune'}
        
        % d'apres les simulations completes ALOHA/C3PO/LUKE
        w1     = min(1 ./ landau,1 ./ npar0 - eps);
        %vth    = phys.c .* sqrt(1 - (1 ./ (phys.e .* tep ./ phys.me ./ phys.c .^ 2 + 1)) .^ 2);
        %nth    = 1 ./ sqrt(1 - (1 ./ (phys.e .* tep ./ phys.me ./ phys.c .^ 2 + 1)) .^ 2);
        %nth    = phys.c ./ (sqrt(2 .* phys.e .* tep ./ phys.me));
        w2     = 1 ./ npar0;
        % lineaire /quasi lineaire (cf. publication LUKE)
        Dpar   = 0.32 .* ((plh  ./ 1e6) * ve) ./ max(width,rmx) ./ Raxe .* sqrt(tep ./ 1e3) ./ (nep ./ 1e19) .^ (3/2) ./ ...
            sqrt(S) .*  landau .^ 2 ./ dlhabs;
        quasi  = 1 + 0.5 .* (1 + tanh(log(max(eps,Dpar)) - log(0.1)));
        %figure(21);clf;plot(x,sum((plh * ve).* quasi,1)./ max(1,sum((plh * ve),1)));drawnow
        % facteur de calibration (a fixer)
        % * 0.5  = ?
        % /2 = quasi
        % *4 = effet zeff (4/(5+zeff))
        % 124e20 ./ lnl = ref : V. P. Ridolfini et al , NF 45 (2005) 1386-1395
        f1 = 31e20 ./ lnl;
        % sans la dependance en zeff
        % cf Decker/ Peysson On Self consitent Simulation of the lower Hybrid Current Drive
        efficiency = (effacc * ve)  .* max(0,f1 .* quasi .* (w2 .^ 2 - w1 .^ 2)./ max(eps,log(w2 ./ w1)));
        jlh        = pabs .* efficiency ./ max(1,nep);
        %figure(21);clf;plot(temps,efficiency);drawnow
        % l'efficacite ne doit pas tomber a 0 sinon la generation de courant ne s'engage pas
        efficiency = trapz(x,pabs .* efficiency .* (vt*x),2) ./ max(eps,trapz(x,pabs .* (vt*x),2));
        
        % grandeurs derivees
        nparu  = trapz(x,abs(vt*x) .* pabs .^ 2 .* landau,2) ./ max(eps,trapz(x,abs(vt*x) .* pabs .^ 2,2));
        dlhabs   = trapz(x,abs(vt*x) .* pabs .* dlhabs,2) ./ max(eps,trapz(x,abs(vt*x) .* pabs,2));
        mask     = pabs == (max(pabs,[],2) * ve);
        nparabs   = sum(landau .* mask,2) ./ max(1, sum(mask,2));
        
    otherwise
        w1     = 1 ./ min(hc,max(npar0,nparabs));
        w2multi     = 1 ./ max(acc,lc);
        w2single    = 1 ./ npar0;
        %rap         = erf(max(0,landau - npar0 - upshift) ./ npar0);
        rap         = min(1,max(0,log(fact ./ 2 + (upshift ./ npar0) - 0.5 + dn0 ./ npar0)));
        w2          = w2multi .* rap + w2single .* (1 - rap);
        % lineaire /quasi lineaire (cf. publication LUKE)
        Dpar   = 0.32 .* ((plh  ./ 1e6) * ve) ./ max(width,rmx) ./ Raxe .* sqrt(tep ./ 1e3) ./ (nep ./ 1e19) .^ (3/2) ./ ...
            sqrt(S) .*  nparabs .^ 2 ./ dlhabs;
        quasi  = 4 + 2 .* (1 + tanh(log(max(eps,Dpar)) - log(0.1)));
        % facteur de calibration
        f1 = 0.5;
        % sans la dependance en zeff
        % cf Decker/ Peysson On Self consitent Simulation f the lower Hybrid Current Drive
        efficiency = max(0,f1 .* quasi .* (w2 .^ 2 - w1 .^ 2)./ max(eps,log(max(w2,w1) ./ w1))) .* 1e20 ./ max(1,fact - 0.15 .* npar0 ./ dn0);
        % l'efficacite ne doit pas tomber a 0 sinon la generation de courant ne s'engage pas
        %efficiency = max(eps,vali) .*  effacc .* trapz(x,pabs .* efficiency .* (vt*x),2) ./ max(eps,trapz(x,pabs .* (vt*x),2));
        efficiency = effacc .* trapz(x,pabs .* efficiency .* (vt*x),2) ./ max(eps,trapz(x,pabs .* (vt*x),2));
        jlh        = 0 .* pabs;
        
end

%%%%%%%%%%%%
% securite
%%%%%%%%%%%%
% repliment
xmem   = x;
pabsmem = pabs;
% securite
indpbad = find(sum(pabs,2) == 0);
pabs(indpbad,:)  = fplh(indpbad,:);


% caclul de elements 0d
mask    = pabs == (max(pabs,[],2) * ve);
xlh     = sum((vt * x) .* mask,2) ./ max(1, sum(mask,2));
% la valeur de dlh en sortie est differente de cette definition
dlh      = max(0.05,sqrt(trapz(x,abs(vt*x) .* pabs .* (vt * x - xlh * ve) .^ 2,2) ./  max(eps,trapz(x,abs(vt*x) .* pabs,2))));
% anti oscillation ?
if transitoire == 1
	fpout    = pabs + sqrt(-1) .* jlh;
else
        fpout    = vt *  mean(pabs,1);
	xlh      = vt .* mean(xlh);
	dlh      = vt .* mean(dlh);
end

% securite faible puissance
indp =find(plh < 1e3);
if ~isempty(indp)
	fpout(indp,:) = 1;
end

% pas de plot
if plotonoff == 0
	return
end

indp = find(plh >1e3);
if isempty(indp)
       return
end

mref =mean(pabs(indp,:),1);
mfp  = mean(fplh(indp,:),1);
% 
if lobe < 0
	h = findobj(0,'type','figure','tag','z0acclhneg');
	if isempty(h)
       		h=figure('tag','z0acclhneg');
	else
       		figure(h);
	end
else
	h = findobj(0,'type','figure','tag','z0acclhpos');
	if isempty(h)
       		h=figure('tag','z0acclhpos');
	else
       		figure(h);
	end
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1)
hold on
plot(-1,-1,'b',-1,-1,'b',-1,-1,'r',-1,-1,'g',-1,-1,'om',-1,-1,'*k',-1,-1,'c');
legend('lower caustic','upper caustic','accessibility','Landau absorbtion','absorbtion point','launch point','LH wave');
plot(xmem,lc(indp,:),'b',xmem,hc(indp,:),'b',xmem,acc(indp,:),'r', ...
     xmem,landau(indp,:),'g',xlh(indp),nparabs(indp),'om',1,npar0,'*k',xmem,npar0 + upshift(indp,:),'c');
if lobe < 0
	title('n// < 0')
else
	title('n// > 0')
end
axis([0,1,0,10])
xlabel('r/a')
ylabel('n_{//}')
%legend('lower caustic','upper caustic','accessibility','Landau absorption','absorption point','launch point');
subplot(2,2,2)
plot(x,pabs(indp,:),'b',x,mref./max(mref),'r')
xlabel('r/a')
ylabel('power deposition (normalized, r=averaged)')
subplot(2,2,3)
if length(indp) == 1
	plot(temps(indp),xlh(indp),'*r',temps(indp),dlh(indp),'ob');
else
	plot(temps(indp),xlh(indp),'r',temps(indp),dlh(indp),'b');
end
legend('xlh','dlh')
set(gca,'ylim',[0,1])
xlabel('time (s)')
ylabel('maximum and width deposition (r/a unit)')
subplot(2,2,4)
if length(indp) == 1
 	plot(temps(indp),nparabs(indp,1),'ob',temps(indp),nparu(indp,1),'+r',temps(indp),dlhabs(indp,1),'xg');
else
 	plot(temps(indp),nparabs(indp,1),'b',temps(indp),nparu(indp,1),'r',temps(indp),dlhabs(indp,1),'g');
end
legend('n_{// landau absorbed}','n_{// wave absorbed} ','absorbtion width')
set(gca,'ylim',[0,10])
xlabel('time (s)')
ylabel('n_{//}')
drawnow

% no graph with slider
if plotonoff < 2
	return
end

% 
if lobe < 0
	h = findobj(0,'type','figure','tag','z0acclhneg_slider');
	if isempty(h)
       		h=figure('tag','z0acclhneg');
	else
       		figure(h);
	end
else
	h = findobj(0,'type','figure','tag','z0acclhpos_slider');
	if isempty(h)
       		h=figure('tag','z0acclhpos');
	else
       		figure(h);
	end
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1)
hold on
plot(-1,-1,'b',-1,-1,'b',-1,-1,'r',-1,-1,'g',-1,-1,'om',-1,-1,'*k',-1,-1,'c');
legend('lower caustic','upper caustic','accessibility','Landau absorbtion','absorbtion point','launch point','LH wave');
zplotprof(gca,temps,xmem,lc,'color','b');
zplotprof(gca,temps,xmem,hc,'color','b');
zplotprof(gca,temps,xmem,acc,'color','r');
zplotprof(gca,temps,xmem,landau,'color','g');
zplotprof(gca,temps,cat(2,xlh,xlh),cat(2,nparabs,nparabs),'color','m','linestyle','none','marker','*');
zplotprof(gca,temps,cat(2,xmem,xmem),cat(2,npar0 + upshift,npar0 + upshift),'color','c');

%plot(xmem,lc(indp,:),'b',xmem,hc(indp,:),'b',xmem,acc(indp,:),'r', ...
%     xmem,landau(indp,:),'g',xlh(indp),nparabs(indp),'om',1,npar0,'*k',xmem,npar0 + upshift(indp,:),'c');
if lobe < 0
	title('n// < 0')
else
	title('n// > 0')
end
axis([0,1,0,10])
xlabel('r/a')
ylabel('n_{//}')
%legend('lower caustic','upper caustic','accessibility','Landau absorption','absorption point','launch point');
subplot(2,2,2)
%plot(x,pabs(indp,:),'b',x,mref./max(mref),'r')
plot(x,mref./max(mref),'r');
hold on
zplotprof(gca,temps,x,pabs,'color','b');
xlabel('r/a')
ylabel('power deposition (normalized, r=averaged)')
legend('time averaged deposition','deposition');

subplot(2,2,3)
if length(indp) == 1
	plot(temps(indp),xlh(indp),'*r',temps(indp),dlh(indp),'ob');
else
	plot(temps(indp),xlh(indp),'r',temps(indp),dlh(indp),'b');
end
legend('xlh','dlh')
addsslide(gca,min(temps),max(temps));
set(gca,'ylim',[0,1])
xlabel('time (s)')
ylabel('maximum and width deposition (r/a unit)')
subplot(2,2,4)
if length(indp) == 1
 	plot(temps(indp),nparabs(indp,1),'ob',temps(indp),nparu(indp,1),'+r',temps(indp),dlhabs(indp,1),'xg');
else
 	plot(temps(indp),nparabs(indp,1),'b',temps(indp),nparu(indp,1),'r',temps(indp),dlhabs(indp,1),'g');
end
legend('n_{// landau absorbed}','n_{// wave absorbed} ','absorbtion width')
addsslide(gca,min(temps),max(temps));
set(gca,'ylim',[0,10])
xlabel('time (s)')
ylabel('n_{//}')
drawnow

% add cursor on each line
function addsslide(ha,tmin,tmax)
hl = findobj(ha,'type','line');
for k=1:length(hl)
    time = get(hl(k),'xdata');
    indok = find((time >= tmin) & (time <= tmax));
    if ~isempty(indok)
        time  = time (indok);
        data = get(hl(k),'ydata');
        data = data(indok);
        cl   = get(hl(k),'color');
        if length(time) > 1
            zplotprof(gca,time(:),cat(2,time(:),time(:)),cat(2,data(:),data(:)),'color',cl,'linestyle','none','marker','s');
        end
    end
end