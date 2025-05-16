% calcul du porfil de densite et des donnees derivees 
function [nep,n0,S0,n0m,S0m,D,V,Vware,n0a,spellet,ge,frac_pellet,sn0fr,taup] = z0profne(cons,zs,profli,option,geo)

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

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);



% calcul du profil de ne si non disponible ou si diffsuion non active
if isfield(profli,'vpr')
	x = profli.xli;
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(zs.nem));
	vpr = profli.vpr;
	spr = profli.spr;
else
	x   = linspace(0,1,21);
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(zs.nem));
	spr   = (2 .* zs.sp) * x;
	vpr   = (2 .* zs.vp) * x;
end

% le profil de ne est calculer ici
ne0   = max(zs.nebord + 1e13,zs.nem .* (1 + zs.ane));
%alpha = max(eps-1,((ne0 - zs.nebord) ./ max(1e13,zs.nem - zs.nebord)) * ve - 1); 
alpha = max(eps-1,((ne0 - zs.nebord) ./ max(0.1 .* zs.nem,zs.nem - zs.nebord)) * ve - 1); 
nep   = ((ne0 - zs.nebord) * ve)  .* (vt * ux)  .^ alpha  + zs.nebord * ve;
nep(:,end) = zs.nebord;
%nep   = ((ne0 - zs.nebord) * ve)  .* (vt * ux)  .^ (max(0,(zs.ane .* zs.nem) ./ max(1e13,zs.nem - zs.nebord)) * ve) + zs.nebord * ve;

%figure(21);subplot(2,1,1);plot(cons.temps,alpha);subplot(2,1,2);plot(cons.temps,ne0, cons.temps,zs.nem, cons.temps, zs.nebord);legend('ne0','nem','nebord');drawnow

% cas mode H
switch option.ne_free
    case 0
        vol_tep = cumtrapz(x,vpr,2);
        r_tep = (geo.a * ve) .* sqrt(vol_tep ./ (max(vol_tep,[],2) * ve));
        if isfield(profli,'qjli')
            q_tep = profli.qjli;
        else
            q_tep   = z0qp(x,zs.q0,zs.qa);
        end
        s_tep = r_tep ./ q_tep .* pdederive(x,q_tep,0,2,2,1) ./ pdederive(x,r_tep,2,2,2,1);
        nout = 1 - 4 ./ 3 ./ (geo.R * ve) .* cumtrapz(x,pdederive(x,r_tep,2,2,2,1) .* (s_tep + 3/8),2);
        % first renormalisation without taking into account nebord
        nem_tep = trapz(x,vpr.* nout,2) ./ trapz(x,vpr,2);
        fact_tep = zs.nem  ./ nem_tep;
        nout     = nout .* (fact_tep * ve);
        % minimal value of nout
        nout(:,end) = zs.nebord;
        nout = max(nout,zs.nebord * ve);
        % final renormalisation 
        nem_tep = max(1,trapz(x,vpr.* (nout - zs.nebord * ve),2) ./ trapz(x,vpr,2));
        fact_tep = max(0,zs.nem - zs.nebord) ./nem_tep;
        nout     = (nout - zs.nebord * ve).* (fact_tep * ve) + zs.nebord * ve;
       
    case 4
        if isfield(profli,'ptot')
            [nout,fbest] = formene_4d(profli.xli,profli.vpr,zs.nem,zs.ane,zs.nebord,zs.negr,option.neped_expo);
        else
            [nout,fbest] = formene_4d(x,vpr,zs.nem,zs.ane,zs.nebord,zs.negr,option.neped_expo);
        end
    otherwise
        if isfield(profli,'ptot')
            [nout,fbest] = formene(profli.xli,profli.vpr,zs.nem,zs.ane,zs.nebord,profli.tep .* profli.nep);
        else
            [nout,fbest] = formene(x,vpr,zs.nem,zs.ane,zs.nebord);
        end
end


switch option.ne_shape
case 'Lmode'
  % rien
case 'Hmode'
  nep    = nout;
  ne0    = nout(:,1);
otherwise
  indh         = find(zs.modeh ~= 0);
  nep(indh,:)  = nout(indh,:);
  ne0(indh)    = nout(indh,1);
end

if isfield(profli,'tip') & (option.ane == 5)
	nep  = max(1,profli.tip - profli.tip(:,end) * ve);
	nemx = trapz(x,vpr .* nep,2)./ zs.vp;
	rap  = max(0,zs.nem -zs.nebord) ./max(eps,nemx);
	nep  = (rap * ve) .* nep + zs.nebord * ve;
	%figure(181);clf;plot(x,nep);drawnow	
end

% cas de donnees experimentales
if isappdata(0,'NE_EXP') & isfield(zs,'temps');
	nes = getappdata(0,'NE_EXP');
	nex = max(1e13,interp1_ex(nes.temps,nes.ne,zs.temps,'nearest','extrap'));
	nep = pchip(nes.x,nex,x);
	nemx = trapz(x,vpr .* nep,2)./ zs.vp;
	rap  = zs.nem ./max(1e13,nemx);
	nep  = (rap * ve) .* nep;	
end

n0         = ones(size(nep));
S0         = ones(size(nep));
n0m         = ones(size(nep));
S0m         = ones(size(nep));
spellet     = zeros(size(nep));
V          = zeros(size(nep));
n0a        = vt;
D          = ones(size(nep));
Vware      = zeros(size(nep));
ge         = ones(size(nep));
frac_pellet = 0 .* vt;
sn0fr      = 0 .* vt;
taup       = zs.taue;


if ~isfield(profli,'nip') | ~isfield(profli,'tip')	
	return
end
% securite sur la forme de ne pour garantit le signe du gradient de Te
%  if isappdata(0,'NE_EXP') & isfield(zs,'temps');
%  	%rien
%  elseif isfield(profli,'ptot') 
%  
%  	indbad = find(profli.tep(:,end - 2) <= (1.1 .* profli.tep(:,end - 1) - 0.1 .* profli.tep(:,end)));
%  	if ~isempty(indbad)
%                  tep   = profli.tep; 
%                  nep_mem = nep;
%  		pep   = tep .* nep;
%  		dne =  0.01 .* (nep(:,end-1) - nep(:,end));
%   		nl    = trapz(x,nep-nep(:,end) * ve,2);
%  		nb    = 101;
%                  racx  = 3;
%                  vxl   = cat(2,1:(length(x)-racx-1),length(x)-1,length(x));
%  	        while(~isempty(indbad)) & (nb > 0)
%   		        nenew  = nep(indbad,:);
%                          nenew(:,end-1) = max(min(nenew(:,end -1) + dne(indbad),nenew(:,end - racx)),nenew(:,end));  
%                          nenew  = pchip(x(vxl),nenew(:,vxl),x);                
%  		        nlnew  = max(1,trapz(x,nenew-nep(indbad,end) *ve,2));
%  		        nenew  = (nenew - nenew(:,end) * ve) .* ((nl(indbad) ./ nlnew) * ve) + nenew(:,end) * ve;		
%  			nep(indbad,:) = nenew;
%  			tep           = pep ./ nep;	
%                          indbad = find(tep(:,end - 2) <= (1.1 .* tep(:,end - 1) - 0.1 .* tep(:,end)));		
%                          nb  = nb - 1;
%  		end
%  
%  	end	
%  
%  end

% energie ponderee pour l'injection de neutres
einj_w = ((option.einj .* max(1,real(zs.pnbi)) + option.einj2 .* imag(zs.pnbi)) ./ max(1, real(zs.pnbi)+ imag(zs.pnbi))) * ve;

%  % calcul de l'effet des itb sur le profil de ne
%  % hypothese pour le transport
if isappdata(0,'NE_EXP') & isfield(zs,'temps');
    %rien
elseif isfield(profli,'ge')
    
    if any(zs.hitb >1)
        inditb = find(zs.hitb >1);
        gest =  max(0,profli.ge);
        gest(:,1) = 0;
        ne_std = - cumtrapz(profli.xli(end:-1:1),gest(:,end:-1:1) ./ profli.xieshape(:,end:-1:1)     .* option.itb_density,2);
        ne_itb = - cumtrapz(profli.xli(end:-1:1),gest(:,end:-1:1) ./ profli.xieshape_itb(:,end:-1:1) .* option.itb_density,2);
        ne_std = ne_std(:,end:-1:1);
        ne_itb = ne_itb(:,end:-1:1);
        % meme valeur au bord
        ne_itb = ne_itb + (ne_std(:,end) - ne_itb(:,end-1)) * ve;
        % on autorise uniquement l'augmentation du piquage du a l'injection de neutre
        dne  = ne_itb - ne_std;
        %% dne must be >=0
        %%dne(dne < 0) = 0;
        % pas d'interference avec le piedestal
        dne(:,end-1:end) = 0;
        % conservation du nombre de particules
        offset = trapz(x,dne .* profli.vpr,2) ./ max(eps,trapz(x,profli.vpr .* (dne ~= 0),2));
        %figure(21);clf;subplot(2,2,1);plot(offset(inditb));
        dne(:,1:end-2) = dne(:,1:end-2)  - offset * ve(1:end-2);
        %hold on; plot(trapz(x,dne .* profli.vpr,2),'r');
        rapd = min(dne ./ (nep - 0.5 .* (nep(:,end-1) + nep(:,end-2)) * ve + eps),[],2);
        fact = - min(-1,rapd) * ve;
        %subplot(2,2,2);plot(fact(inditb));
        dne  = dne ./ fact;
        indbad = find(any(dne>=nep,2));
        dne(indbad,:) = 0;
        %subplot(2,2,3);hold on
        %plot(profli.xli,dne(inditb,:),'b');
        nep(inditb,1:(end - 2)) = nep(inditb,1:(end - 2)) + dne(inditb,1:(end -2));
        %subplot(2,2,4);hold on
        %plot(profli.xli,nep(inditb,:),'r');
        %drawnow
        if any(~isfinite(nep(:))) || any(imag(nep(:)))
            disp('Bad values in electron density profile');
            nep = real(nep);
            nep(~isfinite(nep)) = 1e13;
        end
    end
end

% element utiles
nip   = profli.nip;
n1p   = profli.n1p;
tep   = profli.tep;	
tip   = profli.tip;
zeffp = profli.zeff;

% cas du glacon
% cas donnees par l'utilisateur
nep_mem =nep;
frac_pellet = vt .* max(0,option.pif);
spellet = vt * exp(- (x - option.pix) .^ 2 ./ 2 ./ max(eps,option.piw) .^ 2);
spellet = spellet ./ max(eps,trapz(x,spellet .* profli.vpr,2) * ve);
if option.pif ~= 0
    if (option.pif == 1)
        % detection de l'injection de glacon
        dndt = z0dxdt(cons.nbar,cons.temps)./ cons.nbar;
        sig_pellet = max(0,dndt - medfilt1(dndt,3));
        % le pellet doit partir au temps precedent l'effet sur le plasma
        sig_pellet(1:end-1) = sig_pellet(2:end);
        sig_pellet(1) = 0;
        [tv,frac_pellet] = z0ode(zs.temps,sig_pellet,zs.taue,0);
        frac_pellet = max(0,min(1,frac_pellet));
        frac_pellet = frac_pellet .* (frac_pellet >= 0.05);
        
    end
    if option.piw == 0
        % utilisation du model ngs pour cacluler le depot
        mp = trapz(x,max(1e13,nep) .* profli.vpr,2);
        rp = max(1e-8,1.58e-10 .* mp .^ (1/3));
        %mp = 2 .* phys.avo ./ 20e-6 .* (4 ./ 3 .* pi .* rp .^3);
        % estimation de la vitesse (scaling NGS)
        fsc = 0.079;
        vp = ((1 - option.pix) ./ fsc .* (tep(:,1) ./ 1e3) .^ (5/9) .* (nep(:,1) ./ 1e20) .^ (1/9) ./ (mp ./1e20) .^(5/27)) .^ 3;
        vp = max(1e-3,vp);
        al = 1;
        de = 1.13;
        drdx0 = 3.1e-14 .* (al ^ 2 ./ de .^ (2/3)) .* nep .^(1/3) .* tep .^(5/3) ./ (vp *ve);
        spellet(:,end) = rp;
        for k = length(x):-1:2
            spellet(:,k - 1) = max(0,spellet(:,k)  + (vt * (x(k - 1) - x(k))) .*  drdx0(:,k) ./ (max(1e-9,spellet(:,k)) .^ (2/3))) .* (spellet(:,k) >= 1e-9);
        end
        spellet = 2 .* phys.avo ./ 20e-6 .* (4 ./ 3 .* pi .* spellet .^3);
        spellet = max(0,pdederive(x,spellet,0,2,2,1) ./ max(eps,profli.vpr));
        spellet(:,1) = spellet(:,2);
        spellet = spellet ./ (trapz(x,spellet .* profli.vpr,2) * ve);
    end
	%
    % modification du profil de densite
    %
    netot   = trapz(x,nep .* profli.vpr,2);
    % efficacite de fuelling
    %efficacite          = ones(size(x));
    %efficacite(end)     = 0;
    %efficacite(end - 1) = 1/3;
    %efficacite(end - 2) = 2/3;
    efficacite          = vt * (1 - x) ;
    efficacite          = efficacite .* ((trapz(x,profli.vpr,2) ./ trapz(x,profli.vpr .* efficacite,2)) * ve);
    % profile d'increment de densite brute
    nspellet0 = max(0,spellet .* max(1e-3,zs.taue * ve)  .* efficacite);
    % normalisation
    nspellet0 = nspellet0 ./ max(1,trapz(x,nspellet0 .* profli.vpr,2) * ve);
    npellet0 = ((netot .* frac_pellet) * ve) .* nspellet0;
    % effet continue de l'injection
    npellet1 = 2 .* cumtrapz(x,npellet0,2);
    npellet1 = npellet1(:,end) * ve - npellet1;
    % ajout de l'effet transitoire
    % l'effet continu sur la densite est que differentie
    npellet2 = max(0,z0dxdt(max(1e-3,zs.taue * ve) .* npellet0 .* efficacite,cons.temps));
    npellet2 = ((nep(:,1) ./ max(nep(:,1) ./ 2,max(npellet2,[],2)) ./ 2) * ve) .* npellet2;
    % ajout de l'effet transitoire
    npellet = npellet1 + npellet2;
    % recombinaison
    %npellet = ((nep(:,1) ./ max(nep(:,1),max(npellet,[],2))) * ve) .* npellet;
    %npellet = min(nep(:,1) * ve,npellet);
    ntot_pellet = trapz(x,npellet .* profli.vpr,2);
    ntot_core   = trapz(x,(nep - nep(:,end) * ve) .* profli.vpr,2);
    nep   = (nep - nep(:,end) * ve) .* (max(0,(ntot_core - ntot_pellet) ./ max(1e13,ntot_core)) * ve) + ...
        npellet + nep(:,end) * ve;
    
    %  	figure(51)
    %  	plot(x,npellet0,'b',x,npellet1,'r',x,npellet2,'m',x,npellet,'c',x,nep,'g')
    %  	drawnow
    %  	%keyboard
end


% effect du Tungstene (ionisation dans le plasma)
if (option.W_effect == 1)  && isfield(profli,'nzp')
    if option.Sn_fraction > 0
        zwloc  = z0wavez(profli.tep);
        nzwloc = (1 - option.Sn_fraction) .* profli.nwp .* zwloc;
        dnep_w = nzwloc - nzwloc(:,end) * ve;
        zsnloc  = z0snavez(profli.tep);
        nzsnloc = option.Sn_fraction .* profli.nwp .* zsnloc;
        dnep_sn = nzsnloc - nzsnloc(:,end) * ve;
        nep    = nep + dnep_w + dnep_sn;     
    else
        zwloc  = z0wavez(profli.tep);
        nzwloc = profli.nwp .* zwloc;
        dnep_w = nzwloc - nzwloc(:,end) * ve;
        nep    = nep + dnep_w;
    end
end



% renormalisation finale de nep sur nbar (precision augmentee)
indh = find((zs.modeh > 0.5) &  ...
    (max(nep(:,1:(end - 1)) -nep(:,end - 1) * ve(1:(end - 1)),[],2) > 1e13) & ...
    ((zs.nbar - nep(:,end-1)) > 1e13));
if ~isempty(indh)
    nbarcorenew       = max(1,trapz(profli.xli(1:(end - 1)),nep(indh,1:(end - 1)) -nep(indh,end - 1) * ve(1:(end - 1)) ,2));
    nbarcore          = zs.nbar(indh) - nep(indh,end-1);
    rapcore           = nbarcore ./ nbarcorenew;
    indok             = find((rapcore < 10) & (rapcore > 0.1));
    indh              = indh(indok);
    rapcore 	  = rapcore(indok);
    if ~isempty(indh)
        nepmem = nep;
        nep(indh,1:(end - 1)) = (nep(indh,1:(end - 1))  -  nep(indh,end - 1) * ve(1:(end - 1))) .*  ...
            (rapcore * ve(1:(end - 1))) + nep(indh,end - 1) * ve(1:(end - 1));
        indbad = indh(find(any(nep(indh,:) < 1e13,2)));
        if ~isempty(indbad)
            nep(indh,:) = nepmem(indh,:);
        end
    end
end

fact = (zs.nbar - nep(:,end)) ./ max(1e13,trapz(profli.xli,nep - nep(:,end) * ve,2));
fact = min(10,max(0.1,fact));
nep  = (nep - nep(:,end) * ve) .* (fact * ve) + nep(:,end) * ve;

% cas des DDS
indice_inv_min = 1;
if any(zs.indice_inv > indice_inv_min) & (option.qdds < 0)
    % nouveau profil de densite
    dnedx                  = cat(2,0 * vt,diff(nep - nep(:,end) * ve,1,2));
    mask                   = (vt * (1:length(x))) > (zs.indice_inv * ve);
    dnedx                  = dnedx .* mask;
    nenew                  = cumsum(dnedx,2);
    nenew                  = nenew - nenew(:,end) * ve + nep(:,end) *ve ;
    % normalisation (conservation de la matiere)
    dnetot                 = trapz(x,(nenew - nep) .* profli.vpr,2);
    maskc                   = (vt * (1:length(x))) <= (zs.indice_inv * ve);
    maskc(:,1:2)           = 1;
    dvp                    = trapz(x,maskc .* profli.vpr,2);
    dn                     = min(0,dnetot ./ max(eps,dvp));
    %nep_mem                = nep;
    nep                    = nenew - (dn * ve) .*  (~mask);
    
    %  	    ind_plot_ = find(zs.indice_inv > indice_inv_min);
    %  	    figure(22)
    %  	    plot(x,nep(ind_plot_,:),'r',x,nep_mem(ind_plot_,:),'b')
    %  	    drawnow
    
end


if any(~isfinite(nep(:))) | any(imag(nep(:)))
	%disp('Z0PROFNE : ne profile corrected of imag and NaN')
    fprintf('>');
	nep(~isfinite(nep)) = 1e13;
	nep = real(nep);
elseif any(nep(:) < 0)
  	%disp('Z0PROFNE : ne profile corrected of negative value')
    fprintf('<');
	nep(nep <= 0) = 1e13;
end

%    keyboard
%  

% profils etendus
tex = cat(2,tep,max(0.1,zs.telim));
tix = cat(2,tip,tex(:,end));
zeffx = cat(2,zeffp,zeffp(:,end));
switch option.modeh
case 0
	nex = cat(2,nep,nep(:,end).* tex(:,end) ./tex(:,(end-1)));
	nix = cat(2,nip,nip(:,end).* tex(:,end) ./tex(:,(end-1)));	
	n1x = cat(2,n1p,n1p(:,end).* tex(:,end) ./tex(:,(end-1)));	
otherwise
	nex = cat(2,nep,max(1e13,zs.nelim));
	nix = cat(2,nip,nip(:,end).* nex(:,end) ./nex(:,(end-1)));
	n1x = cat(2,n1p,n1p(:,end).* nex(:,end) ./nex(:,(end-1)));
end
xx   = cat(2,x, 2 .*x(end) - x(end-1));
vex  = ones(size(xx));

rhomax = pchip(x,profli.rmx,xx(end));
grho   = pchip(x,profli.grho,xx);
grho2  = pchip(x,profli.grho2,xx);
vpr    = pchip(x,profli.vpr,xx);
vpr_eq = vpr ./ (rhomax * vex);
		

% la geometrie passe de slab au bord a torique pres du centre 
%(quand le perimetre est proche de la distance a la source externe)
% le neutre venant des molecules sont traites en slab
% ceux issu de l'echange de charge en torique

% constantes
ee  = 1.602176462e-19;
ua  = 1.66053873e-27;
 
% calcul des sections efficaces 
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(tex,tix);

%% equilibre entre 1s et 2s pour le neutres de centre 
alphas = nex .* sss ./ ( nex .* sss + Ass);
% etat d'equilibre 1s/2s
sviss  = svi2s .* alphas + svi1s .* (1 - alphas);

% les vitesses de penetration des neutres
% v0  = sqrt(2 .* ee .* max(13.6,tix) ./ ua ./ (zs.meff * vex));
% utilisation de la vitesse du son (plus stable numeriquement, generalement peu differente)
v0  = sqrt(ee .* (tix + tex .* zeffx)./ ua ./ (zs.meff * vex));
v0th = sqrt(2 .* ee .* max(0.1,(zs.telim*vex)) ./ ua ./ (zs.meff * vex));

% terme de diffusion des neutres (cf. these Laporte)
% ref K. Beckurts et K. Wirtz , Neutron Physics, Springer-Verlay (1964)
lpm  = v0 ./ (nex .* sviss + nix .* sii + n1x .* svcx);
D0   = max(1,lpm .* v0 ./ 3);

% terme de rcombinaison
sa = nex .* n1x .* svrec;

% coeffiecient de transport
rap   = max(0.1,min(10,mean(profli.xii ./ profli.xie,2)))*ve;
tid1  = pdederive(x,tip,0,2,2,1);
nid1  = pdederive(x,nip,0,2,2,1);
ted1  = pdederive(x,tep,0,2,2,1);
ned1  = pdederive(x,nep,0,2,2,1);


%  % essayer le modele de Weiland article de Tockar	
%  % model de Backer (eta = 1/sqrt(2));
%  lite = pdederive(x,tep,0,2,2,1) ./ max(13.6,tep) ./(rhomax * ve);
%  lite(:,1) = lite(:,2);
%  ud = 2 .* pi .* profli.qjli ./ profli.fdia ./ profli.r2i;
%  vd = (1 - 2./3 .* option.ct )  .* option.cq .* (profli.qjli >= option.qdds) .* ...
%       pdederive(x,ud,0,2,2,1) ./ ud ./ (rhomax * ve) + ...
%       option.ct .* lite;
%       
%  % lien si diffusion de la fonction de distribution     
%  D  = profli.xie  ./ max(eps,3 ./ 2  - vd ./ min(-eps,lite));   
%  V  = vd .* D;
%  
%  % securite
%  D  = max(1e-1,min(30,D));
%  V  = max(-10,min(30,V)); 


	 
% pinch de Ware (attention au signe au moment d'ecrire le flux)
lnei   = 15.2 - 0.5.*log(profli.nep ./ 1e20) + log(profli.tep ./ 1e3); % pour H, D et T uniquement
taue   =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me)  .* ...
	                   (phys.e .* profli.tep) .^ (3/2) ./ profli.nip ./ lnei;
% ware : Hiton HAzeltine 76 + these Laporte
nue    = 1./ min(taue,1);
vthe   = sqrt(2 .* phys.e .* profli.tep  ./ phys.me);
nues   = nue .* profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthe;
k1     = (0.53 + profli.zeff ) ./ profli.zeff ./ (1 + 1.32 .* profli.zeff);
g      = 1 + sqrt((0.44 + 0.11 .* profli.zeff) .* nues)    + (0.55 + 0.225 .* profli.zeff) .* nues;
fzeff  = (1 + k1  - k1 ./ g .* profli.ftrap) ./ g; 
emax   = max(1e-3,abs(10 .* zs.vloop ./ 2 ./ pi ./ profli.Raxe(:,end))) * ve;
epar   = min(emax,max(-emax,profli.epar));
Vware = fzeff .* epar ./ max(eps,profli.bpol) .* max(0,min(1,profli.ftrap));
%Vware(:,1) = Vware(:,2);
Vware(:,1) = 0;
Vware  = max(-30,min(30,Vware)); 

% calcul de taup par inversion de la formule donnant le nebord
% Wesson 9.3.8
% vitesse du son
switch option.gaz
    case 1
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1);
    case 2
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2);
    case {3,5}
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
            (2 + 3.* cons.iso) .* (1 + cons.iso));
    case 4
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4);
    case 11
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
                            (1 + 11 .* cons.iso) .* (1 + cons.iso));
end
cs   = cs0 .* sqrt(profli.tep(:,end));  
ntot = trapz(profli.xli,profli.vpr .* profli.nep,2);
if isfield(zs,'dsol')
	dsol = zs.dsol;
elseif option.sol_lscale  == 0
    dsol = profli.Raxe(:,end) .* profli.epsi(:,end) ./ 100;
elseif option.sol_lscale  > 0
    dsol = profli.Raxe(:,end) .* profli.epsi(:,end) .* option.sol_lscale;
else
    dsol = - profli.Raxe(:,end) .* option.sol_lscale;
end
%dsol = profli.Raxe(:,end) .* profli.epsi(:,end) ./ 100;
% compatible avec Stangeby 4.31 car pour les rapport d'aspect que l'on considere : R*sin(ut) ~ a.
% dsol est largeur caracteristique de la SOL : lambda_n = dsol.
% la formule du Wesson sous estime le temps de confinement de la matiere
% celle de Stangeby semble plus proche de l'experience (il y a un facteur 4 entre les 2).
% Wesson : taup = ntot ./ zs.nebord ./ profli.spr(:,end) ./ dsol ./ cs;
% METIS initiale : taup = ntot ./ zs.nebord ./ (4 .* pi * a) ./ dsol ./ cs;
% Stangeby : taup =  2 .* ntot ./ zs.nebord ./ (4 .* pi * R * sin(ut)) ./ dsol ./ cs;
taup = option.ftaup .* ntot ./ zs.nebord ./ profli.spr(:,end) ./ dsol ./ cs;

% estamation de na0
fpe     = min(1,taup ./ max(taup,zs.taue)) .* frac_pellet;
%n0ge = cumtrapz(profli.xli,((z0dxdt(profli.rmx(:,end),zs.temps) ./ profli.rmx(:,end)) * profli.xli) .* ...
%       pdederive(profli.xli,nep .* profli.vpr,0,2,2,1),2) ./ max(eps,profli.vpr_tor .* profli.grho2);
n0ge = cumtrapz(profli.xli,((z0dxdt(profli.rmx(:,end),zs.temps) ./ profli.rmx(:,end)) * profli.xli) .* ...
       pdederive(profli.xli,nep .* profli.vpr,0,2,2,1),2);
n0a    = max(1e13,zs.nem) .* max(1e-3,zs.vp) ./ max(1e-3,taup) +  ...
	 max(0,z0dxdt(max(1e13,zs.nem) .* max(1e-3,zs.vp),cons.temps))  - n0ge(:,end);
n0a    = max(1,n0a - trapz(x,profli.pnbi ./ (einj_w .* ee) .* profli.vpr,2)) .* max(0.1,1 - fpe);
% fraction qui retourne dans le plasma sous forme de neutres
switch option.configuration
case {0,1}
    n0a    = option.fn0a .* n0a;
case {2,3}
    n0a    = (zs.xpoint .* option.fn0a_div + (~zs.xpoint) .* option.fn0a) .* n0a;
otherwise
    n0a    = option.fn0a_div .* n0a;
end
% n0a(1) = n0a(2);
n0a    = max(1,n0a);
if all(zs.n0a == 0)
	zs.n0a(:) = n0a;
end

% external data for neutral source
if isappdata(0,'NEUTRAL_EXP')
    neutral_exp = getappdata(0,'NEUTRAL_EXP');
    n0a = import_external_n0a(neutral_exp.temps(:),neutral_exp.n0a(:),n0a,zs.temps);
end

% calcul de la source moleculaire (probleme slab)
% a la dsmf 50% en 2s pour les neutres issu des molecule
n0m = int_n0(xx,rhomax,vpr_eq,grho,nex,v0th,(svi1s + svi2s) ./ 2 + (svcx .* (n1x ./ nex) + sii .* (nix ./ nex)) ,0,zs.n0a);
% la source d'eletrons associee (uniquement les termes de ionisation fournissent des electrons)
S0m = max(1,((svi1s + svi2s) ./ 2 + sii .* (nix ./ nex)) .* nex .* n0m);
% normalisation
ind_in = find(xx <= 1);
nnn    = zs.n0a ./ trapz(xx(ind_in),S0m(:,ind_in) .* vpr_eq(:,ind_in),2) ./ rhomax;
S0m    = (nnn * vex) .* S0m;
n0m    = (nnn * vex) .* n0m;
% securite la source de neutres est moins dense que le plasma
n0m = min(n0m,nex);
% la source associee
S0m = max(1,((svi1s + svi2s) ./ 2 + sii .* (nix ./ nex)) .* nex .* n0m);
% normalisation
nnn    = zs.n0a ./ trapz(xx(ind_in),S0m(:,ind_in) .* vpr_eq(:,ind_in),2) ./ rhomax;
S0m    = (nnn * vex) .* S0m;
n0m    = (nnn * vex) .* n0m;

% convergence entre neutres chauds et froids
nnn_mem = nnn;
n0m = 0.5 .* n0m;
for lk = 1:31
    % calcul de la nouvelle valeur de n0 (echange de charge)
    [n0,n0d1] = difloc(xx,vpr_eq .* grho2 .* D0,vpr_eq .* (rhomax * vex) .^ 2 , ...
	        -(sviss + sii.* (nix ./ nex)) .* nex, sa + svcx .* n1x .* n0m);
    %[n0_,n0d1_] = difloc2(xx,vpr_eq .* grho2 .* D0,vpr_eq .* (rhomax * vex) .^ 2 , ...
	%	    -(sviss + sii.* (nix ./ nex)) .* nex, sa + svcx .* n1x .* n0m);
    % eventually slower than original version.
    %disp([max(max(abs(n0_ - n0))),max(max(abs(n0d1_ - n0d1)))]);
    n0 = max(0,n0);


    %  calcul dela source d'eletrons associee :
    % uniquement les termes de ionisation fournissent des electrons, 
    % le terme de recombinaison est un terme de perte
    % le terme d'echange de charge n'intervient pas dans la source puisqu'il ne fournit pas d'electron
    S0      = (sviss + sii.* (nix ./ nex)) .* nex .* n0  - sa;
    % securite en cas de non convergence du solver
    indbad = find((n0(:,1)>= max(n0,[],2)) & (n0(:,1) >= (max(n0m,[],2) / 3)));
    n0(indbad,:) = 0;
    S0(indbad,:) = 0;

   % normalisation finale
    s0tot   = S0m(:,1:(end-1)) + S0(:,1:(end-1));
    n0stot  = trapz(x,s0tot .* profli.vpr,2);
    s0alt   = S0m(:,1:(end-1));
    n0salt  = trapz(x,s0alt .* profli.vpr,2);
    indbad  = find(n0stot <= 1);
    s0tot(indbad,:) = s0alt(indbad,:);
    S0(indbad,:) = 0;
    %nnn     = zs.n0a ./ max(1,trapz(x,s0tot .* profli.vpr,2));
    nnn     = n0a ./ max(1,trapz(x,s0tot .* profli.vpr,2));
    s0tot   = (nnn * ve)  .* s0tot;
    S0      = (nnn * vex) .* S0;
    S0m     = (nnn * vex) .* S0m;
    n0m     = (nnn * vex) .* n0m;
    n0      = (nnn * vex) .* n0;
  
    %err =sqrt(sum((nnn - nnn_mem).^2) ./ sum(nnn .^ 2))
    err =max(abs(nnn - nnn_mem) ./ max(eps,nnn));
    if err < 1e-6
        break
    end
end

% prise en compte de l'injection de glacon
if option.pif > 0
    fpe     = min(1,taup ./ max(taup,zs.taue)) .* frac_pellet;
    ns0tot  = trapz(x,s0tot .* profli.vpr,2) ./ max(0.1,1 - fpe);
    
    spellet = ((ns0tot .* fpe) * ve) .* spellet;
    s0tot   = s0tot + spellet;
else
    spellet  =  0 .* (vt * ve);
end

% calul de la source totale (neutre de bord +nbi)
snbi  = profli.pnbi ./ (einj_w .* ee);
stot  = max(1, s0tot + snbi);

% calcul du flux
ge =   1 ./ max(eps,profli.vpr_tor .* profli.grho2) .*  ...
       cumtrapz(profli.xli, profli.vpr .* stot - z0dxdt(nep .* profli.vpr,zs.temps) +   ...
      ((z0dxdt(profli.rmx(:,end),zs.temps) ./ profli.rmx(:,end)) * profli.xli) .* ...
       pdederive(profli.xli,nep .* profli.vpr,0,2,2,1),2);
ge(:,1) = 0;

% calcul de Dn et Vn
if (option.kishape == 0) && strcmp(option.density_model,'default');
        % la variation radial ne fonctionne que sur certains chocs JET
	D = (vt * linspace(1,0.3,length(profli.xli))) .* (profli.xie .* profli.xii) ./ max(1e-3,profli.xie + profli.xii);
 	D(D > 30) = 30;
	D(D < 1e-3) = 1e-3;
	V = - (ge + D .* ned1 ./ (profli.rmx(:,end) * ve)) ./ nep;
	% separation du pinch neoclassique et anormale
	V  = V - Vware;
else
	switch option.density_model
        case 'curvature'

		s   = pdederive(vt*profli.xli, profli.qjli,0,2,2,1) ./ profli.qjli .* (vt*profli.xli);
		VsD = 2 .*  (0.25  + 2 ./ 3 .* max(0,s)) .* profli.ri;
                %D donne par le pinch de courbure + nep (valeur en l'absence de source)
                D = abs((ge + Vware .* nep) ./ (ned1 ./ (profli.rmx(:,end) * ve) + VsD .* nep));
                 % calcul de V a partir du flux
                V = - (ge + D .* ned1 ./ (profli.rmx(:,end) * ve)) ./ nep;
		% separation du pinch neoclassique et anormale
		V  = V - Vware;

        case 'control'

                warning off
                D = - (Vware .* nep + ge) ./ ned1  .* (profli.rmx(:,end) * ve);
		warning on
                D_backup = abs(profli.xie .* profli.xii) ./ max(1e-3,profli.xie + profli.xii) ./ profli.qjli;
                D(~isfinite(D)) = D_backup(~isfinite(D));
                Doff = max(0.1 .* mean(D(:)),0.001);
                D = D - min(0,min(D - Doff,[],2) * ve);
                  % calcul de V a partir du flux
                V = - (ge + D .* ned1 ./ (profli.rmx(:,end) * ve)) ./ nep ;
		% separation du pinch neoclassique et anormale
		V  = V - Vware;

        case 'minconv'

                coef_min = 1e-3;
                % from BgBs model
                Dbgbs    = max(coef_min,profli.xie) .* max(coef_min,profli.xii) ./ max(2 .* coef_min,profli.xie + profli.xii) ./ profli.qjli;
                Dbgbs(zs.modeh ~= 0,end) = coef_min;
                Dbgbs(zs.modeh ~= 0,end - 1) = coef_min;
                Dbgbs(zs.modeh ~= 0,end - 2) = 0.5 .* Dbgbs(zs.modeh ~= 0,end - 2) + 0.5 .* coef_min;

                warning off
                D = - (Vware .* nep + ge ) ./ ned1 .* (profli.rmx(:,end) * ve);
 		warning on
                D(~isfinite(D)) = sqrt(eps);
		D = max(Dbgbs,D);
                % calcul de V a partir du flux
                V = - (ge + D .* ned1 ./ (profli.rmx(:,end) * ve)) ./ nep;
		% separation du pinch neoclassique et anormale
		V  = V - Vware;
                % signe > 0
		V  = max(0,V);
                % calcul final de D
                warning off
                D = - ((V + Vware) .* nep + ge ) ./ ned1 .* (profli.rmx(:,end) * ve);
 		warning on
                D(~isfinite(D)) = sqrt(eps);

	otherwise

		D = max(1e-3,min(profli.xie,profli.xii));
		D(D > 30) = 30;
		D(D < 1e-3) = 1e-3;
		V = - (ge + D .* ned1 ./ (profli.rmx(:,end) * ve)) ./ nep;
		% separation du pinch neoclassique et anormale
		V  = V - Vware;
	end
end


if 1>2
	warning off
	figure(117);
	clf
	subplot(2,3,1)
	plot(x,ge,'b');
	ylabel('ge');
	subplot(2,3,2)
	plot(x,nep,'b');
	ylabel('ne (b)');
	subplot(2,3,3)
	semilogy(xx,D0,'b')
	ylabel('D0');
	subplot(2,3,4)
	plot(x,D,'r',x,V,'b',x,Vware,'g');
	ylabel('D (r), V (b)  &  Vware (g)');
	subplot(2,3,5)
	semilogy(xx,S0,'r',xx,S0m,'m',x,snbi,'g',x,stot,'b');
	ylabel('S0 (r), S0m (m), Snbi (g) & Stot (b)');
	subplot(2,3,6)
	semilogy(xx,n0,'r',x,profli.n0,'m',xx,n0m,'b',x,profli.n0m,'c');
	ylabel('n0 (r), n0 old (m), n0m (b), n0m old (c)');
	drawnow
	warning on
	%save contexte_profne0d_loc
	keyboard
end


% resultat pour les temps qui ont converges
% mise a dimension
n0  = n0(:,1:(end-1));
S0  = S0(:,1:(end-1));
n0m = n0m(:,1:(end-1));
S0m = S0m(:,1:(end-1));


% =============================
% source de freinage au bord
% =============================
if ~isfield(profli,'omega')
	return
end

% palsma de fond
switch option.gaz
case {1,11}
   zj = 1;
   aj = 1;
case {2,5}
   zj = 1;
   aj = 2;
case 3
   zj = 1;
   aj = mean((2  + 3 .* cons.iso)  ./  (1+ cons.iso));
case 4
   zj = 2;
   aj = 4;
end

% impurete principale
zimp = option.zimp;
aimp = ceil(zimp .* (7/3));

% 2ieme impurete
zmax = option.zmax;
amax = ceil(zmax .* (7/3));


% pour chaque espece d'ions
nDp   = max(1e13,profli.n1p .* ((zs.nDm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nTp   = max(1e13,profli.n1p .* ((zs.nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nHp   = max(1e13,profli.n1p - nTp - nDp);
nhep  = max(1e13,profli.nhep);
nz1p  = max(1e13,profli.nzp);
nz2p  = max(1e11,profli.nzp .* option.rimp);   
% masse
Mtor    = phys.mp .*  (nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p);
% unite sn0fr :  kg  m^2 s^-2 = N.m
% l'echange de charge ne porte que sur un electron quel que soit la charge de l'ion
fact    = max(1,profli.nep ./ profli.nip);
sn0fr   =  - trapz(x, profli.vpr .* fact .* svcx(:,1:(end-1)) .* profli.n0m .*  ...
                          Mtor .* profli.Raxe .^ 2 .* profli.omega,2);

%  figure(21)
%  plot(x,profli.nip .* profli.n0m .* svcx(:,1:(end-1)),'r',x,stot,'b')
%  drawnow

% integration de la source de neutres
% version matricielle
function [n0,S0] = int_n0(x,rhomax,vpr,grho,ne,v0,s0,Sa,ct)


% vectorisaition
ve     = ones(size(x));
vt     = ones(size(rhomax));
rhomax = rhomax * ve;
vpr(:,1) = vpr(:,2);
ct     = ct * ve;

% la derivee de vpr pose probleme au centre
int1 =  - rhomax .* cumtrapz(x, pdederive(x,vpr,0,2,2,1) ./ vpr ./ rhomax  + ...
                                pdederive(x,grho,0,2,2,1) ./ grho ./ rhomax  + ...
                                pdederive(x,v0,0,2,2,1) ./ v0 ./ rhomax  - ...
			        s0 .* ne ./ grho ./ v0,2);
int1(:,1) = int1(:,2);
int1      = int1 - int1(:,end) * ve;
eint1 = exp(int1);
eint1(~isfinite(eint1)) = 0;

int2  = rhomax .* cumtrapz(x,(pdederive(x,vpr .* grho .* v0,0,2,2,1) ./ rhomax - ...
                            s0 .* ne .* vpr) ./ (vpr .* grho .* v0),2);
int2(:,1) = int2(:,2);
int2      = int2 - int2(:,end) * ve;

eint2     = exp(int2);
eint2(~isfinite(eint2)) = 0;

inter3    = Sa .* eint2 ./ grho ./ v0;
inter3(~isfinite(inter3)) = 0;

int3      = - rhomax .* cumtrapz(x,inter3,2);
int3(:,1) = int3(:,2);
int3      = int3 - int3(:,end) * ve;
			    
n0  = max(1,(ct + int3) .* eint1);	
n0(~isfinite(n0)) = 0;

S0  = max(1,pdederive(x,vpr .* grho .* v0 .* n0 ,0,2,2,1) ./ rhomax ./ vpr);
S0(~isfinite(S0)) = 1;

% correction pb derivee au centre
S0  = max(1,pchip(cat(2,0,x(4:end)),cat(2,ones(size(S0,1),1),S0(:,4:end)),x));
S0(~isfinite(S0)) = 1;


function [n0,n0d1] = difloc(x,A,C,S0,S1)


n0   = zeros(size(A));
n0d1 = zeros(size(A));
B  = pdederive(x,A,0,2,2,1);
A(:,1) = A(:,2) ./ 10;
pas = -1;
debut = length(x);
fin   = 2;
%
n0d1(:,end) = - C(:,end) .* S1(:,end) ./ B(:,end);
%
for k= debut:pas:fin
        xk = x(k);
	xkm1 = x(k+pas);
	dx   = (xk - xkm1);
	indk = [k,k+pas];
	
	[nd1,ud1] = attn0(xk,n0d1(:,k),n0(:,k),x(indk),C(:,indk),S0(:,indk),S1(:,indk),A(:,indk),B(:,indk));
	nk1       = dx .* nd1;
	uk1       = dx .* ud1;
	
	[nd1,ud1] = attn0(xk + dx ./ 2,n0d1(:,k) + uk1 ./ 2,n0(:,k) +  nk1 ./ 2, ...
	            x(indk),C(:,indk),S0(:,indk),S1(:,indk),A(:,indk),B(:,indk));
	nk2       = dx .* nd1;
	uk2       = dx .* ud1;
	
	[nd1,ud1] = attn0(xk + dx ./ 2,n0d1(:,k) + uk2 ./ 2,n0(:,k) +  nk2 ./ 2, ...
	            x(indk),C(:,indk),S0(:,indk),S1(:,indk),A(:,indk),B(:,indk));
	nk3       = dx .* nd1;
	uk3       = dx .* ud1;
	
	[nd1,ud1] = attn0(xkm1,n0d1(:,k) + uk3,n0(:,k) +  nk3, ...
	            x(indk),C(:,indk),S0(:,indk),S1(:,indk),A(:,indk),B(:,indk));
	nk4       = dx .* nd1;
	uk4       = dx .* ud1;	
		
	n0(:,k+pas)     = n0(:,k)   + (nk1 ./ 6 + nk2 ./ 3 + nk3 ./ 3 + nk4 ./ 6);
	n0d1(:,k+pas)   = n0d1(:,k) + (uk1 ./ 6 + uk2 ./ 3 + uk3 ./ 3 + uk4 ./ 6);
		
end

function [nd1,ud1]= attn0(x,u,n,xx,C,S0,S1,A,B)



% interpolation
C  = interpn0(xx,C,x);
S0 = interpn0(xx,S0,x);
S1 = interpn0(xx,S1,x);
A  = interpn0(xx,A,x);
B  = interpn0(xx,B,x);


% evaluation
nd1 = u;
ud1 = -(C .* (S0 .* n + S1) + B .* u) ./ A;
%ud1 = cos(x.*pi);


function y=interpn0(xx,yy,x)

y = yy(:,1) + (yy(:,2) - yy(:,1) ) ./  ...
    (xx(:,2) - xx(:,1)) .* (x - xx(:,1));


function [n0,n0d1] = difloc2(x,A,C,S0,S1)


n0   = zeros(size(A));
n0d1 = zeros(size(A));
B  = pdederive(x,A,0,2,2,1);
A(:,1) = A(:,2) ./ 10;
pas = -1;
debut = length(x);
fin   = 2;
%
n0d1(:,end) = - C(:,end) .* S1(:,end) ./ B(:,end);
%
for k= debut:pas:fin
    xk = x(k);
    xkm1 = x(k+pas);
    dx   = (xk - xkm1);
    indk = [k,k+pas];
    
    % compute one time interpolant
    Ci  = interp_n0_fact(x(indk),C(:,indk));
    S0i = interp_n0_fact(x(indk),S0(:,indk));
    S1i = interp_n0_fact(x(indk),S1(:,indk));
    Ai  = interp_n0_fact(x(indk),A(:,indk));
    Bi  = interp_n0_fact(x(indk),B(:,indk));
    
    [nd1,ud1] = attn0_fast(xk,n0d1(:,k),n0(:,k),x(indk),Ci,S0i,S1i,Ai,Bi);
    nk1       = dx .* nd1;
    uk1       = dx .* ud1;
    
    [nd1,ud1] = attn0_fast(xk + dx ./ 2,n0d1(:,k) + uk1 ./ 2,n0(:,k) +  nk1 ./ 2, ...
        x(indk),Ci,S0i,S1i,Ai,Bi);
    nk2       = dx .* nd1;
    uk2       = dx .* ud1;
    
    [nd1,ud1] = attn0_fast(xk + dx ./ 2,n0d1(:,k) + uk2 ./ 2,n0(:,k) +  nk2 ./ 2, ...
        x(indk),Ci,S0i,S1i,Ai,Bi);
    nk3       = dx .* nd1;
    uk3       = dx .* ud1;
    
    [nd1,ud1] = attn0_fast(xkm1,n0d1(:,k) + uk3,n0(:,k) +  nk3, ...
        x(indk),Ci,S0i,S1i,Ai,Bi);
    nk4       = dx .* nd1;
    uk4       = dx .* ud1;
    
    n0(:,k+pas)     = n0(:,k)   + (nk1 ./ 6 + nk2 ./ 3 + nk3 ./ 3 + nk4 ./ 6);
    n0d1(:,k+pas)   = n0d1(:,k) + (uk1 ./ 6 + uk2 ./ 3 + uk3 ./ 3 + uk4 ./ 6);
    
end

function [nd1,ud1]= attn0_fast(x,u,n,xx,C,S0,S1,A,B)



% interpolation
C  = interpn0_val(xx,C,x);
S0 = interpn0_val(xx,S0,x);
S1 = interpn0_val(xx,S1,x);
A  = interpn0_val(xx,A,x);
B  = interpn0_val(xx,B,x);


% evaluation
nd1 = u;
ud1 = -(C .* (S0 .* n + S1) + B .* u) ./ A;
%ud1 = cos(x.*pi);


function y=interpn0_val(xx,in,x)

y = in.offset + in.slope .* (x - xx(:,1));


function out = interp_n0_fact(xx,yy)

out.offset = yy(:,1);
out.slope  = (yy(:,2) - yy(:,1) ) ./  ...
    (xx(:,2) - xx(:,1));


function val_out = import_external_n0a(time_in,val_in,val_origin,time_out)

val_out          = interp1_ex(time_in,val_in,time_out,'linear','extrap');
indnok           = find(~isfinite(val_out));
val_out(indnok)  = val_origin(indnok);
val_out          = max(1,val_out);


