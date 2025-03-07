function [tite,tauee,tauii,pei,tauei,profli]= ....
         zeqie0prof(vp,a,k,te0,nem,pel,pion,wth,ae,ane,ate,A,temps,tite_in,pei_in,...
	 nDm,nTm,n1m,nhem,nimpm,zimp,rimp,zmax,pped,tebord,nebord,xiioxie,grad_ped, ...
     hitb,qdds,indice_inv,profli,fact_confinement,coef_shape,creux,exp_shape,kishape, ...
     Sn_fraction,telim,extended_qei,frhe0,nboronm)
	
% frhe0 is 0 if option.gaz ~= 5
% nboronm is 0 if option.gaz ~= 11

% backward compatibility
if nargin < 34
  coef_shape = 'bgb';
end
if nargin < 35
  creux = 0;
end
if nargin < 36
  exp_shape = 0;
end

% constantes
gg   = 0.64e15;
ee   = 0.1602176462e-18;
me   = 9.10938188e-31;
mp   = 1.6726485e-27;
c    =   2.99792458e8;
mu0  =   4*pi*1e-7;
epsi0 = 1./c.^2./mu0;

% excurtion autour de 1
rr  = 10;
%telim = 1.2e5;

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - zimp);
mask = (dd == min(dd));
Aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(Aimp)
    Aimp = 7/3 .* zimp;
end
dd   = abs(Z_el - zmax);
mask = (dd == min(dd));
Amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(Amax)
    Amax = 7/3 .* zmax;
end

% minimal temperature
if nargin >= 33
    t_sub = min(13.6,(1 - fact_confinement) .* tebord + fact_confinement .* 13.6);
end

% control du transport dans le piedestal
if imag(xiioxie)
    xiioxie_ped = imag(xiioxie);
    xiioxie     = real(xiioxie);
else
    xiioxie_ped = [];
end

% initialisation
profli.xie  = ones(size(te0,1),size(profli.xli,2));
profli.xii  = ones(size(te0,1),size(profli.xli,2));
profli.qe   = zeros(size(te0,1),size(profli.xli,2));
profli.qi   = zeros(size(te0,1),size(profli.xli,2));
profli.qei  = zeros(size(te0,1),size(profli.xli,2));

% cas de donnees experimentales
if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') && isfield(profli,'qjli')

	% profil utiles
	x   = profli.xli;
	vpr = profli.vpr;
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(te0));
	nep = profli.nep;
	nip = profli.nip;
	tep = profli.tep;
	tip = profli.tip;

	% profil de densites
    if frhe0 > 0 % He4 is account separtly if option.gaz == 5
        nhep3 = profli.nhep;
        nhep = profli.nhep;
    else
        nhep3 = zeros(size(profli.nhep));
        nhep  = profli.nhep;        
    end
	nHp  = profli.n1p .* ((max(0,n1m - nDm - nTm + nboronm) ./ n1m) * ve); % no tritium but boron if option.gaz == 11
	nDp  = profli.n1p .* ((nDm ./ n1m) * ve);
	nTp  = profli.n1p .* (((nTm - nboronm) ./ n1m) * ve); % no tritium but boron if option.gaz == 11
	nzp  = profli.nzp;
	if isfield(profli,'nwp')
		nwp  = profli.nwp;
	else
		nwp  = 0 .* nzp;		
    end
    nbp  = profli.n1p .* ((nboronm ./ n1m) * ve); % no tritium but boron if option.gaz == 11
	
	% element de geometrie complementaire
        if isfield(profli,'grho2')
		grho2 = profli.grho2;
	else
		grho2 = ones(size(vpr));	
	end
	rhomax = profli.rmx(:,end);


	tes = getappdata(0,'TE_EXP');
	tis = getappdata(0,'TI_EXP');
	tex = max(1,interp1_ex(tes.temps,tes.te,temps,'nearest','extrap'));
	tix = max(1,interp1_ex(tis.temps,tis.ti,temps,'nearest','extrap'));
	tep = pchip(tes.x,tex,x);
	tip = pchip(tis.x,tix,x);

	% recopie
	if exp_shape == 0
	  profli.tep = max(1,real(tep));
	  profli.tip = max(1,real(tip));
	else
	    tep = max(1,real(tep));
	    tip = max(1,real(tip));
	    tea = tep(:,end);
	    tia = tip(:,end);
	    wth_exp  = max(eps,1.602176462e-19 .* (3/2) .* trapz(x,(tep .* profli.nep + tip .* profli.nip) .* vpr,2));
	    wth_0    = max(eps,1.602176462e-19 .* (3/2) .* trapz(x,((tea * ve) .* profli.nep + (tia * ve).* profli.nip) .* vpr,2));
            warning off
            fact_exp =(wth - wth_0) ./ (wth_exp - wth_0);
            warning on
            fact_exp(~isfinite(fact_exp)) = 1;
            fact_exp(fact_exp <= 0) = 1;
            fact_exp = fact_exp * ve;
            %%figure(21);clf;plot(fact_exp);drawnow
            profli.tep = min(telim,max(1,fact_exp  .* (tep - tea *ve) + tea *ve)); 
            profli.tip = min(telim,max(1,fact_exp  .* (tip - tia *ve) + tia *ve)); 
	end

	
	% terme source
    % equipartition
    switch  extended_qei
        case 'on'
            [qeib0,qeib] = equipartition_full(tep,tip,nep,nHp,nDp,nTp,nhep,nzp,rimp .* nzp,(1 - Sn_fraction) .* nwp, Sn_fraction .* nwp,zimp,zmax,nbp,nhep3);
        otherwise
            warning off
            lnei          =  15.2 - 0.5 .* log(nep ./ 1e20) + log(tep ./ 1e3);
            warning on
            ind = find(~isfinite(lnei) | (lnei <10));
            if ~isempty(ind)
                lnei(ind) = 10 .* ones(1,length(ind));
            end
            taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (epsi0 .^ 2 ./ ee .^ 4 .* sqrt(me))  .* ...
                ((ee .* max(t_sub * ve,tep)) .^ (3/2) ./ nep ./ lnei);
            if Sn_fraction > 0
                factor_w_sn = (1 - Sn_fraction) .* z0wavez(tep) .^ 2 ./ 183.84 .* nwp + ...
                    Sn_fraction.* z0snavez(tep) .^ 2 ./ 118.71 .* nwp;
                qeib0    = 3 .* me ./ mp ./ taues .* (nHp + nDp ./ 2 + nTp ./ 3 + nhep + nhep3 .* (4/3.02)  +  ...
                    (zimp .^ 2 ./ Aimp) .* nzp + rimp .*  (zmax .^ 2 ./ Amax) .* nzp + factor_w_sn + 25/11 .* nbp);
                
            else
                qeib0    = 3 .* me ./ mp ./ taues .* (nHp + nDp ./ 2 + nTp ./ 3 + nhep + nhep3  .* (4/3.02)  +  ...
                    (zimp .^ 2./ Aimp) .* nzp + rimp .*  (zmax .^ 2 ./ Amax) .* nzp + z0wavez(tep) .^ 2 ./ 183.84 .* nwp + 25/11 .* nbp);
            end
            qeib     = qeib0 .* ee .* (tep - tip);
    end
	% les flux de puissance
	Qe = max(0,cumtrapz(x,profli.source_el .* vpr,2));
	Qi = max(0,cumtrapz(x,profli.source_ion .* vpr,2));
	Qei = cumtrapz(x,qeib .* vpr,2);

	% calcul des coefficients de transport associes
	tepd1       = pdederive(x,tep,0,2,2,1);
	tipd1       = pdederive(x,tip,0,2,2,1);
	% attention vpr = V' * rhomax !
	xiea         = - (Qe - Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tepd1)) ./ ...
	                max(1e13,nep) ./ grho2 .* (rhomax * ve) .^ 2;
	xiia         = - (Qi  + Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tipd1)) ./  ...
	                max(1e13,nip) ./ grho2.* (rhomax * ve) .^ 2;
	%  prolongation
	xiea         = pchip(x(2:end),real(xiea(:,2:end)),x);
	xiia         = pchip(x(2:end),real(xiia(:,2:end)),x);
	
	% limitation
	profli.xie  = max(1e-3,min(30,xiea));
	profli.xii  = max(1e-3,min(30,xiia));
	
	% flux de chaleur en sortie :
	profli.qe   = 	(Qe - Qei);
	profli.qi   = 	(Qi + Qei);
	profli.qei   = 	Qei;

	% calcul de tite
	tite = max(1./rr,min(rr,trapz(x,tip .* vpr,2) ./ trapz(x,tep .* vpr,2)));

	% pei
	pei  = Qei(:,end);

	% temps de confinement
	Qmin = max(1,1e-3 .* (abs(Qe(:,end)) + abs(Qei(:,end)) + abs(Qi(:,end))) ./ 3);
	tauee = ee .* trapz(x,vpr .* nep .* tep,2) ./ max(Qmin,Qe(:,end) - Qei(:,end));
	tauii = ee .* trapz(x,vpr .* nip .* tip,2) ./ max(Qmin,Qi(:,end) + Qei(:,end));
	tauei = ee .* trapz(x,vpr .* nip .* tip,2) ./ max(Qmin,abs(pei));

	wel  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
	wion = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);


% suite du clacul du profil de ti et te complet
elseif isfield(profli,'source_el')

	% profil utiles
	x   = profli.xli;
	vpr = profli.vpr;
	ux  = (1  - x .^ 2);
	ve  = ones(size(x));
	vt  = ones(size(te0));
	nep = profli.nep;
	nip = profli.nip;
	tep = profli.tep;
	tip = profli.tip;

	% profil de densites
    if frhe0 > 0 % He4 is account separtly if option.gaz == 5
        nhep3 = profli.nhep;
        nhep = profli.nhep;
    else
        nhep3 = zeros(size(profli.nhep));
        nhep  = profli.nhep;
    end
	nHp  = profli.n1p .* ((max(0,n1m - nDm - nTm + nboronm) ./ n1m) * ve); % no tritium but boron if option.gaz == 11
	nDp  = profli.n1p .* ((nDm ./ n1m) * ve);
	nTp  = profli.n1p .* (((nTm - nboronm)./ n1m) * ve); % no tritium but boron if option.gaz == 11
	nzp  = profli.nzp;
	if isfield(profli,'nwp')
		nwp  = profli.nwp;
	else
		nwp  = 0 .* nzp;		
	end
    nbp  = profli.n1p .* ((nboronm ./ n1m) * ve); % no tritium but boron if option.gaz == 11
        
    % forme BgB
    %if all(profli.xieshape(:,end) == 1)
    if kishape == 0
        opt_bgb = 1;
    else
        opt_bgb = 0;
    end

	% element de geometrie complementaire
        if isfield(profli,'grho2')
		grho2 = profli.grho2;
	else
		grho2 = ones(size(vpr));	
	end
	rhomax = profli.rmx(:,end);

	% hypothese pour le transport
    % forme BgB
    if opt_bgb
        switch coef_shape
            case {'cdbm','cdbm+neo'}
                [kie_bgb,kie_itb_bgb,kii_bgb,kii_itb_bgb] = z0cdbm(A,profli);
            case {'stiff','stiff+neo','stiff_limited','stiff_limited+neo'}
                [kie_bgb,kie_itb_bgb,kii_bgb,kii_itb_bgb] = z0stiff(A,profli);
            case {'alpha','alpha+neo'}
                [kie_bgb,kie_itb_bgb,kii_bgb,kii_itb_bgb] = z0alpha_limit(A,profli);
            otherwise
                [kie_bgb,kie_itb_bgb,kii_bgb,kii_itb_bgb] = z0bgbs(A,profli);
        end
        switch coef_shape
            case {'bgb+neo','cdbm+neo','stiff+neo','alpha+neo','stiff_limited+neo'}
                chii_neo    = neo_Hinton(profli,A);
                kie_bgb     = kie_bgb + chii_neo;
                kie_itb_bgb = kie_itb_bgb + chii_neo;
                kii_bgb     = kii_bgb + chii_neo;
                kii_itb_bgb = kii_itb_bgb + chii_neo;
        end
        xie     = min(10 .* kie_bgb, max(0.1 .* kie_bgb,kie_bgb .* ((hitb <= 1) * ve) + kie_itb_bgb .* ((hitb > 1) * ve)));
        %xii_bgb = min(10 .* kii_bgb, max(0.1 .* kii_bgb,kii_bgb .* ((hitb <= 1) * ve) + kii_itb_bgb .* ((hitb > 1) * ve)));
		if any(indice_inv > 1) && (qdds < 0)
			xie_st = min(10 .* profli.xieshape, max(0.1 .* profli.xieshape,profli.xieshape .*  ...
            			((hitb <= 1) * ve) + profli.xieshape_itb .* ((hitb > 1) * ve)));
		   % modification of xie due to ST
		   for ln = 1:length(temps)
		        if indice_inv(ln) > 1.5
                    indice_modif = max(fix(indice_inv(ln)),max(find(profli.qjli(ln,:) <= abs(qdds))));
					xie(ln,1:indice_modif) = xie(ln,1:indice_modif) + xie_st(ln,1:indice_modif) - xie_st(ln,indice_modif);
				else
				    indice_modif = max(find(profli.qjli(ln,:) <= abs(qdds)));
					if ~isempty(indice_modif)
						xie(ln,1:indice_modif) = xie(ln,1:indice_modif) + xie_st(ln,1:indice_modif) - xie_st(ln,indice_modif);
					end
				end
		   end
		end
        xie = max(0.1,min(30,xie));
    else
        xie     = min(10 .* profli.xieshape, max(0.1 .* profli.xieshape,profli.xieshape .*  ...
            ((hitb <= 1) * ve) + profli.xieshape_itb .* ((hitb > 1) * ve)));
    end
    %figure(121);plot(x,min(10 .* profli.xieshape, max(0.1 .* profli.xieshape,profli.xieshape .*  ...
    %    ((hitb <= 1) * ve) + profli.xieshape_itb .* ((hitb > 1) * ve))),'b',x,xie,'r');drawnow
%     if any(indice_inv > 1.5) && (qdds < 0)        
%         figure(25);
%         clf;
%         plot(x,xie(indice_inv > 1.5,:));
%         drawnow
%     end
	% traitement du bord
	p0               = tebord .* (nep(:,end) + nip(:,end)) .* ee;
	w0               = (3/2) .* trapz(x,vpr .* (p0*ve),2);
	% securite energie minimale
	wmin             = 3 .* nebord .* 30 .* ee .*trapz(x,vpr,2);
	wth              = max(wmin,max(wth,min(w0 + 1,w0 .* 1.1)));
	ppied            = pped * ve;
	ppied(:,end)     = p0;
	wpied            = min(0.9 .* wth,(3/2) .* trapz(x,vpr .* ppied,2));
	modehp           = (pped>1);
	indhp            = find(modehp);
	indnhp           = find(~modehp);
	
	
	% temperature ionique de bord
	titebord = max(0.3,min(3,tip(:,end) ./ tep(:,end)));
	tibord   = tebord .* titebord;

	% valeur initiale de teped et tiped
	teped            = tep(:,end-1);
	tiped            = tip(:,end-1);
	
        ppedin           = ee .* teped .* nep(:,end-1) + ee .* tiped .* nip(:,end-1);

	% les profil initiaux doivent etre solution du probleme en energie
	% pression initiale
	pep = ee .* tep .* nep;
	pip = ee .* tip .* nip;

	% valeur de bord
	if ~isempty(indnhp)
		pep(indnhp,:) = pep(indnhp,:)  + (ee .* tebord(indnhp) .* nep(indnhp,end) - pep(indnhp,end)) * ve;
		pip(indnhp,:) = pip(indnhp,:)  + (ee .* tibord(indnhp) .* nip(indnhp,end) - pip(indnhp,end)) * ve;
	end
	if ~isempty(indhp)
		pep(indhp,1:(end-1)) = pep(indhp,1:(end-1))  + (ee .* teped(indhp) .* nep(indhp,end - 1)- pep(indhp,end-1)) * ve(1:(end-1));
		pep(indhp,end) = ee .* tebord(indhp) .*  nep(indhp,end);
		pip(indhp,1:(end-1)) = pip(indhp,1:(end-1))  + (ee .* tiped(indhp).* nip(indhp,end - 1) - pip(indhp,end-1)) * ve(1:(end-1));
		pip(indhp,end) = ee .* tibord(indhp).*  nip(indhp,end);
	end
	
	if (grad_ped == 2) || (grad_ped == 3)

		% cas du mode L
		if ~isempty(indnhp)
			% temperature  electronique
			tep0 = tep - tebord * ve;
			pep0 = ee .* nep .* tep0;
	
			% temperature  ionique
			tip0 = tip - tibord * ve;
			pip0 = ee .* nip .* tip0;
	
			% calcul de la constante de normalisation
			wloc_L    = (3/2) .* trapz(x,vpr .*(pep0 + pip0),2);
			w0loc_L   = ee .*(3/2) .* trapz(x,vpr .*(nep .* (tebord * ve) + nip .* (tibord * ve)),2);
	
			fact(indnhp)  = max(0,wth(indnhp) - w0loc_L(indnhp)) ./ max(eps,wloc_L(indnhp));
	
			tep(indnhp,:) = (fact(indnhp)' * ve) .* tep0(indnhp,:) + tebord(indnhp) * ve;
			pep(indnhp,:) = ee .* nep(indnhp,:) .* tep(indnhp,:);
	
			tip(indnhp,:) = (fact(indnhp)' * ve) .* tip0(indnhp,:) + tibord(indnhp) * ve;
			pip(indnhp,:) = ee .* nip(indnhp,:) .* tip(indnhp,:);
	
		end
	
		% cas du mode H
		if ~isempty(indhp)
			% temperature  electronique
			tep0 = tep(:,1:(end-1)) - teped * ve(1:(end-1));
			pep0 = ee .* nep(:,1:(end-1)) .* tep0;
	
			% temperature  ionique
			tip0 = tip(:,1:(end-1)) - tiped * ve(1:(end-1));
			pip0 = ee .* nip(:,1:(end-1)) .* tip0;
	
			% calcul de la constante de normalisation
			wloc_H    = (3/2) .* trapz(x(1:end-1),vpr(:,1:end-1) .*(pep0 + pip0),2);
			vpx       = cumtrapz(x,vpr,2);
			p0        = ee .* nep(:,end) .* tebord + ee .* nip(:,end) .* tibord;
	
			w0loc_H   = ee .*(3/2) .* trapz(x(1:(end-1)),vpr(:,1:(end-1)) .*(nep(:,1:(end-1)) .* (teped * ve(1:(end-1))) +  ...
							nip(:,1:(end-1)) .* (tiped * ve(1:(end-1)))),2) + ...
					(3/2) .* (pped + p0)  ./ 2 .* (vpx(:,end) - vpx(:,end-1));
	
			fact(indhp)  = max(0,wth(indhp) - w0loc_H(indhp)) ./ max(eps,wloc_H(indhp));
	
			tep(indhp,1:(end-1)) = (fact(indhp)' * ve(1:(end-1))) .* tep0(indhp,:) + teped(indhp) * ve(1:(end-1));
			tep(indhp,end) =  tebord(indhp);
			pep(indhp,:) = ee .* nep(indhp,:) .* tep(indhp,:);
	
			tip(indhp,1:(end-1)) = (fact(indhp)' * ve(1:(end-1))) .* tip0(indhp,:) + tiped(indhp) * ve(1:(end-1));
			tip(indhp,end) =  tibord(indhp);
			pip(indhp,:) = ee .* nip(indhp,:) .* tip(indhp,:);
	
		end
	
	else

		% energie non normalise
		wloc = (3/2) .* trapz(x,vpr .* (pep + pip),2);
		ppedloc = pep+pip;
		ppedloc(:,1:(end-1)) =(pep(:,end-1) + pip(:,end-1)) * ve(1:(end-1));
		wpiedloc  =  (3/2) .* trapz(x,vpr .* ppedloc,2);
		w0loc  = (3/2) .* (pep(:,end) + pip(:,end)) .* trapz(x,vpr,2);
		if ~isempty(indhp)
			fact(indhp)  = (wth(indhp) - wpiedloc(indhp)) ./ max(1,wloc(indhp) - wpiedloc(indhp));
			pep(indhp,1:(end-1)) = (pep(indhp,1:(end-1)) - (pep(indhp,(end-1)) * ve(1:(end-1)) )) .*  ...
					(fact(indhp)' * ve(1:(end-1))) + (pep(indhp,end-1) * ve(1:(end-1)) );
			pip(indhp,1:(end-1)) = (pip(indhp,1:(end-1)) - (pip(indhp,(end-1)) * ve(1:(end-1)) )) .* ...
					(fact(indhp)' * ve(1:(end-1))) + (pip(indhp,end-1) * ve(1:(end-1)) );
		end
		if ~isempty(indnhp)
			fact(indnhp)  = (wth(indnhp) - w0loc(indnhp)) ./ max(1,wloc(indnhp) - w0loc(indnhp));
			pep(indnhp,:) = (pep(indnhp,:)  - (pep(indnhp,end) * ve )) .*  ...
					(fact(indnhp)' * ve) + (pep(indnhp,end) * ve );
			pip(indnhp,:) = (pip(indnhp,:)  - (pip(indnhp,end) * ve )) .* ...
					(fact(indnhp)' * ve) + (pip(indnhp,end) * ve );
		end

	end

	% securite
	tep = min(telim,max(t_sub * ve,real(pep ./ nep ./ ee)));
	tip = min(telim,max(t_sub * ve,real(pip ./ nip ./ ee)));

	%wapres =  (3/2) .* trapz(x,vpr .*(pep + pip),2);
	%figure(131);plot(temps,wth,temps,wapres,'o');keyboard
		
	
	% boucle de convergence interne
	deltaei   = Inf;
	f         = 0.3;
	tepmem    = tep;
	tipmem    = tip;
	xiioxie_mem = 0;
	Qeimem    = zeros(size(tep));
	klmax     = 101;
	%figure(71);clf
	
	for kl = 1:klmax

		% si xiioxie == 0 => calcul de la valeur
		if xiioxie <= 0                
		    if xiioxie < 0

				xiioxie = z0xiioxie(A,profli,tep,tip,abs(xiioxie));
                                tloc = 1:length(xiioxie);
                                %figure(113);hold on;
   			        %plot(tloc,xiioxie,'color',rand(1,3));
                                %drawnow 
                                xiioxie = xiioxie * ve;
 				% hypothese pour le transport
				xii     = xiioxie .* xie;
		    else
                        % attention : un seul appel par convergence (passe pas la convergence globale).
			% calcul de la balance TEM/ITG
			% calcul des longeurs gradients nromalisees
			lmin   = profli.rmx(:,end) ./ 21;
			lmax   = 18 .* pi .* profli.Raxe(:,end);
			ene    = min(1./lmin,max(1./lmax,2 .* abs(nep(:,5) - nep(:,end-5)) ./ max(nep(:,1),1e13) ./ profli.rmx(:,end))) .* profli.Raxe(:,end);
			eti    = min(1./lmin,max(1./lmax,2 .* abs(tip(:,5) - tip(:,end-5)) ./ max(tip(:,1),30) ./  profli.rmx(:,end))) .*  profli.Raxe(:,end);
			ete    = min(1./lmin,max(1./lmax,2 .* abs(tep(:,5) - tep(:,end-5)) ./ max(tep(:,1),30) ./  profli.rmx(:,end))) .*  profli.Raxe(:,end);
			Kt     = profli.ftrap(:,end) ./ (1 - profli.ftrap(:,end));
			% d'apres le travail de E. Asp
			% la partie en s et q n'est pas considerer ici
			eticr  = 2 ./ 3 .* ene +  20 ./ 9 .* tip(:,1)./ max(30,tep(:,1));
			etecr  = 20 ./ Kt + 2 ./ 3 .* ene  - Kt ./ 2  .* (1 - ene ./ 2 ) .^ 2;
			% il y a toujours de la turbulence sur les electron (niveau neoclassique ionique)
			dte    = max(eps,ete - etecr);
			% le niveau neoclssique sur les ions et le niveau de base sur les eletrons est le meme
			dti    = max(eps,eti - eticr);
			% la constante est ajustee sur TS
			xiioxie = min(10,max(0.1,(7/2) .* sqrt(tip(:,1)./ max(30,tep(:,1))) .* (1 + tanh(dti ./ ete))))*ve;
			% cette expression suppose que la turbulence electronique hors ITB n'est pas stabilisee
			% l'echelle choisiee est celle du gradient des electrons 
			xiioxie = 0.3 .* xiioxie + 0.7 .* xiioxie_mem;
			xiioxie_mem = xiioxie;
			
%  			tloc = 1:length(xiioxie);
%  			figure(113);hold on;
%  			plot(tloc,xiioxie,'color',rand(1,3));
%  			subplot(2,2,1);plot(tloc,xiioxie,'b');
%  			subplot(2,2,2);plot(tloc,dte,'r',tloc,dti,'g');
%  			subplot(2,2,3);plot(tloc,ete,'r',tloc,etecr,'b');
%  			subplot(2,2,4);plot(tloc,eti,'r',tloc,eticr,'b');
%			drawnow	
 			% hypothese pour le transport
			xii     = xiioxie .* xie;
                    end
                else
			% hypothese pour le transport
			xii     = xiioxie .* xie;
		end
                if isempty(xiioxie_ped)
		      xiioxie_ped = xiioxie;    
                end

		% uniquement pour le mode H car Xie et Xii en 0.95 conribue dans l'integrale sur grad(Te) pour calculer Te 
  		if (grad_ped >= 1) && isfield(profli,'xie') && ~isempty(indhp)
			rapxe                = mean(xie(:,1:end-2) ./ max(eps,profli.xie(:,1:end-2) .* profli.nep(:,1:end-2)),2);
			rapxi                = mean(xii(:,1:end-2) ./ max(eps,profli.xii(:,1:end-2) .* profli.nip(:,1:end-2)),2);
            %xie_ = xie;
            %xii_ = xii;
  			xie(indhp,end-1:end)= cat(2,rapxe(indhp),rapxe(indhp)) .* profli.xie(indhp,end-1:end) .* profli.nep(indhp,end-1:end);
  			xii(indhp,end-1:end)= cat(2,rapxi(indhp),rapxi(indhp)) .* profli.xii(indhp,end-1:end) .* profli.nip(indhp,end-1:end);
            %figure(41);plot(x,xie_,'r',x,xie,'b');figure(42);plot(x,xii_,'r',x,xii,'b');drawnow;keyboard
        end
		if (grad_ped == 2) || (grad_ped == 3)
                xii = xii .* nip ./ nep;
        end

        % terme source
        % equipartition
        switch  extended_qei
            case 'on'
                qeib0 = equipartition_full(tep,tip,nep,nHp,nDp,nTp,nhep,nzp,rimp .* nzp,(1 - Sn_fraction) .* nwp, Sn_fraction .* nwp,zimp,zmax,nbp,nhep3);
            otherwise
                warning off
                lnei          =  14.9 - 0.5.*log(nep ./ 1e20) + log(tep ./ 1e3);
                warning on
                ind = find(~isfinite(lnei) | (lnei <10));
                if ~isempty(ind)
                    lnei(ind) = 10 .* ones(1,length(ind));
                end
                taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (epsi0 .^ 2 ./ ee .^ 4 .* sqrt(me))  .* ...
                    ((ee .* max(t_sub * ve,tep)) .^ (3/2) ./ nep ./ lnei);
                if Sn_fraction > 0
                    factor_w_sn = (1 - Sn_fraction) .* z0wavez(tep) .^ 2 ./ 183.84 .* nwp + ...
                        Sn_fraction.* z0snavez(tep) .^ 2 ./ 118.71 .* nwp;
                    qeib0    = 3 .* me ./ mp ./ taues .* (nHp + nDp ./ 2 + nTp ./ 3 + nhep + nhep3 .* (4/3.02) +  ...
                        (zimp .^ 2 ./ Aimp) .* nzp + rimp .*  (zmax .^ 2 ./ Amax) .* nzp + factor_w_sn + 25/11 .* nbp);
                else
                    qeib0    = 3 .* me ./ mp ./ taues .* (nHp + nDp ./ 2 + nTp ./ 3 + nhep + nhep3 .* (4/3.02)  +  ...
                        (zimp .^ 2 ./ Aimp) .* nzp + rimp .*  (zmax .^ 2 ./ Amax) .* nzp + z0wavez(tep) ./ 183.84 .* nwp + 25/11 .* nbp);
                end
        end
        if isappdata(0,'TE_EXP') && ~isappdata(0,'TI_EXP') && isfield(profli,'qjli')
            tes = getappdata(0,'TE_EXP');
            tex = max(1,interp1_ex(tes.temps,tes.te,temps,'nearest','extrap'));
            tep_e  = pchip(tes.x,tex,x);
            qeib     = qeib0 .* ee .* (tep_e - tip);
        elseif     ~isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') && isfield(profli,'qjli')
            tis = getappdata(0,'TI_EXP');
            tix = max(1,interp1_ex(tis.temps,tis.ti,temps,'nearest','extrap'));
            tip_e = pchip(tis.x,tix,x);
            qeib     = qeib0 .* ee .* (tep - tip_e);
        else
            qeib     = qeib0 .* ee .* (tep - tip);
        end
% this is now taken into account in sources (controlled with parameter
% cx_ion in SOL section)
% 		% terme supplementaire pour l'equipartition du a la recombinaison
% 		if (grad_ped >= 2)
% 			trec    = [0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
% 			srec    = [7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 1e-6;% m^3/s;
% 			rrec    = reshape(exp(pchip(log(trec),log(srec),log(profli.tep(:)))),size(tep));
% 			nu_rec  = nep .* profli.n1p .* rrec;
% 			qeib    = qeib + nu_rec .* max(0,(3/2) .* tep - t_sub  * ve) .* 1.6022e-19; % passage en W/m^3
% 			%figure(51);clf;plot(x,qeib,'r',x,qeib - nu_rec .* max(0,tep - t_sub * ve) .* 1.6022e-19,'b',x,nu_rec .* max(0,tep - t_sub * ve) .* 1.6022e-19,'g');drawnow
% 			% manque les termes sur les autres impuretes
% 		end

		% les flux de puissance
        % si profil creux, le flux peut s'inverser
        if (creux == 1)
            Qe = cumtrapz(x,profli.source_el .* vpr,2);
            Qi = cumtrapz(x,profli.source_ion .* vpr,2);
        else
            Qe = max(0,cumtrapz(x,profli.source_el .* vpr,2));
            Qi = max(0,cumtrapz(x,profli.source_ion .* vpr,2));
        end
        Qei = cumtrapz(x,qeib .* vpr,2);
        if creux == 1
            Qei = min(1,0.99 .* abs(Qe + Qi) ./ (abs(Qei) + eps) ) .* Qei;
        end
        % stabilisation
        if kl > 1
            Qei = 0.3 .* Qei + 0.7 .* Qeimem;
        end
        Qeimem = Qei;
        % securite
        if (creux == 0)
            Qei = max(-0.9 .* Qi , min(0.9 .* Qe,Qei));
        end
        
                 %figure(51);clf;plot(x,mean(Qe,1),'r',x,mean(Qe-Qei,1),'m',x,mean(Qi,1),'b',x,mean(Qi+Qei,1),'c',x,mean(Qei,1),'g');drawnow
		
		% calcul de tibord/tebord (hypothese des flux)
		titebord    = max(0.5,min(2,nep(:,end) ./ nip(:,end) .* max(1,Qi(:,end) + Qei(:,end)) ./ max(1,Qe(:,end) - Qei(:,end)) ./ xiioxie_ped(:,end)));
		titebord    = max(0.5,min(2,0.1 .* titebord + 0.9 .* tip(:,end) ./ tep(:,end)));
		tibord      = max(min(13.3,t_sub),tebord .* titebord);
		
		% calcul de teped et tiped		
		%pepedr           = pped .* max(0.25,min(0.75,max(1,Qe(:,end-1) - Qei(:,end-1)) ./ ...
		%		   max(1,Qe(:,end-1) + Qi(:,end-1))));
		pipedr           = pped .* max(0.25,min(0.75,max(1,Qi(:,end-1) + Qei(:,end-1)) ./ ...
				   max(1,Qe(:,end-1) + Qi(:,end-1)) ./ xiioxie_ped(:,end)));
		pepedr           = pped - pipedr;				   
		teped            = max(tebord,min(0.9 .* tep(:,1),pepedr ./ nep(:,end-1) ./ ee));
		tiped            = max(tibord,min(0.9 .* tip(:,1),pipedr ./ nip(:,end-1) ./ ee));				
		%teped   	 = max(tebord,min(0.9 .* tep(:,1),0.1 .* teped + 0.9 .* tep(:,end-1)));
		%tiped   	 = max(tibord,min(0.9 .* tip(:,1),0.1 .* tiped + 0.9 .* tip(:,end-1)));
		
		if (grad_ped == 2) || (grad_ped == 3)

			%figure(71);hold on;plot(temps,teped .* nep(:,end-1) .* ee + tiped .* nip(:,end-1) .* ee,'r',temps,pped,'b',temps,ppedin,'g');drawnow
			% cas du mode L
			if ~isempty(indnhp)
				% temperature  electronique
				% si profil creux, le flux peut s'inverser
				if (creux == 1) 
				    inte     = - (Qe - Qei) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
				    inte_ext = - (Qe + Qi) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
				    inte     = max(inte,inte_ext);
				else
				    inte = -max(0,(Qe - Qei) ./ max(eps,xie)./ max(eps,vpr)) ./ grho2;
				end
				inte(:,1) = 0;
				tep0 = real(cumtrapz(x(:,end:-1:1),inte(:,end:-1:1),2));
				tep0 = tep0(:,end:-1:1);
				pep0 = ee .* nep .* tep0;
	
				% temperature  ionique
				if (creux == 1) 
				    inti     = -(Qi + Qei) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
				    inti_ext = -(Qe + Qi) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
				    inti     = max(inti,inti_ext);
				else
				    inti = -max(0,(Qi + Qei) ./ max(eps,xii)./ max(eps,vpr)) ./ grho2;
				end
				inti(:,1) = 0;
				tip0 = real(cumtrapz(x(:,end:-1:1),inti(:,end:-1:1),2));
				tip0 = tip0(:,end:-1:1);
				pip0 = ee .* nip .* tip0;
	
				% calcul de la constante de normalisation
				wloc_L    = (3/2) .* trapz(x,vpr .*(pep0 + pip0),2);
				w0loc_L   = ee .*(3/2) .* trapz(x,vpr .*(nep .* (tebord * ve) + nip .* (tibord * ve)),2);
	
				fact(indnhp)  = max(0,wth(indnhp) - w0loc_L(indnhp)) ./ max(eps,wloc_L(indnhp));
	
				tep(indnhp,:) = (fact(indnhp)' * ve) .* tep0(indnhp,:) + tebord(indnhp) * ve;
				pep(indnhp,:) = ee .* nep(indnhp,:) .* tep(indnhp,:);
	
				tip(indnhp,:) = (fact(indnhp)' * ve) .* tip0(indnhp,:) + tibord(indnhp) * ve;
				pip(indnhp,:) = ee .* nip(indnhp,:) .* tip(indnhp,:);
	
			end
		
			% cas du mode H
			if ~isempty(indhp)
				% temperature  electronique
				if (creux == 1) 
				  inte     = - (Qe - Qei) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
				  inte_ext = - (Qe + Qi) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
				  inte     = max(inte,inte_ext);
				else
				  inte = -max(0,(Qe - Qei) ./ max(eps,xie)./ max(eps,vpr)) ./ grho2;
				end
				inte(:,1) = 0;
				tep0 = real(cumtrapz(x((end-1):-1:1),inte(:,(end-1):-1:1),2));
				tep0 = tep0(:,end:-1:1);
				pep0 = ee .* nep(:,1:(end-1)) .* tep0;
	
				% temperature  ionique                
				if (creux == 1) && (kl > 3)
				  inti     = - (Qi + Qei) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
				  inti_ext = - (Qe + Qi) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
				  inti     = max(inti,inti_ext);
				else
				  inti = -max(0,(Qi + Qei) ./ max(eps,xii)./ max(eps,vpr)) ./ grho2;
				end
				inti(:,1) = 0;
				tip0 = real(cumtrapz(x(:,(end-1):-1:1),inti(:,(end-1):-1:1),2));
				tip0 = tip0(:,end:-1:1);
				pip0 = ee .* nip(:,1:(end-1)) .* tip0;
	
				% calcul de la constante de normalisation
				wloc_H    = (3/2) .* trapz(x(1:end-1),vpr(:,1:end-1) .*(pep0 + pip0),2);
				vpx       = cumtrapz(x,vpr,2);
				p0        = ee .* nep(:,end) .* tebord + ee .* nip(:,end) .* tibord;
	
				w0loc_H   = ee .*(3/2) .* trapz(x(1:(end-1)),vpr(:,1:(end-1)) .*(nep(:,1:(end-1)) .* (teped * ve(1:(end-1))) +  ...
								nip(:,1:(end-1)) .* (tiped * ve(1:(end-1)))),2) + ...
					(3/2) .* (pped + p0)  ./ 2 .* (vpx(:,end) - vpx(:,end-1));
	
				fact(indhp)  = max(0,wth(indhp) - w0loc_H(indhp)) ./ max(eps,wloc_H(indhp));
	
				tep(indhp,1:(end-1)) = (fact(indhp)' * ve(1:(end-1))) .* tep0(indhp,:) + teped(indhp) * ve(1:(end-1));
				tep(indhp,end) =  tebord(indhp);
				pep(indhp,:) = ee .* nep(indhp,:) .* tep(indhp,:);
	
				tip(indhp,1:(end-1)) = (fact(indhp)' * ve(1:(end-1))) .* tip0(indhp,:) + tiped(indhp) * ve(1:(end-1));
				tip(indhp,end) =  tibord(indhp);
				pip(indhp,:) = ee .* nip(indhp,:) .* tip(indhp,:);
	
			end

		else

			% temperature  electronique
			if (creux == 1) 
			    inte     = - (Qe - Qei) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
		            inte_ext = - (Qe + Qi) ./ max(eps,xie)./ max(eps,vpr) ./ grho2;
		            inte     = max(inte,inte_ext);
			else
			    inte = -max(0,(Qe - Qei) ./ max(eps,xie)./ max(eps,vpr)) ./ grho2;
			end 
			inte(:,1) = 0;
			tep = real(cumtrapz(x(:,end:-1:1),inte(:,end:-1:1),2));
			pep = nep .* tep(:,end:-1:1) ./ (max(nep,[],2) * ve);
			if ~isempty(indnhp)
				pep(indnhp,:) = pep(indnhp,:)  + (ee .* tebord(indnhp) .* nep(indnhp,end) - pep(indnhp,end)) * ve;
			end
			if ~isempty(indhp)
				pep(indhp,1:(end-1)) = pep(indhp,1:(end-1))  + (ee .* teped(indhp) .* nep(indhp,end - 1)- pep(indhp,end-1)) * ve(1:(end-1));
				pep(indhp,end) = ee .* tebord(indhp) .*  nep(indhp,end);
			end
			% temperature  ionique
			if (creux == 1) 
			    inti     = - (Qi + Qei) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
		            inti_ext = - (Qe + Qi) ./ max(eps,xii)./ max(eps,vpr) ./ grho2;
			    inti     = max(inti,inti_ext);
			else
			    inti = -max(0,(Qi + Qei) ./ max(eps,xii)./ max(eps,vpr)) ./ grho2;
			end
			inti(:,1) = 0;
			tip = real(cumtrapz(x(:,end:-1:1),inti(:,end:-1:1),2));
			pip = nip .* tip(:,end:-1:1) ./ (max(nip,[],2)*ve);
			if ~isempty(indnhp)
				pip(indnhp,:) = pip(indnhp,:)  + (ee .* tibord(indnhp) .* nip(indnhp,end) - pip(indnhp,end)) * ve;
			end
			if ~isempty(indhp)
				pip(indhp,1:(end-1)) = pip(indhp,1:(end-1))  + (ee .* tiped(indhp).* nip(indhp,end - 1) - pip(indhp,end-1)) * ve(1:(end-1));
				pip(indhp,end) = ee .* tibord(indhp).*  nip(indhp,end);
			end
	
	
			% energie non normalise
	
			wloc = (3/2) .* trapz(x,vpr .*(pep + pip),2);
			ppedloc = pep+pip;
			ppedloc(:,1:(end-1)) =(pep(:,end-1) + pip(:,end-1)) * ve(1:(end-1));
			wpiedloc  =  (3/2) .* trapz(x,vpr .* ppedloc,2);
			w0loc  =  (3/2) .* (pep(:,end) + pip(:,end)) .* trapz(x,vpr,2);
			if ~isempty(indhp)
				%fact(indhp)  = (wth(indhp) - wpied(indhp)) ./ max(1,wloc(indhp) - wpiedloc(indhp));
				fact(indhp)  = (wth(indhp) - wpiedloc(indhp)) ./ max(eps,wloc(indhp) - wpiedloc(indhp));
				pep(indhp,1:(end-1)) = (pep(indhp,1:(end-1)) - (pep(indhp,(end-1)) * ve(1:(end-1)) )) .*  ...
						(fact(indhp)' * ve(1:(end-1))) + (pep(indhp,end-1) * ve(1:(end-1)) );
				pip(indhp,1:(end-1)) = (pip(indhp,1:(end-1)) - (pip(indhp,(end-1)) * ve(1:(end-1)) )) .* ...
						(fact(indhp)' * ve(1:(end-1))) + (pip(indhp,end-1) * ve(1:(end-1)) );
			end
			if ~isempty(indnhp)
				%fact(indnhp)  = (wth(indnhp) - w0(indnhp)) ./ max(1,wloc(indnhp) - w0loc(indnhp));
				fact(indnhp)  = (wth(indnhp) - w0loc(indnhp)) ./ max(eps,wloc(indnhp) - w0loc(indnhp));
				pep(indnhp,:) = (pep(indnhp,:)  - (pep(indnhp,end) * ve )) .*  ...
						(fact(indnhp)' * ve) + (pep(indnhp,end) * ve );
				pip(indnhp,:) = (pip(indnhp,:)  - (pip(indnhp,end) * ve )) .* ...
						(fact(indnhp)' * ve) + (pip(indnhp,end) * ve );
			end

		end


		% temperature  electronique
		% securite
		tep = min(telim,max(max(tebord * ve,t_sub * ve),real(pep ./ nep ./ ee)));
		tip = min(telim,max(max(tibord * ve,t_sub * ve),real(pip ./ nip ./ ee)));
        
%         % just to testing
%         % calcul des coefficients de transport associes
%         tepd1_       = pdederive(x,tep,0,2,2,1);
%         tipd1_       = pdederive(x,tip,0,2,2,1);
%         % attention vpr = V' * rhomax !
%         xiea_         = - (Qe - Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tepd1_)) ./ ...
%             max(1e13,nep) ./ grho2 .* (rhomax * ve) .^ 2;
%         xiia_        = - (Qi  + Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tipd1_)) ./  ...
%             max(1e13,nip) ./ grho2.* (rhomax * ve) .^ 2;
%         ki_ = round(size(xie,1) /2);
%         figure(21);plot(x,xie(ki_,:),'r',x,xii(ki_,:),'b-.',x,xiea_(ki_,:),'m',x,xiia_(ki_,:),'c-.');drawnow
		
		% effet des dents de scie resolu en temps
		% maintenant l'information est transmise par xieshape
		% cas des DDS sur Pe et Pion
		if any(indice_inv > 1.5) && (qdds < 0) && (size(nep,1) == 4) 
			% pression finale
			pep = ee .* tep .* nep;
			pip = ee .* tip .* nip;
			% nouveau profil de pression electronique
			dpepdx                 = cat(2,0 * vt,diff(pep - pep(:,end) * ve,1,2));
			mask                   = (vt * (1:length(x))) > (indice_inv * ve); 
			dpepdx                  = dpepdx .* mask;
			pepnew                  = cumsum(dpepdx,2);
			pepnew                  = pepnew - pepnew(:,end) * ve + pep(:,end) *ve ;
			% normalisation (conservation de l'energie)
			dwe                    = trapz(x,(pepnew - pep) .* profli.vpr,2);
			maskc                   = (vt * (1:length(x))) <= (indice_inv * ve); 
			maskc(:,1:2)           = 1;
			dvp                    = trapz(x,maskc .* profli.vpr,2);
			dpep                   = dwe ./ max(eps,dvp);
			pep_mem                = pep;
			pep                    = min(pep_mem,pepnew - (dpep * ve) .*  (~mask));
			%  
			% nouveau profil de pression ionique
			dpipdx                 = cat(2,0 * vt,diff(pip - pip(:,end) * ve,1,2));
			dpipdx                  = dpipdx .* mask;
			pipnew                  = cumsum(dpipdx,2);
			pipnew                  = pipnew - pipnew(:,end) * ve + pip(:,end) *ve ;
			% normalisation (conservation de l'energie)
			dwi                    = trapz(x,(pipnew - pip) .* profli.vpr,2);
			dvp                    = trapz(x,maskc .* profli.vpr,2);
			dpip                   = dwi ./ max(eps,dvp);
			%nep_mem                = nep;
			pip_mem                = pip;
			pip                    = min(pip_mem,pipnew - (dpip * ve) .*  (~mask));
			%		
			%figure(21)
			%plot(1:size(pep,1),trapz(x,pep .* profli.vpr,2)-trapz(x,pep_mem .* profli.vpr,2),'r', ...
			%	1:size(pip,1),trapz(x,pip .* profli.vpr,2)-trapz(x,pip_mem .* profli.vpr,2),'b');
			% temperature finale
            tep_mem = tep;
			tep = min(tep_mem,max(max(tebord * ve,t_sub * ve),real(pep ./ nep ./ ee)));
            tip_mem = tip;
			tip = min(tip_mem,max(max(tibord * ve,t_sub * ve),real(pip ./ nip ./ ee)));
            
%             figure(21)
% 			clf
%             subplot(2,2,1)
% 			plot(x,pep,'r',x,pep_mem,'b')
%             subplot(2,2,2)
% 			plot(x,pip,'r',x,pip_mem,'b')
%             subplot(2,2,3)
% 			plot(x,tep,'r',x,tep_mem,'b')
%             subplot(2,2,4)
% 			plot(x,tip,'r',x,tip_mem,'b')
% 			drawnow

		end

	        % petite perturbation entre 2 appel pour la convergence
                %		f   = 1 ./ (kl + 2);
                
% 			wel  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
% 			wion    = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);
%             figure(24);plot(temps,wth,'r',temps,wel,'b',temps,wion,'k',temps,wel+wion,'g');drawnow

			
		% amorti
	    tep =  min(telim,max(max(tebord * ve,t_sub * ve),real(f .* tep + tepmem .* (1-f))));		
		tip =  min(telim,max(max(tibord * ve,t_sub * ve),real(f .* tip + tipmem .* (1-f))));

		deltae  =  max(abs(tep(:) - tepmem(:))) ./ mean(tepmem(:));
		deltai  =  max(abs(tip(:) - tipmem(:))) ./ mean(tipmem(:));
		deltaeimem = deltaei;
		deltaei =  max(deltae,deltai);

		tepmem = tep;
		tipmem = tip;


		% pression finale
		pep = ee .* tep .* nep;
		pip = ee .* tip .* nip;
		
		% diminution de f
		f = 0.7 .* f;

		% fin boucle de convergence
		%disp([kl,f,deltaei]);
		if (deltaei < 1e-3) & (kl > 11)
			break
		end
		
	end
	if kl >= klmax
		fprintf('e');
    end
    
    % memorisation
    tep_d = tep;
    tip_d = tip;
    
    if isappdata(0,'TE_EXP') & ~isappdata(0,'TI_EXP') & isfield(profli,'qjli');
        tes = getappdata(0,'TE_EXP');
        tex = max(1,interp1_ex(tes.temps,tes.te,temps,'nearest','extrap'));
        tep  = pchip(tes.x,tex,x);
        
        % recompute Ti to match wth
        wel  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
        wion_th = max(eps,wth - wel);
        wion    = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);
        %figure(23);plot(temps,wth,'r',temps,wel,'b',temps,wion,'k',temps,wel+wion,'m',temps,wion_th,'c');drawnow
                
        % cas du mode L
        flim = 0.001;
        if ~isempty(indnhp)
            wion0   = (3/2) .* ee .* trapz(x,vpr .* nip .* (tip(:,end) * ve),2);
            % compute coherent tibord with available energy
            tibord = max(t_sub,tip(:,end) .* min(1, (0.5 .* wion_th ./ wion0))); 
            wion0   = (3/2) .* ee .* trapz(x,vpr .* nip .* (tibord * ve),2);
            fact_e = max(eps,wion_th(indnhp) - wion0(indnhp)) ./ max(eps,wion(indnhp) - wion0(indnhp));
            fact_e = max(flim,min(1/flim,fact_e));
            %figure(33);plot(temps(indnhp),fact_e);drawnow
            tip(indnhp,1:end-1) =  (fact_e * ve(1:end-1))  .* (tip(indnhp,1:end-1) - tibord(indnhp) * ve(1:end-1)) + ...
                tibord(indnhp) * ve(1:end-1);
            tip(indnhp,end) = tibord(indnhp);
        end
        % cas du mode H
        if ~isempty(indhp)
            tip_star = tip(:,end-1) * ve;
            tip_star(:,end) = tip(:,end);
            wion0  = (3/2) .* ee .* trapz(x,vpr .* nip .* tip_star,2);
            % compute coherent tibord and tiped with available energy
            tip_star(:,end-1:end) = max(t_sub*ones(1,2),tip_star(:,end-1:end) .* (min(1, (0.8 .* wion_th ./ wion0))*ones(1,2))); 
            wion0  = (3/2) .* ee .* trapz(x,vpr .* nip .* tip_star,2);
            fact_e = max(eps,wion_th(indhp) - wion0(indhp)) ./ max(eps,wion(indhp) - wion0(indhp));
            %figure(33);plot(temps(indhp),fact_e);drawnow
            fact_e = max(flim,min(1/flim,fact_e));
            tip(indhp,1:end-2) =  (fact_e * ve(1:end-2))  .* (tip(indhp,1:end-2) - tip_star(indhp,end-1) * ve(1:end-2)) + ...
                tip_star(indhp,end-1) * ve(1:end-2);
            tip(indhp,end-1:end) = tip_star(indhp,end-1:end);
        end
        
        wion   = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);
        %figure(22);plot(temps,wth,'r',temps,wel,'b',temps,wion,'k',temps,wel+wion,'g');drawnow
        
    end
    
    if ~isappdata(0,'TE_EXP') & isappdata(0,'TI_EXP') & isfield(profli,'qjli');
        tis = getappdata(0,'TI_EXP');
        tix = max(1,interp1_ex(tis.temps,tis.ti,temps,'nearest','extrap'));
        tip = pchip(tis.x,tix,x);
        
        % recompute Te to match wth
        wion = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);

        wel_th = max(eps,wth - wion);
        wel    = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
        %figure(23);plot(temps,wth,'r',temps,wel,'b',temps,wion,'k',temps,wel+wion,'m',temps,wel_th,'c');drawnow
        
        % cas du mode L
        flim = 0.001;
        if ~isempty(indnhp)
            wel0   = (3/2) .* ee .* trapz(x,vpr .* nep .* (tep(:,end) * ve),2);
            fact_e = max(eps,wel_th(indnhp) - wel0(indnhp)) ./ max(eps,wel(indnhp) - wel0(indnhp));
            fact_e = max(flim,min(1/flim,fact_e));
            %figure(33);plot(temps(indnhp),fact_e);drawnow
            tep(indnhp,1:end-1) =  (fact_e * ve(1:end-1))  .* (tep(indnhp,1:end-1) - tep(indnhp,end) * ve(1:end-1)) + ...
                tep(indnhp,end) * ve(1:end-1);
        end
        % cas du mode H
        if ~isempty(indhp)
            tep_star = tep(:,end-1) * ve;
            tep_star(:,end) = tep(:,end);
            wel0  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep_star,2);
            fact_e = max(eps,wel_th(indhp) - wel0(indhp)) ./ max(eps,wel(indhp) - wel0(indhp));
            %figure(34);plot(temps(indhp),fact_e);drawnow
            fact_e = max(flim,min(1/flim,fact_e));
            tep(indhp,1:end-2) =  (fact_e * ve(1:end-2))  .* (tep(indhp,1:end-2) - tep(indhp,end-1) * ve(1:end-2)) + ...
                tep(indhp,end-1) * ve(1:end-2);
        end
        
        wel   = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
        %figure(23);plot(temps,wth,'r',temps,wel,'b',temps,wion,'k',temps,wel+wion,'g');drawnow
        
        
    end

    % update Qei with external data
    if (isappdata(0,'TE_EXP') & ~isappdata(0,'TI_EXP') & isfield(profli,'qjli')) || ...
            (~isappdata(0,'TE_EXP') & isappdata(0,'TI_EXP') & isfield(profli,'qjli'))
        qeib     = qeib0 .* ee .* (tep - tip);
        
% this is now taken into account in sources (controlled with parameter
% cx_ion in SOL section)
%         % terme supplementaire pour l'equipartition du a la recombinaison
%         if (grad_ped >= 2)
%             trec    = [0.1   ,1       ,8     ,20     ,50     ,100     ,300   ,500   ,950   ,1e5]; % eV
%             srec    = [7e-13 ,1.7e-13 ,4e-14 ,2e-14  ,8e-15  ,4e-15   ,1e-15 ,5e-16 ,2e-16 ,1e-19] .* 1e-6;% m^3/s;
%             rrec    = reshape(exp(pchip(log(trec),log(srec),log(profli.tep(:)))),size(tep));
%             nu_rec  = nep .* profli.n1p .* rrec;
%             qeib    = qeib + nu_rec .* max(0,(3/2) .* tep - t_sub  * ve) .* 1.6022e-19; % passage en W/m^3
%             %figure(51);clf;plot(x,qeib,'r',x,qeib - nu_rec .* max(0,tep - t_sub * ve) .* 1.6022e-19,'b',x,nu_rec .* max(0,tep - t_sub * ve) .* 1.6022e-19,'g');drawnow
%             % manque les termes sur les autres impuretes
%         end
        
        % les flux de puissance
        % si profil creux, le flux peut s'inverser
        Qei = cumtrapz(x,qeib .* vpr,2);
        if creux == 1
            Qei = min(1,0.99 .* abs(Qe + Qi) ./ (abs(Qei) + eps) ) .* Qei;
        else
            Qei = max(-0.9 .* Qi , min(0.9 .* Qe,Qei));
        end
    end
    
%  disp([kl,f,deltaei]);
%  	figure(51);plot(temps,tebord,'b',temps,tep(:,end),'ob', ...
%  	                temps,tibord,'r',temps,tip(:,end),'or', ...
%  	                temps,teped,'c',temps,tep(:,end-1),'oc', ...
%  	                temps,tiped,'m',temps,tip(:,end-1),'om', ...
%  	                temps,tep(:,1),'r',temps,tip(:,1),'b');hold;drawnow;	
%  	figure(51);plot(temps,tebord,'b',temps,tep(:,end),'ob', ...
%  	                temps,tibord,'r',temps,tip(:,end),'or');set(gca,'xlim',[0 3]);drawnow
%  	                
%  			
%  	figure(71);clf; plot(x,tep,'r',x,tip,'b',x,tep_d,'m',x,tip_d,'c'); drawnow
	
	% calcul des coefficients de transport associes
	tepd1       = pdederive(x,tep,0,2,2,1);
	tipd1       = pdederive(x,tip,0,2,2,1);
	% attention vpr = V' * rhomax !
	xiea         = - (Qe - Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tepd1)) ./ ...
	                max(1e13,nep) ./ grho2 .* (rhomax * ve) .^ 2;
	xiia         = - (Qi  + Qei) ./ max(eps,vpr) ./ (ee .* min(-30,tipd1)) ./  ...
	                max(1e13,nip) ./ grho2.* (rhomax * ve) .^ 2;
	%  prolongation
	xiea         = pchip(x(2:end),real(xiea(:,2:end)),x);
	xiia         = pchip(x(2:end),real(xiia(:,2:end)),x);

%  figure(21);clf
%  plot(x,mean(xie,1),'r',x,mean(xii,1),'b',x,mean(xiea,1),'m',x,mean(xiia,1),'c')
%  drawnow
	
	% limitation
	profli.xie  = max(1e-3,min(30,xiea));
	profli.xii  = max(1e-3,min(30,xiia));
	
	% flux de chaleur en sortie :
	%profli.qe   = 	(Qe - Qei);
	%profli.qi   = 	(Qi + Qei);
        % vrai sans les securites
        profli.qe   = 	cumtrapz(x,profli.source_el .* vpr,2) - Qei;
        profli.qi   = 	cumtrapz(x,profli.source_ion .* vpr,2) + Qei;
	profli.qei  = 	Qei;

	% calcul de tite
	tite = max(1./rr,min(rr,trapz(x,tip .* vpr,2) ./ trapz(x,tep .* vpr,2)));

	% pei
	pei  = Qei(:,end);

	% temps de confinement
	Qmin = max(1,1e-3 .* (abs(Qe(:,end)) + abs(Qei(:,end)) + abs(Qi(:,end))) ./ 3);
	tauee = ee .* trapz(x,vpr .* nep .* tep,2) ./ max(Qmin,Qe(:,end) - Qei(:,end));
	tauii = ee .* trapz(x,vpr .* nip .* tip,2) ./ max(Qmin,Qi(:,end) + Qei(:,end));
	tauei = ee .* trapz(x,vpr .* nip .* tip,2) ./ max(Qmin,abs(pei));

	wel  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
	wion = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);
	 
	devw = abs((wel+wion-wth)./max(1,max(wth)))> 0.01;


	if any(devw)
%  		indok = find(~devw);
%  		indnok = find(devw);
		if creux ==1
		  indok = find(~devw);
		  indnok = find(devw);
		  if (length(indok) > 3)
		      tepf = sgolayfilt(tep(indok,:),1,3,[],1);
		      tipf = sgolayfilt(tip(indok,:),1,3,[],1);
		      %figure(27);plot(temps(indok),tepf,'r',temps,tep,'b');
		      tep(indnok,:) = interp1(temps(indok),tepf,temps(indnok),'linear','extrap');
		      tip(indnok,:) = interp1(temps(indok),tipf,temps(indnok),'linear','extrap');
		      
		      % si non corrige
		      wel  = (3/2) .* ee .* trapz(x,vpr .* nep .* tep,2);
		      wion = (3/2) .* ee .* trapz(x,vpr .* nip .* tip,2);
	 
	              if any(abs((wel+wion-wth)./max(1,max(wth)))> 0.01)
			  fprintf('E');		  	              
	              end		      
		  else
		    fprintf('E');		  
		  end
		else
		    fprintf('E');
		end
%		if (length(indok) > 3)
%  			tepnew   = pchip(temps(indok)',tep(indok,:)',temps')';
%  			indnn    = find(all(tepnew > 13.6,2));
%  			tep      = profli.tep;
%  			tep(indnn,:) = tepnew(indnn,:);
%  			
%  			tipnew   = pchip(temps(indok)',tip(indok,:)',temps')';
%  			indnn    = find(all(tipnew > 13.6,2));
%  			tip      = profli.tip;
%  			tip(indnn,:) = tipnew(indnn,:);	
			
%		else
%			fprintf('E');
%  			indmax = find(max(devw)==devw,1);
%  			fprintf('Warning : energy mismatch in zeqie0prof (max = %g J @ %g s) \n', ...
%  		       	abs(wel(indmax)+wion(indmax)-wth(indmax)),temps(indmax));
%  			figure(61);clf;
%  			plot(temps,wel,'r',temps,wion,'b',temps,wel+wion,'g',temps,wth,'ko', ...
%  		     		temps(indnok),wth(indnok),'k*');
%  			indmax = find(max(devw)==devw);
%  			title(sprintf('Warning : energy mismatch in zeqie0prof (max = %g J @ %g s)', ...
%  		       	abs(wel(indmax)+wion(indmax)-wth(indmax)),temps(indmax)));
%  			drawnow
%    			keyboard
%			
%  		end	
	end
	
	% recopie
	profli.tep = max(1,real(tep));
	profli.tip = max(1,real(tip));
	

else

    % just for initialisation of the convergence loop
    % calcul du facteur z^2/a
    if frhe0 > 0
        z2sa = max(0,n1m - nDm - nTm) + (1/2) .* nDm + (1/3) .* nTm + nhem .* (4/3.02) + nem .* frhe0  + (zimp .^ 2 ./ Aimp) .* nimpm + rimp .* (zmax .^ 2 ./ Aimp) .* nimpm + (25/11) .* nboronm;
        z2sa = z2sa ./ (n1m + (1 + rimp) .* nimpm + nhem + nem .* frhe0 + nboronm);        
    else
        z2sa = max(0,n1m - nDm - nTm) + (1/2) .* nDm + (1/3) .* nTm + nhem  + (zimp .^ 2 ./ Aimp) .* nimpm + rimp .* (zmax .^ 2 ./ Aimp) .* nimpm;
        z2sa = z2sa ./ (n1m + (1 + rimp) .* nimpm + nhem + nboronm);
    end
    indnok = find(~isfinite(z2sa));
    repli = ( 1 + (A == 4)) .^ 2 ./ A;
    z2sa(indnok) = repli(indnok);
    
	if isfield(profli,'tep')
		x   = profli.xli;
		vpr = profli.vpr;
		ux  = (1  - x .^ 2);
		ve  = ones(size(x));
		vt  = ones(size(te0));
		nep = profli.nep;
		nip = profli.nip;
		tep = profli.tep;
		tip = profli.tip;

	else
		x   = linspace(0,1,21);
		vpr   = (2 .* vp) * x;
		ux  = (1  - x .^ 2);
		ve  = ones(size(x));
		vt  = ones(size(te0));
		nep = max(1e13,((nem .* (1 + ane)-nebord) * ve)  .* (vt * ux)  .^ (ane * ve) + nebord * ve);
                nep(:,end) = nebord;
		tep = ((te0 - tebord) * ve)  .* (vt * ux)  .^ (ate *ve) + tebord * ve;
		nip = nep .* (ae * ve);
		tip = tep .* (tite_in * ve);
    end

    
    

	% recalcul local de wth pour coherence
	wth = 3./2 .* ee.* trapz(x,(nep .*tep +nip .* tip) .* vpr,2);

	% securite sur pei
	pei  = max(-pion .* (1-1/rr),min(pel .* (1-1/rr),pei_in));


	% calcul de taue tauee et tauii coherent
	pth  = max(1,pel + pion);
	taue = wth ./ pth;
	% hypothese pour le transport
	if xiioxie <= 0
		xiioxie = 1;
	end
	aa      = 1./xiioxie; % du modele bgbs (seule hypothese reglable)
	%aa     = sqrt((1 + ane) ./ (1 + ate)); % dependance en piquage
	%figure(51);clf;plot(aa);drawnow;
	% par definition we = tauee * Pel et wion = tauii * Pion : wth = we +wion
	tauee   = wth ./ max(1, aa .* (pion + pei) + (pel-pei));
	tauii   = aa .* tauee;
	wion    = tauii .* (pion+pei);
	we      = tauee .* (pel-pei);

    
    % juste to initialise profile on first loop
    % no need of relativistic correction
	warning off
	lnei          =  15.2 - 0.5 .* log(nep ./ 1e20) + log(tep ./ 1e3);
	warning on
	ind = find(~isfinite(lnei) | (lnei <10));
	if ~isempty(ind)
		lnei(ind) = 10 .* ones(1,length(ind));
	end
	taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (epsi0 .^ 2 ./ ee .^ 4 .* sqrt(me))  .* ...
			((ee .* tep) .^ (3/2) ./ nep ./ lnei);
	qei0    = 3 .* me ./ mp ./ taues .* nip .* (z2sa * ve) .* (ee .* tep);
	pei0    = trapz(x,qei0 .* vpr,2);


	% calcul de wion
	kinv  = 3/2 .* ee  .* trapz(x,tep .* nip .* vpr,2);
	[wion,dwiondt] = zdwdt0(wion,pion+pei,tauii,temps,vp);
	% temps d'equipartition
	% puissance minimale sur un canal
	pmin   = max(1,abs(pei))./rr;
	tauei  = wion ./ max(pmin,abs(pei - dwiondt));
	% tite filtre
	tite   = tauii .* (pion + pei0) ./ kinv ./ ( 1+ tauii .* pei0 ./ kinv);

	% securite tite
	tite = min(rr,max(1/rr,tite));


	% calcul final de pei
	pei  = pei0 .* (1 -tite);
	% securite sur pei
	pei  = max(-pion .* (1-1/rr),min(pel .* (1-1/rr),pei));
	
	
end

function te_out = compute_temp(xli,ted1_out,tep)

ve     = ones(size(xli));

te_out = fliplr(cumtrapz(fliplr(xli),fliplr(ted1_out),2));
te_out = te_out  + (tep(:,end - 1) - te_out(:,end-1)) * ve;
te_out(:,end) = tep(:,end);
te_out(:,1) = te_out(:,2);


%ref: Wesson
function chii_neo = neo_Hinton(profil,Amain)

ve     = ones(size(profil.xli));
vt     = ones(size(profil.ptot,1),1);
epsi   = profil.epsi;
delta  = profil.Raxe - profil.Raxe(:,end) * ve;
a      = (profil.Raxe(:,end) .* epsi(:,end)) * ve;
delta_prim = pdederive(profil.xli,delta,0,2,2,1) ./ a;

f1     = (1 + 3/2 .* (epsi .^ 2 +  epsi .* delta_prim) + 3/8 .* epsi .^ 3 .* delta_prim) ./ ( 1 + 1/2 .* epsi .*  delta_prim);
f2     = sqrt(1 - epsi .^ 2) .* (1 + epsi .* delta_prim ./ 2) ./ (1 + delta_prim ./ epsi .* (sqrt(1 - epsi .^ 2)  - 1));

nI_ZI2    = (profil.nep .* profil.zeff - profil.n1p);

alpha  = nI_ZI2 ./ profil.nip;

lnldi   = 17.3 - 0.5 .* log(profil.nep./1e20) + log(profil.tip ./1e3);
tau_i   = 6.60e17 .* sqrt(Amain * ve) .* (profil.tip ./1e3).^(3/2) ./ profil.nip ./ lnldi;
vth_i   = 3.09e5  .* sqrt(profil.tip ./1e3 ./ (Amain*ve));
B       = sqrt(profil.bpol .^ 2 + (profil.fdia .* profil.ri) .^ 2);
rho_i   = 4.57e-3 .* sqrt((Amain * ve).* profil.tip ./1e3) ./ B;
nustar  = profil.Raxe .* profil.qjli ./ (epsi .^ (3/2) .* vth_i .* tau_i);
mustar_i  = nustar .* (1 + 1.54 .* alpha);

part1   = 0.66 .* (1 + 1.54 .* alpha) + (1.88 .* sqrt(epsi) - 1.54 .* epsi) .* (1 + 3.75 .* alpha);
part1   = part1 ./ (1 + 1.03 .* sqrt(mustar_i) + 0.31 .* mustar_i);

part2   = 0.59 .* mustar_i .* epsi ./ (1 + 0.74 .* mustar_i .* epsi .^ (3/2));
part2   = part2 .* (1 + (1.33 .* alpha .* (1 + 0.60 .* alpha)) ./ (1 + 1.79 .* alpha));

chii_neo = epsi .^(-3/2) .* profil.qjli .^ 2 .* rho_i .^ 2 ./ tau_i .* ...
    (part1 .* f1 + part2 .* (f1 -f2));


% chii_neo(:,1)  = chii_neo(:,2);
% figure(124); 
% plot(profil.xli,chii_neo);
% drawnow
% 
