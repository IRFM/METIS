% ZERODEVOLUTION : version avec evolution temporelle de METIS
%-------------------------------------------------------------
% fonction Matlab 7; zerodevolution.m -> zerodevolution
%
% Cette fonction implemente la version avec evolution temporelle de METIS et elle
% peux etre implemente dans une boucle avec un asservissement ou sous simulink
%
% syntaxe : 
%
%   initialisation :       
%       [zs,profil,z0dstruct] = zerodevolution([],option,temps_ini,cons,geo,[exp0d,sepa,xdur]);
%
%   calcul :
%       [zs,profil,z0dstruct] = zerodevolution(z0dstruct,option,temps,cons,geo,[exp0d,sepa,xdur]);
%
% entrees :
%
%     z0dstruct               = structure de calcul au temps precedents 
%     option                  = structure d'option de METIS (z0dinput.option)
%     temps                   = temps du caclul 
%     cons                    = structure de consignes de METIS (z0dinput.cons)
%     geo                     = structure pour la geometrie de METIS (z0dinput.geo)
%     exp0d                   = structure de donnees experimentales de METIS pour la comparaison (z0dinput.exp0d, optionnelle)
%     sepa                    = structure donnant la separatrice
%     xdur                    = structure donnant le profil xdur
%        
%
% sorties :
%      zs                     = structure de donnees scalaire de METIS pour le  temps "temps".             
%      profil                 = structure de profils de METIS pour le  temps "temps". 
%      z0dstruct              = structure de calcul completee pour le temps "temps"
%
% remarque : si cons et geo sont des structure a 1 temps, la valeur de cons.temps doit etre identique a temps.
%
% test
%      [option,cons1t,geo1t] = zerodevolution;
%      [zs,profil,z0dstruct] = zerodevolution([],option,1,cons1t,geo1t);
%      [zs,profil,z0dstruct] = zerodevolution(z0dstruct,option,2,cons1t,geo1t);
%
% fonction ecrite par J-F Artaud
% version CVS (creation le 28/08/2007)
%-----------------------------------------------------------------------
%
function [zs,profil,z0dstruct] = zerodevolution(z0dstruct,option,time,cons1t,geo1t,exp0d_in,sepa1t,xdur1t)

% test des entrees
if nargin == 0
	z0dinput 	= zerod_init(-2,0,2,0);
	zs              = z0dinput.option;
	profil          = z0dinput.cons;
	z0dstruct       = z0dinput.geo;
	noms = fieldnames(z0dstruct);
	for k=1:length(noms)
		if isempty(z0dstruct.(noms{k}))
			z0dstruct.(noms{k}) = NaN;
		end
	end
	return
end
if nargin < 6
	exp0d_in = [];
end
if nargin < 7 
	sepa1t = [];
end
if nargin < 8 
	xdur1t = [];
end

if ~isfield(option,'ploss_exp')
	option.ploss_exp ='with_prad';
end
%  if ~isfield(option,'inloop')
%  	option.inloop = 0;
%  end
%  inloop = option.inloop;

% nombre de tempes interne 
nbte = 3;

% selon la structure d'entree
% initialisation
if isempty(z0dstruct)
	% initilisation
        t0              = time - 0.01 .* geo1t.a(1) .* geo1t.R(1);
	temps    	= linspace(t0,time,nbte)';
	temps           = cat(1,temps,2 .* temps(end) - temps(end-1));
	vt              = ones(size(temps));
	z0dinput 	= zerod_init(-2,0,option.gaz,temps);
	if ~isempty(option)
		noms = fieldnames(option);
		for k=1:length(noms)
			if isfield(z0dinput.option,noms{k})
				z0dinput.option.(noms{k}) = option.(noms{k});
			end
		end
	end
	%z0dinput.option = option;
	noms = fieldnames(cons1t);
	cons1t.temps    = temps;
	for k=1:length(noms)
		if isfield(z0dinput.cons,noms{k})
			z0dinput.cons.(noms{k})(:) = vt * cons1t.(noms{k})(1);
		end
	end
	z0dinput.cons.temps = temps;
	%
	if ~isempty(geo1t)
		noms = fieldnames(geo1t);
		for k=1:length(noms)
			if isfield(z0dinput.geo,noms{k})
                		if ~isempty(geo1t.(noms{k})) 
					if  ~isempty(z0dinput.geo.(noms{k}))
                    				z0dinput.geo.(noms{k})(:) = vt *  geo1t.(noms{k})(1);
					else
                    				z0dinput.geo.(noms{k}) = vt *  geo1t.(noms{k})(1);					
					end
                		end
			end
		end
	end
	
	% mise a jour de la structure experimentale vide
	noms = fieldnames(z0dinput.zsinfo);
	exp0d  = [];
	texp   = z0dinput.cons.temps;
	exp0d.temps = texp;
	vtnan  = texp .* NaN;
	for k = 1:length(noms)
		nomc = noms{k};
		if isfield(exp0d_in,nomc)
			var = getfield(exp0d_in,nomc);
			if ~isempty(var) && all(isfinite(exp0d_in.temps))
				warning off
				var = interp1(exp0d_in.temps,var,texp,'nearest');
				warning on
			else
				var = vtnan;
			end
			exp0d  = setfield(exp0d,nomc,var);
		else
			exp0d = setfield(exp0d,nomc,vtnan);
		end
	end

	% donnees experimentale
	z0dinput.exp0d = exp0d;

	
	% utilisation de la sepratrice si disponible
	if isfield(sepa1t,'Rsepa') & isfield(sepa1t,'Zsepa')
		z0dinput.exp0d.Rsepa = vt * sepa1t.Rsepa(1,:);
		z0dinput.exp0d.Zsepa = vt * sepa1t.Zsepa(1,:);
		fprintf('$')
	    fprintf('$')
	elseif isfield(exp0d_in,'Rsepa') & isfield(exp0d_in,'Zsepa')
		z0dinput.exp0d.Rsepa = vt * exp0d_in.Rsepa(1,:);
		z0dinput.exp0d.Zsepa = vt * exp0d_in.Zsepa(1,:);
		fprintf('$')
	    fprintf('$')
	end

	% activation XDUR si disponible
	if isfield(xdur1t,'XDURx') & isfield(xdur1t,'XDURv')
		z0dinput.exp0d.XDURt = temps;
		z0dinput.exp0d.XDURx = vt * xdur1t.XDURx(1,:);
		z0dinput.exp0d.XDURv = vt * xdur1t.XDURv(1,:);	
	elseif isfield(exp0d_in,'XDURx') & isfield(exp0d_in,'XDURv')
		z0dinput.exp0d.XDURt = temps;
		indx   = find(temps(end)<= exp0d_in.XDURt,1);
		if isempty(indx)
			indx = 1;
		end
		if size(exp0d_in.XDURx,1) > 1
		  z0dinput.exp0d.XDURx = exp0d_in.XDURx(indx,:);
		else
		  z0dinput.exp0d.XDURx = exp0d_in.XDURx;
		end
		z0dinput.exp0d.XDURv = vt * exp0d_in.XDURv(indx,:);	
			
	end
		
	% calcul initial 
	z0dinput.option.init_geo  = 1;
	z0dinput.option.evolution = 0;	
	z0dinput.option.vloop     = 0;	
	fprintf('E INIT ->')
	[z0dstruct.zs,z0dstruct.info,z0dstruct.profil] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d,1e-3);
	% seul les end-1 temps sont valides 
	noms = fieldnames(z0dstruct.zs);
	for k = 1:length(noms)
		if length(z0dstruct.zs.(noms{k})) > 1
			z0dstruct.zs.(noms{k}) = z0dstruct.zs.(noms{k})(1:(end-1));
		end
	end
	noms = fieldnames(z0dstruct.profil);
	for k = 1:length(noms)
		if size(z0dstruct.profil.(noms{k}),1) > 1
			z0dstruct.profil.(noms{k}) = z0dstruct.profil.(noms{k})(1:(end-1),:);
		end
	end
	noms = fieldnames(z0dinput.cons);
	for k = 1:length(noms)
		if length(z0dinput.cons.(noms{k})) > 1
			  z0dinput.cons.(noms{k}) = z0dinput.cons.(noms{k})(1:(end-1));
		end
	end
	noms = fieldnames(z0dinput.geo);
	for k = 1:length(noms)
		if length(z0dinput.geo.(noms{k})) > 1
			  z0dinput.geo.(noms{k}) = z0dinput.geo.(noms{k})(1:(end-1));
		end
	end
	noms = fieldnames(z0dinput.exp0d);
	for k = 1:length(noms)
		if size(z0dinput.exp0d.(noms{k}),1) > 1
			z0dinput.exp0d.(noms{k}) = z0dinput.exp0d.(noms{k})(1:(end-1),:);
		end
	end
	
	z0dstruct.z0dinput = z0dinput;
	% calage du flux  au bord par definition au temps initial de la simulation
	z0dstruct.profil.psi  = z0dstruct.profil.psi -  z0dstruct.profil.psi(end, end) ;

elseif (z0dstruct.zs.temps(end) >= (time - 1e-9))

    % recherche des indices valides
    indval = min(length(z0dstruct.zs.temps),max(find(z0dstruct.zs.temps < (time - 1e-9))) + 1); 
    if isempty(indval) || (indval == 1)
        indval = [1,1,1];
    elseif indval  == 2
        indval = [1,1,2];
    else
        indval = 1:indval;
    end
    if length(indval) < nbte
        nbte = length(indval);
    end
	if length(indval)  ~= length(z0dstruct.zs.temps)
		noms = fieldnames(z0dstruct.zs);
		for k=1:length(noms)
			if size(z0dstruct.zs.(noms{k}),1) > 1
				z0dstruct.zs.(noms{k}) = z0dstruct.zs.(noms{k})(indval);
			end
		end	
		noms = fieldnames(z0dstruct.profil);
		for k=1:length(noms)
			if size(z0dstruct.profil.(noms{k}),1) > 1
				z0dstruct.profil.(noms{k}) = z0dstruct.profil.(noms{k})(indval,:);
			end
		end	
		noms = fieldnames(z0dstruct.z0dinput.cons);
		for k=1:length(noms)
			if size(z0dstruct.z0dinput.cons.(noms{k}),1) > 1
				z0dstruct.z0dinput.cons.(noms{k}) = z0dstruct.z0dinput.cons.(noms{k})(indval);
			end
		end	
		noms = fieldnames(z0dstruct.z0dinput.geo);
		for k=1:length(noms)
			if size(z0dstruct.z0dinput.geo.(noms{k}),1) > 1
				z0dstruct.z0dinput.geo.(noms{k}) = z0dstruct.z0dinput.geo.(noms{k})(indval);
			end
		end	
		noms = fieldnames(z0dstruct.z0dinput.exp0d);
		for k=1:length(noms)
			if size(z0dstruct.z0dinput.exp0d.(noms{k}),1) > 1
				z0dstruct.z0dinput.exp0d.(noms{k}) = z0dstruct.z0dinput.exp0d.(noms{k})(indval,:);
			end
		end	
                % retour en arriere		
		fprintf('R');	
		inloop = 0;

	else
		% meme temps => inloop = 1;
		fprintf('Z');
		inloop = 1;
    end
    tloc = z0dstruct.z0dinput.cons.temps;
    if any(diff(tloc)<= 0)
            dtloc_ref =  0.01 .* geo1t.a(1) .* geo1t.R(1);
            dtloc     = diff(tloc);
            dtloc(dtloc <= 0) = dtloc_ref;
            tloc = cumsum(cat(1,0,dtloc));
            if all(indval == 1)
                tloc = tloc - tloc(end) + z0dstruct.z0dinput.cons.temps(end);
            else
                tloc = tloc - tloc(end - 1) + z0dstruct.z0dinput.cons.temps(end -1);
            end
            tloc(end) = time;
            z0dstruct.zs.temps             = tloc;
            z0dstruct.profil.temps         = tloc;
            z0dstruct.z0dinput.cons.temps  = tloc;
            z0dstruct.z0dinput.exp0d.temps = tloc;
    end
   
	% temps suivant
	z0dinput 	= z0dstruct.z0dinput;
	%z0dinput.option = option;
	if ~isempty(option)
		noms = fieldnames(option);
		for k=1:length(noms)
			if isfield(z0dinput.option,noms{k})
				z0dinput.option.(noms{k}) = option.(noms{k});
			end
		end
	end
	if length(cons1t.temps) > 1
		indloc = find(cons1t.temps >= time,1);
		if isempty(indloc)
			indloc = 1;
		end
	else
		indloc = 1;
	end
	cons1t.temps(indloc)    = time;
	noms = fieldnames(cons1t);
	for k=1:length(noms)
		if isfield(z0dinput.cons,noms{k})
			z0dinput.cons.(noms{k})(end) = cons1t.(noms{k})(indloc);
		end
	end
	if ~isempty(geo1t)
		noms = fieldnames(geo1t);
		for k=1:length(noms)
			if isfield(z0dinput.geo,noms{k}) 
                		if ~isempty(geo1t.(noms{k}))
                    			z0dinput.geo.(noms{k})(end) = geo1t.(noms{k})(indloc);
                		end
			end
		end
	end

	
	% mise a jour de la structure experimentale vide
	noms = fieldnames(z0dinput.zsinfo);
	exp0d  = z0dinput.exp0d;
	texp   = z0dinput.cons.temps;
	%exp0d.temps = texp;
	vtnan  = texp .* NaN;
	for k = 1:length(noms)
		nomc = noms{k};
		if isfield(exp0d_in,nomc)
			var = getfield(exp0d_in,nomc);
			if ~isempty(var) && all(isfinite(exp0d_in.temps))
				warning off
				var = interp1(exp0d_in.temps,var,texp,'nearest');
				warning on
			else
				var = vtnan;
			end
			exp0d  = setfield(exp0d,nomc,var);
		else
			exp0d = setfield(exp0d,nomc,vtnan);
		end
	end
	exp0d.temps = texp;	
	% donnees experimentale
	z0dinput.exp0d = exp0d;
	
	
	% utilisation de la sepratrice si disponible
	if isfield(sepa1t,'Rsepa') && isfield(sepa1t,'Zsepa')
		if isfield(z0dinput.exp0d,'Rsepa')
			z0dinput.exp0d.Rsepa(end,:)  = sepa1t.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa(end, :) = sepa1t.Zsepa(indloc,:);
		else
			vtloc = ones(size(z0dinput.exp0d.temps));
			z0dinput.exp0d.Rsepa  = vtloc * sepa1t.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa  = vtloc * sepa1t.Zsepa(indloc,:);
		    fprintf('$')
		end
		fprintf('$')
	elseif isfield(exp0d_in,'Rsepa') && isfield(exp0d_in,'Zsepa')
		if isfield(z0dinput.exp0d,'Rsepa')
			z0dinput.exp0d.Rsepa(end,:) = exp0d_in.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa(end,:) = exp0d_in.Zsepa(indloc,:);
		else
			vtloc = ones(size(z0dinput.exp0d.temps));
			z0dinput.exp0d.Rsepa  = vtloc * exp0d_in.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa  = vtloc * exp0d_in.Zsepa(indloc,:);
		    fprintf('$')
		end
		fprintf('$')
	end

	% activation XDUR si disponible
	if isfield(xdur1t,'XDURx') && isfield(xdur1t,'XDURv')
		if size(xdur1t.XDURv,1) == 1
		   indx = 1;
		else
		    indx   = find(time <= xdur1t.XDURt,1);
		    if isempty(indx)
			    indx = length(xdur1t.XDURt);
		    end
		end
		if isfield(z0dinput.exp0d,'XDURv')
		    z0dinput.exp0d.XDURt(end) = time;
		    if size(xdur1t.XDURx,1) == 1
			  z0dinput.exp0d.XDURx(end,:) = xdur1t.XDURx(1,:);
		    else
			  z0dinput.exp0d.XDURx(end,:) = xdur1t.XDURx(indx,:);
	            end
		    z0dinput.exp0d.XDURv(end,:) = xdur1t.XDURv(indx,:);
		else
		    vtloc = ones(size(z0dinput.exp0d.temps));
		    z0dinput.exp0d.XDURt = z0dinput.exp0d.temps;
		    if size(xdur1t.XDURx,1) == 1
			  z0dinput.exp0d.XDURx = vtloc * xdur1t.XDURx(1,:);
		    else
			  z0dinput.exp0d.XDURx = vtloc * xdur1t.XDURx(indx,:);
	            end
		    z0dinput.exp0d.XDURv = vtloc * xdur1t.XDURv(indx,:);
		end
	elseif isfield(exp0d_in,'XDURx') & isfield(exp0d_in,'XDURv')
		indx   = find(time <= exp0d_in.XDURt,1);
		if isempty(indx)
			indx = 1;
		end
		if isfield(z0dinput.exp0d,'XDURv')
		    z0dinput.exp0d.XDURt(end) = time;
		    if size(exp0d_in.XDURx,1) == 1
			z0dinput.exp0d.XDURx(end,:) = exp0d_in.XDURx;
		    else
			z0dinput.exp0d.XDURx(end,:) = exp0d_in.XDURx(indx,:);
		    end
		    z0dinput.exp0d.XDURv(end,:) = exp0d_in.XDURv(indx,:);	
		else
		    vtloc = ones(size(z0dinput.exp0d.temps));
		    z0dinput.exp0d.XDURt = z0dinput.exp0d.temps;
		    if size(exp0d_in.XDURx,1) == 1
			  z0dinput.exp0d.XDURx = vtloc * exp0d_in.XDURx;
	            else
			  z0dinput.exp0d.XDURx = vtloc * exp0d_in.XDURx(indx,:);
		    end
		    z0dinput.exp0d.XDURv = vtloc * exp0d_in.XDURv(indx,:);	
		
		end
	end
	% 
	% la diffusion du courant est figee sur le temps 1:nbte
	%
% 	if inloop == 0
% 		z0dstruct.profil.psi(1:nbte,:) = ones(nbte,1) * z0dstruct.profil.psi(1,:);
% 	end
	% calcul  
	z0dinput.option.evolution = 1;
	z0dinput.option.init_geo  = 1;
	[last_zs,z0dstruct.info,last_profil] = ...
	      zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d, ...
	      [],z0dstruct.zs,z0dstruct.profil,nbte,inloop);
	z0dstruct.z0dinput = z0dinput;

%  	figure(168);
%  	plot(z0dstruct.profil.temps(end),z0dstruct.profil.psi(end,end),'or',last_profil.temps(end-1),last_profil.psi(end-1,end),'+b')
%  	hold on
	%mise en place du dernier temps
	noms = fieldnames(last_zs);
	for l=1:length(noms)
		nomc = noms{l};
		if size(last_zs.(nomc),1) > 1
			z0dstruct.zs.(nomc)(end) = last_zs.(nomc)(end-1);
		else
			z0dstruct.zs.(nomc) = last_zs.(nomc);		
		end
	end
	noms = fieldnames(last_profil);
	for l=1:length(noms)
		nomc = noms{l};
		if isfield(z0dstruct.profil,nomc)
			if size(last_profil.(nomc),1) > 1
				z0dstruct.profil.(nomc)(end,:) = last_profil.(nomc)(end-1,:);
			else
				z0dstruct.profil.(nomc) = last_profil.(nomc);		
			end
		end
	end
 	
 	switch option.dwdt_method

	case 'v4.2'

	    % les derivees temporelles de l'energie ne sont pas precise si caculees au vol
	    % calcul a posteriri pour plus de precision et moins de bruits
	    %pth                  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
	    %				  z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
	    switch z0dstruct.z0dinput.option.ploss_exp
	    case 'no_prad'
			pth  = z0dstruct.zs.pin;
	    case 'max_power'
			pth  =  max(eps,max(z0dstruct.profil.qe + z0dstruct.profil.qi,[],2));		
	    case 'max(pel)+max(pion)'
			pth  =  max(eps,max(z0dstruct.profil.qe,[],2) + max(z0dstruct.profil.qi,[],2));	
	    otherwise
	    	        pth  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
			       z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad + z0dstruct.zs.pioniz);
	    end
	    plim                 = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
	    % clacul de la derivee filtree
	    [wth,dwthdt]  = zdwdt0(z0dstruct.zs.wth,pth,z0dstruct.zs.taue,z0dstruct.zs.temps,z0dstruct.zs.vp);
	    %indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - mean(z0dstruct.zs.taue((end-nbte):end)))));
	    %indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - ...
	    %                    mean(z0dstruct.zs.wth((end-nbte):end) ./ max(1,z0dstruct.zs.pin((end-nbte):end))))));
	    %if isempty(indb)
	    %	indb = 1;
	    %elseif indb == length(z0dstruct.zs.temps)
	    %	indb = length(z0dstruct.zs.temps) - 1;
	    %end
	    %dwthdt        = (wth(end) - wth(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
	    dwthdt = dwthdt(end);
	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
	    z0dstruct.zs.dwthdt(end) = dwthdt;
	    z0dstruct.zs.wth(end)   = wth(end);
	    % derivee de W
	    z0dstruct.zs.w(end)       = wth(end) + (z0dstruct.zs.esup_fus(end) + z0dstruct.zs.esup_icrh(end) +  ...
						    real(z0dstruct.zs.esup_nbi(end)) + imag(z0dstruct.zs.esup_nbi(end)) + z0dstruct.zs.esup_lh(end));
	    %z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
	    z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(end - 1)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(end - 1));

	case {'old','explicit','mixed','freebie'}

	    %pth    = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
	    switch z0dstruct.z0dinput.option.ploss_exp
	    case 'no_prad'
			pth  = z0dstruct.zs.pin;
	    case 'max_power'
			pth  =  max(eps,max(z0dstruct.profil.qe + z0dstruct.profil.qi,[],2));		
	    case 'max(pel)+max(pion)'
			pth  =  max(eps,max(z0dstruct.profil.qe,[],2) + max(z0dstruct.profil.qi,[],2));	
	    otherwise
	    	        pth  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
			       z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad + z0dstruct.zs.pioniz);
	    end
	    plim   = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
	    wth           = last_zs.wth;
	    wth(1:2)      = z0dstruct.zs.wth(end-1:end);
	    [wth,dwthdt]  = zdwdt0(wth,last_zs.pth,last_zs.taue,last_zs.temps,last_zs.vp);
	    dwthdt(1:2) = z0dstruct.zs.dwthdt(end-1:end);
	    %figure(21);clf;plot(last_zs.temps,wth,'r',last_zs.temps,dwthdt,'b');hold on
	    pp          = polyfit(last_zs.temps(2:end),dwthdt(2:end),1);
	    dwthdt      = polyval(pp,z0dstruct.zs.temps(end)); 
	    %pp          = polyfit(last_zs.temps(2:end),wth(2:end),2);
	    %dwthdt      = polyval(polyder(pp),z0dstruct.zs.temps(end)); 
	    %plot(z0dstruct.zs.temps(end),dwthdt,'ob');drawnow
	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
	    z0dstruct.zs.dwthdt(end) = dwthdt;
	    esup = last_zs.esup_fus + last_zs.esup_icrh + real(last_zs.esup_nbi) + imag(last_zs.esup_nbi) + last_zs.esup_lh;
	    desupdt = pdederive(last_zs.temps,esup,2,2,1,1);
	    z0dstruct.zs.dwdt(end)    =  z0dstruct.zs.dwthdt(end) + desupdt(end-1);

%  	case 'explicit'
%  
%  	    pth    = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
%  	    plim   = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
%  	    [wth,dwthdt]  = zdwdt0(last_zs.wth,last_zs.pth,last_zs.taue,last_zs.temps,last_zs.vp);
%  	    pp   = polyfit(last_zs.temps,wth,2);
%              dwthdt = polyval(polyder(pp),z0dstruct.zs.temps(end-1)); 
%  	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
%   	    z0dstruct.zs.dwthdt(end) = dwthdt;
%  	    esup = last_zs.esup_fus + last_zs.esup_icrh + real(last_zs.esup_nbi) + imag(last_zs.esup_nbi) + last_zs.esup_lh;
%  	    pp   = polyfit(last_zs.temps,esup,2);
%  	    desupdt = polyval(polyder(pp),z0dstruct.zs.temps(end-1)); 
%  	    z0dstruct.zs.dwdt(end)    =  z0dstruct.zs.dwthdt(end) + desupdt;

       end

else
	% temps suivant
	z0dinput 	= z0dstruct.z0dinput;
	%z0dinput.option = option;
	if ~isempty(option)
		noms = fieldnames(option);
		for k=1:length(noms)
			if isfield(z0dinput.option,noms{k})
				z0dinput.option.(noms{k}) = option.(noms{k});
			end
		end
	end
	if length(cons1t.temps) > 1
		indloc = find(cons1t.temps >= time,1);
		if isempty(indloc)
			indloc = 1;
		end
	else
		indloc = 1;
	end
	cons1t.temps(indloc)    = time;
	noms = fieldnames(cons1t);
	for k=1:length(noms)
		if isfield(z0dinput.cons,noms{k})
			z0dinput.cons.(noms{k})(end + 1) = cons1t.(noms{k})(indloc);
		end
	end
	if ~isempty(geo1t)
		noms = fieldnames(geo1t);
		for k=1:length(noms)
			if isfield(z0dinput.geo,noms{k}) 
                if ~isempty(geo1t.(noms{k}))
                    z0dinput.geo.(noms{k})(end + 1) = geo1t.(noms{k})(indloc);
                end
			end
		end
	end

	
	% mise a jour de la structure experimentale vide
	noms = fieldnames(z0dinput.zsinfo);
	exp0d  = z0dinput.exp0d;
	texp   = z0dinput.cons.temps;
	%exp0d.temps = texp;
	vtnan  = texp .* NaN;
	for k = 1:length(noms)
		nomc = noms{k};
		if isfield(exp0d_in,nomc)
			var = getfield(exp0d_in,nomc);
			if ~isempty(var) && all(isfinite(exp0d_in.temps))
				warning off
				var = interp1(exp0d_in.temps,var,texp,'nearest');
				warning on
			else
				var = vtnan;
			end
			exp0d  = setfield(exp0d,nomc,var);
		else
			exp0d = setfield(exp0d,nomc,vtnan);
		end
	end
	exp0d.temps = texp;	
	% donnees experimentale
	z0dinput.exp0d = exp0d;
	
	
	% utilisation de la sepratrice si disponible
%  	if isfield(sepa1t,'Rsepa') & isfield(sepa1t,'Zsepa')
%  		z0dinput.exp0d.Rsepa(end + 1,:) = sepa1t.Rsepa(indloc,:);
%  		z0dinput.exp0d.Zsepa(end + 1, :) = sepa1t.Zsepa(indloc,:);
%  		fprintf('$')
%  	elseif isfield(exp0d_in,'Rsepa') & isfield(exp0d_in,'Zsepa')
%  		z0dinput.exp0d.Rsepa(end + 1,:) = exp0d_in.Rsepa(indloc,:);
%  		z0dinput.exp0d.Zsepa(end + 1,:) = exp0d_in.Zsepa(indloc,:);
%  		fprintf('$')
%  	end
	if isfield(sepa1t,'Rsepa') && isfield(sepa1t,'Zsepa')
		if isfield(z0dinput.exp0d,'Rsepa')
			z0dinput.exp0d.Rsepa(end + 1,:)  = sepa1t.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa(end + 1, :) = sepa1t.Zsepa(indloc,:);
		else
			vtloc = ones(size(z0dinput.exp0d.temps));
			z0dinput.exp0d.Rsepa  = vtloc * sepa1t.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa  = vtloc * sepa1t.Zsepa(indloc,:);
		    fprintf('$')
		end
		fprintf('$')
	elseif isfield(exp0d_in,'Rsepa') && isfield(exp0d_in,'Zsepa')
		if isfield(z0dinput.exp0d,'Rsepa')
			z0dinput.exp0d.Rsepa(end + 1,:) = exp0d_in.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa(end + 1,:) = exp0d_in.Zsepa(indloc,:);
		else
			vtloc = ones(size(z0dinput.exp0d.temps));
			z0dinput.exp0d.Rsepa  = vtloc * exp0d_in.Rsepa(indloc,:);
			z0dinput.exp0d.Zsepa  = vtloc * exp0d_in.Zsepa(indloc,:);
		    fprintf('$')
		end
		fprintf('$')
	end

%  	% activation XDUR si disponible
%  	if isfield(xdur1t,'XDURx') & isfield(xdur1t,'XDURv')
%  		z0dinput.exp0d.XDURt(end + 1 ) = time;
%  		z0dinput.exp0d.XDURx(end + 1,:) = xdur1t.XDURx(indloc,:);
%  		z0dinput.exp0d.XDURv(end + 1,:) = xdur1t.XDURv(indloc,:);	
%  	elseif isfield(exp0d_in,'XDURx') & isfield(exp0d_in,'XDURv')
%  		z0dinput.exp0d.XDURt(end + 1) = time;
%  		indx   = find(time <= exp0d_in.XDURt,1);
%  		if isempty(indx)
%  			indx = 1;
%  		end
%  		z0dinput.exp0d.XDURx = exp0d_in.XDURx;
%  		z0dinput.exp0d.XDURv(end + 1,:) = exp0d_in.XDURv(indx,:);	
%  			
%  	end
	% activation XDUR si disponible
	if isfield(xdur1t,'XDURx') && isfield(xdur1t,'XDURv')
		if size(xdur1t.XDURv,1) == 1
		   indx = 1;
		else
		    indx   = find(time <= xdur1t.XDURt,1);
		    if isempty(indx)
			    indx = length(xdur1t.XDURt);
		    end
		end
		if isfield(z0dinput.exp0d,'XDURv')
		    z0dinput.exp0d.XDURt(end + 1) = time;
		    if size(xdur1t.XDURx,1) == 1
			  z0dinput.exp0d.XDURx(end + 1,:) = xdur1t.XDURx(1,:);
		    else
			  z0dinput.exp0d.XDURx(end + 1,:) = xdur1t.XDURx(indx,:);
	            end
		    z0dinput.exp0d.XDURv(end + 1,:) = xdur1t.XDURv(indx,:);
		else
		    vtloc = ones(size(z0dinput.exp0d.temps));
		    z0dinput.exp0d.XDURt = z0dinput.exp0d.temps;
		    if size(xdur1t.XDURx,1) == 1
			  z0dinput.exp0d.XDURx = vtloc * xdur1t.XDURx(1,:);
		    else
			  z0dinput.exp0d.XDURx = vtloc * xdur1t.XDURx(indx,:);
	            end
		    z0dinput.exp0d.XDURv = vtloc * xdur1t.XDURv(indx,:);
		end
	elseif isfield(exp0d_in,'XDURx') && isfield(exp0d_in,'XDURv')
		indx   = find(time <= exp0d_in.XDURt,1);
		if isempty(indx)
			indx = 1;
		end
		if isfield(z0dinput.exp0d,'XDURv')
		    z0dinput.exp0d.XDURt(end) = time;
		    if size(exp0d_in.XDURx,1) == 1
			z0dinput.exp0d.XDURx(end + 1,:) = exp0d_in.XDURx;
		    else
			z0dinput.exp0d.XDURx(end + 1,:) = exp0d_in.XDURx(indx,:);
		    end
		    z0dinput.exp0d.XDURv(end + 1,:) = exp0d_in.XDURv(indx,:);	
		else
		    vtloc = ones(size(z0dinput.exp0d.temps));
		    z0dinput.exp0d.XDURt = z0dinput.exp0d.temps;
		    if size(exp0d_in.XDURx,1) == 1
			  z0dinput.exp0d.XDURx = vtloc * exp0d_in.XDURx;
	            else
			  z0dinput.exp0d.XDURx = vtloc * exp0d_in.XDURx(indx,:);
		    end
		    z0dinput.exp0d.XDURv = vtloc * exp0d_in.XDURv(indx,:);	
		
		end
	end
	% 
	% la diffusion du courant est figee sur le temps 1:nbte
	%
	%z0dstruct.profil.psi(1:nbte,:) = ones(nbte,1) * z0dstruct.profil.psi(1,:);
	% calcul  
	z0dinput.option.evolution = 1;
	z0dinput.option.init_geo  = 1;
	[last_zs,z0dstruct.info,last_profil] = ...
	      zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d, ...
	      [],z0dstruct.zs,z0dstruct.profil,nbte);
	z0dstruct.z0dinput = z0dinput;

%  	figure(168);
%  	plot(z0dstruct.profil.temps(end),z0dstruct.profil.psi(end,end),'or',last_profil.temps(end-1),last_profil.psi(end-1,end),'+b')
%  	hold on
	
	%mise enplace du dernier temps
	noms = fieldnames(last_zs);
	for l=1:length(noms)
		nomc = noms{l};
		if size(last_zs.(nomc),1) > 1
			z0dstruct.zs.(nomc)(end + 1) = last_zs.(nomc)(end-1);
		else
			z0dstruct.zs.(nomc) = last_zs.(nomc);		
		end
	end
	noms = fieldnames(last_profil);
	for l=1:length(noms)
		nomc = noms{l};
		if isfield(z0dstruct.profil,nomc)
			if size(last_profil.(nomc),1) > 1
				z0dstruct.profil.(nomc)(end + 1,:) = last_profil.(nomc)(end-1,:);
			else
				z0dstruct.profil.(nomc) = last_profil.(nomc);		
			end
		end
	end
  		
%  	% les derivees temporelles de l'energie ne sont pas precise si caculees au vol
%  	% calcul a posteriri pour plus de precision et moins de bruits
%  	pth                  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
%    				       z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
%    	plim                 = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
%    	% clacul de la derivee filtree
%  	[wth,dwthdt]  = zdwdt0(z0dstruct.zs.wth,pth,z0dstruct.zs.taue,z0dstruct.zs.temps,z0dstruct.zs.vp);
%  	%indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - mean(z0dstruct.zs.taue((end-nbte):end)))));
%  	%indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - ...
%  	%                    mean(z0dstruct.zs.wth((end-nbte):end) ./ max(1,z0dstruct.zs.pin((end-nbte):end))))));
%  	%if isempty(indb)
%  	%	indb = 1;
%  	%elseif indb == length(z0dstruct.zs.temps)
%  	%	indb = length(z0dstruct.zs.temps) - 1;	
%  	%end
%          %dwthdt        = (wth(end) - wth(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
%    	dwthdt = dwthdt(end);
%  	dwthdt = max(-plim(end),min(plim(end),dwthdt));
%  	z0dstruct.zs.dwthdt(end) = dwthdt;
%  	z0dstruct.zs.wth(end)   = wth(end);
%  	% derivee de W	
%   	z0dstruct.zs.w(end)       = wth(end) + (z0dstruct.zs.esup_fus(end) + z0dstruct.zs.esup_icrh(end) +  ...
%  	                                        real(z0dstruct.zs.esup_nbi(end)) + imag(z0dstruct.zs.esup_nbi(end)) + z0dstruct.zs.esup_lh(end));
%   	%z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
%   	z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(end - 1)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(end - 1));

 	switch option.dwdt_method

	case 'v4.2'

	    % les derivees temporelles de l'energie ne sont pas precise si caculees au vol
	    % calcul a posteriri pour plus de precision et moins de bruits
	    %pth                  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
	    %				  z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
	    switch z0dstruct.z0dinput.option.ploss_exp
	    case 'no_prad'
			pth  = z0dstruct.zs.pin;
	    case 'max_power'
			pth  =  max(eps,max(z0dstruct.profil.qe + z0dstruct.profil.qi,[],2));		
	    case 'max(pel)+max(pion)'
			pth  =  max(eps,max(z0dstruct.profil.qe,[],2) + max(z0dstruct.profil.qi,[],2));	
	    otherwise
	    	        pth  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
			       z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad + z0dstruct.zs.pioniz);
	    end
	    plim                 = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
	    % clacul de la derivee filtree
	    [wth,dwthdt]  = zdwdt0(z0dstruct.zs.wth,pth,z0dstruct.zs.taue,z0dstruct.zs.temps,z0dstruct.zs.vp);
	    %indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - mean(z0dstruct.zs.taue((end-nbte):end)))));
	    %indb          = max(find(z0dstruct.zs.temps <= (z0dstruct.zs.temps(end) - ...
	    %                    mean(z0dstruct.zs.wth((end-nbte):end) ./ max(1,z0dstruct.zs.pin((end-nbte):end))))));
	    %if isempty(indb)
	    %	indb = 1;
	    %elseif indb == length(z0dstruct.zs.temps)
	    %	indb = length(z0dstruct.zs.temps) - 1;
	    %end
	    %dwthdt        = (wth(end) - wth(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
	    dwthdt = dwthdt(end);
	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
	    z0dstruct.zs.dwthdt(end) = dwthdt;
	    z0dstruct.zs.wth(end)   = wth(end);
	    % derivee de W
	    z0dstruct.zs.w(end)       = wth(end) + (z0dstruct.zs.esup_fus(end) + z0dstruct.zs.esup_icrh(end) +  ...
						    real(z0dstruct.zs.esup_nbi(end)) + imag(z0dstruct.zs.esup_nbi(end)) + z0dstruct.zs.esup_lh(end));
	    %z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(indb)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(indb));
	    z0dstruct.zs.dwdt(end)    = (z0dstruct.zs.w(end) - z0dstruct.zs.w(end - 1)) ./ (z0dstruct.zs.temps(end) - z0dstruct.zs.temps(end - 1));

	case {'old','explicit','mixed','freebie'}

	    %pth    = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
	    switch z0dstruct.z0dinput.option.ploss_exp
	    case 'no_prad'
			pth  = z0dstruct.zs.pin;
	    case 'max_power'
			pth  =  max(eps,max(z0dstruct.profil.qe + z0dstruct.profil.qi,[],2));		
	    case 'max(pel)+max(pion)'
			pth  =  max(eps,max(z0dstruct.profil.qe,[],2) + max(z0dstruct.profil.qi,[],2));	
	    otherwise
	    	        pth  = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + ...
			       z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad + z0dstruct.zs.pioniz);
	    end
	    plim   = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
	    wth           = last_zs.wth;
	    wth(1:2)      = z0dstruct.zs.wth(end-1:end);
	    [wth,dwthdt]  = zdwdt0(wth,last_zs.pth,last_zs.taue,last_zs.temps,last_zs.vp);
	    dwthdt(1:2) = z0dstruct.zs.dwthdt(end-1:end);
	    %figure(21);clf;plot(last_zs.temps,wth,'r',last_zs.temps,dwthdt,'b');hold on
	    pp          = polyfit(last_zs.temps(2:end),dwthdt(2:end),1);
	    dwthdt      = polyval(pp,z0dstruct.zs.temps(end)); 
	    %pp          = polyfit(last_zs.temps(2:end),wth(2:end),2);
	    %dwthdt      = polyval(polyder(pp),z0dstruct.zs.temps(end)); 
	    %plot(z0dstruct.zs.temps(end),dwthdt,'ob');drawnow
	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
	    z0dstruct.zs.dwthdt(end) = dwthdt;
	    esup = last_zs.esup_fus + last_zs.esup_icrh + real(last_zs.esup_nbi) + imag(last_zs.esup_nbi) + last_zs.esup_lh;
	    desupdt = pdederive(last_zs.temps,esup,2,2,1,1);
	    z0dstruct.zs.dwdt(end)    =  z0dstruct.zs.dwthdt(end) + desupdt(end-1);

%  	case 'explicit'
%  
%  	    pth    = z0dstruct.zs.pin  - (z0dstruct.zs.pbrem + z0dstruct.zs.pcyclo + z0dstruct.z0dinput.option.fprad .* z0dstruct.zs.prad);
%  	    plim   = max(1,min(z0dstruct.zs.vp .* 1e6, max(0.8 .* pth, 0.2 .* abs(z0dstruct.zs.pin))));
%  	    [wth,dwthdt]  = zdwdt0(last_zs.wth,last_zs.pth,last_zs.taue,last_zs.temps,last_zs.vp);
%  	    pp   = polyfit(last_zs.temps,wth,1);
%              dwthdt = polyval(polyder(pp),z0dstruct.zs.temps(end-1)); 
%  	    dwthdt = max(-plim(end),min(plim(end),dwthdt));
%   	    z0dstruct.zs.dwthdt(end) = dwthdt;
%  	    esup = last_zs.esup_fus + last_zs.esup_icrh + real(last_zs.esup_nbi) + imag(last_zs.esup_nbi) + last_zs.esup_lh;
%  	    pp   = polyfit(last_zs.temps,esup,1);
%  	    desupdt = polyval(polyder(pp),z0dstruct.zs.temps(end-1)); 
%  	    z0dstruct.zs.dwdt(end)    =  z0dstruct.zs.dwthdt(end) + desupdt;

         end

end

% sorties
noms = fieldnames(z0dstruct.zs);
zs = [];
for k=1:length(noms)
	if ~isfield(z0dstruct.zs,noms{k})
		% rien
	else
		zs.(noms{k}) = z0dstruct.zs.(noms{k})(end);		
	end
end
noms = fieldnames(z0dstruct.profil);
for k=1:length(noms)
		profil.(noms{k}) = z0dstruct.profil.(noms{k})(end,:);		
end
