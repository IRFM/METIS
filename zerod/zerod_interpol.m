% ZEROD_INTERPOL : simule zerodevolution en interpolant les donnees du workspace (structure post)
%------------------------------------------------------------------------------------------------
% fonction Matlab 2012b; zerod_interpol.m -> zerod_interpol
%
% Cette fonction  simule zerodevolution en interpolant les donnees du workspace (structure post)
%
% syntaxe : 
%
%   initialisation :       
%       [zs,profil,z0dstruct] = zerod_interpol([],option,temps_ini,cons,geo,[exp0d,sepa,xdur]);
%
%   calcul :
%       [zs,profil,z0dstruct] = zerod_interpol(z0dstruct,option,temps,cons,geo,[exp0d,sepa,xdur]);
%
% entrees :
%
%     z0dstruct               = structure de calcul au temps precedents 
%     option                  = structure d'option de METIS (z0dinput.option)
%     temps                   = temps pour l'interpolation
%
%  ---- optionnelles ----
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
%      [option,cons1t,geo1t] = zerod_interpol;
%      [zs,profil,z0dstruct] = zerod_interpol([],z0dinput.option,1);
%      [zs,profil,z0dstruct] = zerod_interpol(z0dstruct,z0dinput.option,2);
%
% fonction ecrite par J-F Artaud
% version CVS (creation le 28/09/2015)
%-----------------------------------------------------------------------
%
function [zs,profil,z0dstruct] = zerod_interpol(z0dstruct,option,time,cons1t,geo1t,exp0d_in,sepa1t,xdur1t)

% lecture des donnees dans le workspace
post            = evalin('base','post');
post.z0dinput.geo.temps = post.z0dinput.cons.temps;

% test des entrees
if nargin == 0
	z0dinput 	= zerod_interp1t(post.z0dinput,min(post.z0dinput.cons.temps));
	zs              = z0dinput.option;
	profil          = z0dinput.cons;
	z0dstruct       = z0dinput.geo;
	noms = fieldnames(z0dstruct);
	for k=1:length(noms)
		if isempty(z0dstruct.(noms{k}))
			z0dstruct = rmfield(z0dstruct,noms{k});
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

if nargin < 4
      cons1t = zerod_interp1t(post.z0dinput.cons,time);
end
if nargin < 5
      geo1t = zerod_interp1t(post.z0dinput.geo,time);
end

if ~isfield(option,'ploss_exp')
	option.ploss_exp ='with_prad';
end
% nombre de temps interne 
nbte = 3;

% selon la structure d'entree
% initialisation
if isempty(z0dstruct)
	% initilisation
        t0              = time - max(1e-6,min(1,0.01 .* geo1t.a(1) .* geo1t.R(1)));
	temps    	= linspace(t0,time,nbte)';
	temps           = cat(1,temps,2 .* temps(end) - temps(end-1));
	vt              = ones(size(temps));
	z0dinput 	= zerod_interp1t(post.z0dinput,temps);
	if ~isempty(option)
		noms = fieldnames(option);
		for k=1:length(noms)
			if isfield(z0dinput.option,noms{k})
				z0dinput.option.(noms{k}) = option.(noms{k});
			end
		end
	end
	z0dinput.option = option;
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
	fprintf('Ii-')
	z0dstruct.zs       = zerod_interp1t(post.zerod,temps);
	z0dstruct.info     = post.z0dinput.zsinfo;
	z0dstruct.profil   = zerod_interp1t(post.profil0d,temps);
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
		fprintf('Ri-');	
		inloop = 0;

    else
		% meme temps => inloop = 1;
		fprintf('Si-');
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
	z0dinput.option = option;
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
			if ~isempty(var)
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
	if isfield(sepa1t,'Rsepa') & isfield(sepa1t,'Zsepa')
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
	elseif isfield(exp0d_in,'Rsepa') & isfield(exp0d_in,'Zsepa')
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
	if isfield(xdur1t,'XDURx') & isfield(xdur1t,'XDURv')
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

        % calcul  
	z0dinput.option.evolution = 1;
	z0dinput.option.init_geo  = 1;
	last_zs     = zerod_interp1t(post.zerod,[time - sqrt(eps),time]');
	last_profil = zerod_interp1t(post.profil0d,[time - sqrt(eps),time]');
	z0dstruct.z0dinput = z0dinput;
	
	% cas de la separatrice
	if isfield(last_profil,'Rsepa') 
	    if ~isfield(z0dstruct.profil,'Rsepa') 
		  z0dstruct.profil.Rsepa = ones(size(z0dstruct.profil.temps)) * last_profil.Rsepa(end,:);
		  z0dstruct.profil.Zsepa = ones(size(z0dstruct.profil.temps)) * last_profil.Zsepa(end,:);
	    else
                  nbp           = size(z0dstruct.profil.Zsepa,2);
                  [last_profil.Rsepa,last_profil.Zsepa] = resample_sepa(last_profil.Rsepa,last_profil.Zsepa,nbp);
	    end	
	end

	%mise en place du dernier temps
	noms = fieldnames(last_zs);
	for l=1:length(noms)
		nomc = noms{l};
		if size(last_zs.(nomc),1) > 1
			z0dstruct.zs.(nomc)(end) = last_zs.(nomc)(end);
		else
			z0dstruct.zs.(nomc) = last_zs.(nomc);		
		end
	end
	noms = fieldnames(last_profil);
	for l=1:length(noms)
		nomc = noms{l};
		if isfield(z0dstruct.profil,nomc)
			if size(last_profil.(nomc),1) > 1
				z0dstruct.profil.(nomc)(end,:) = last_profil.(nomc)(end,:);
			else
				z0dstruct.profil.(nomc) = last_profil.(nomc);		
			end
		end
	end
 	
else
	% temps suivant
	z0dinput 	= z0dstruct.z0dinput;
	z0dinput.option = option;
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
	
	if isfield(sepa1t,'Rsepa') & isfield(sepa1t,'Zsepa')
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
	elseif isfield(exp0d_in,'Rsepa') & isfield(exp0d_in,'Zsepa')
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

	% activation XDUR si disponible
	if isfield(xdur1t,'XDURx') & isfield(xdur1t,'XDURv')
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
	elseif isfield(exp0d_in,'XDURx') & isfield(exp0d_in,'XDURv')
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
	% calcul  
	z0dinput.option.evolution = 1;
	z0dinput.option.init_geo  = 1;
	last_zs     = zerod_interp1t(post.zerod,[time - sqrt(eps),time]');
	last_profil = zerod_interp1t(post.profil0d,[time - sqrt(eps),time]');
	z0dstruct.z0dinput = z0dinput;
	fprintf('Ti-');	
	
	% cas de la separatrice
	if isfield(last_profil,'Rsepa') 
	    if ~isfield(z0dstruct.profil,'Rsepa') 
		  z0dstruct.profil.Rsepa = ones(size(z0dstruct.profil.temps)) * last_profil.Rsepa(end,:);
		  z0dstruct.profil.Zsepa = ones(size(z0dstruct.profil.temps)) * last_profil.Zsepa(end,:);
	    else
                  nbp           = size(z0dstruct.profil.Zsepa,2);
                  [last_profil.Rsepa,last_profil.Zsepa] = resample_sepa(last_profil.Rsepa,last_profil.Zsepa,nbp);
	    end	
	end

	%mise enplace du dernier temps
	noms = fieldnames(last_zs);
	for l=1:length(noms)
		nomc = noms{l};
		if size(last_zs.(nomc),1) > 1
			z0dstruct.zs.(nomc)(end + 1) = last_zs.(nomc)(end);
		else
			z0dstruct.zs.(nomc) = last_zs.(nomc);		
		end
	end
	noms = fieldnames(last_profil);
	for l=1:length(noms)
		nomc = noms{l};
		if isfield(z0dstruct.profil,nomc)
			if size(last_profil.(nomc),1) > 1
				z0dstruct.profil.(nomc)(end + 1,:) = last_profil.(nomc)(end,:);
			else
				z0dstruct.profil.(nomc) = last_profil.(nomc);		
			end
		end
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



% interpolation 1 temps d'une structure
function struct_out = zerod_interp1t(struct_in,temps)

if ~isstruct(struct_in)
    if isempty(struct_in)
        struct_out = struct_in;        
    elseif isnumeric(struct_in)
	if size(struct_in,1) == 1
	    struct_out = struct_in;
	else
	    struct_out = struct_in(end)  .* ones(size(temps));
	end
    elseif iscell(struct_in)
         struct_out = struct_in{end};
    else
        struct_out = struct_in;    
    end
else
    noms = fieldnames(struct_in);
    if isfield(struct_in,'temps')
	temps_loc = struct_in.temps;
    else
	temps_loc = [];
    end
    for k = 1:length(noms)
	if isstruct(struct_in.(noms{k}))
	    struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps); 
	elseif isempty(temps_loc)
	    struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps); 
	elseif length(struct_in.temps) == size(struct_in.(noms{k}),1)
	    struct_out.(noms{k}) = interp1(temps_loc,struct_in.(noms{k}),temps,'linear',NaN);
	    near = interp1(temps_loc,struct_in.(noms{k}),temps,'nearest','extrap');
	    struct_out.(noms{k})(~isfinite(struct_out.(noms{k}))) = near(~isfinite(struct_out.(noms{k})));
	else
	    struct_out.(noms{k})  = zerod_interp1t(struct_in.(noms{k}),temps); 
	end
    end
    if ~isempty(temps_loc)
      struct_out.temps = temps;
    end
 end
 
 
 
function [Rsepa,Zsepa] = resample_sepa(R,Z,nbp)

% reservation memoire
Rsepa = NaN * ones(size(R,1),nbp);
Zsepa = NaN * ones(size(Z,1),nbp);

% mise sur les memes angles
Rc   = (max(R,[],2) + min(R,[],2)) ./ 2;
Zc   = (max(Z,[],2) + min(Z,[],2)) ./ 2; 
th   = linspace(0,2.*pi,nbp);
vh   = ones(1,size(R,2));
vt   = ones(size(R,1),1);
cx   = (R - Rc*vh) + sqrt(-1) .* (Z-Zc*vh);
thx  = unwrap(angle(cx),[],2);
thx = angle(cx);
thx = thx .* (thx>=0) + (2 .* pi + thx) .* (thx <0);
rhox = abs(cx);
thx(thx<0) = thx(thx<0) + 2 .* pi;
for k = 1:size(thx,1)
	[thk,indice] = unique(thx(k,:));
	rhok         = rhox(k,indice);
	[thk,indice] = sort(thk); 
	rhok         = rhok(indice);
	rhok = cat(2,rhok,rhok(:,2:end-1),rhok);
	thk = cat(2,thk -2.*pi,thk(:,2:end-1),thk+2.*pi);
	ind_bad = find(diff(thk) <= 0);
	if ~isempty(ind_bad)
	    thk(ind_bad) = [];
	    rhok(ind_bad) = [];
	    if length(thk) < 3
		error('no enougth distinct points to resample LCFS');
	    end
	end
	rhoo   = pchip(thk,rhok,th);
	R      = Rc(k)  + rhoo .* cos(th);
	Z      = Zc(k)  + rhoo .* sin(th);
	Rsepa(k,:) = R;
	Zsepa(k,:) = Z;
end
