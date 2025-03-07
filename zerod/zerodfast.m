% calcul rapide avec extraction de point pertinent
function [zs,info,profli]= zerodfast(option,cons,geo,exp0d,thyb,last_zs,last_profil,nbte,inloop)

if nargin < 5
	thyb = [];
end

if isfield(option,'tol0d')
	tol0d = option.tol0d;
else 
	tol0d = 0;
end

%MATLAB:interp1:NaNinY
warning off
% compatibilite
if ~isfield(cons,'xece') & isfield(option,'xece')
 	cons.xece = option.xece .* ones(size(cons.temps));
end

if nargin <= 8
	inloop = 0;	
elseif isempty(inloop)
	inloop = 0;
end

% handling external equilibrium
% when double call of zerod (fast and then full)
if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
    if isappdata(0,'METIS_EXTERNAL_CURDIF_EQUI')
        equi_ext_full = getappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
    else
        equi_ext_full = [];
    end
else
    equi_ext_full = [];
end
    


if inloop == 1

       % pour le mode evolution
	times = cons.temps;	
	if length(cons.temps) >= nbte
		temps = cons.temps(end-nbte+1:end);
		indsel  = (length(cons.temps) - nbte + 1):length(cons.temps);
	else
		temps = cons.temps;	
	        indsel  = 1:length(cons.temps);
	end		
	% indication de reduction
	fprintf('L %g ->',temps(end));		
	% modification des donnees
	noms = fieldnames(cons);
	for l=1:length(noms)
		nomc = noms{l};
		val  = cons.(nomc);
		valn = val(indsel);
		valn(end+1) = val(end);
		cons.(nomc) = valn;
	end
        % maximum tau_E = 30 s
	cons.temps(end) = min( cons.temps(end - 1) + 30, ...
                          max(cons.temps(end - 1) + last_zs.wth(end-1) ./last_zs.pin(end-1),  ...
                          max(cons.temps(end - 1) + last_zs.wth(end) ./last_zs.pin(end),  ...
                          2 .* cons.temps(end - 1) - cons.temps(end-2))));


	noms = fieldnames(geo);
	for l=1:length(noms)
   		nomc = noms{l};
   		val  = geo.(nomc);
   		if ~isempty(val)
      			valn  = val(indsel);
 			valn(end+1) = val(end);
     			geo.(nomc) = valn;
   		end
	end

	noms = fieldnames(exp0d);
	tref = exp0d.temps;
	for l=1:length(noms)
   		nomc = noms{l};
   		val  = exp0d.(nomc);
  		if ~isempty(val) && (size(val,1) > 1)
			tloc = tref;
			indok = find(isfinite(tloc) & all(isfinite(val),2));
			val = val(indok,:);
			tloc = tloc(indok);
   			if length(tloc) >= 2;
      				if (size(val,1) == length(times(indok)))
					warning off
	         			valn  = interp1(tloc,val,cons.temps,'linear');		
					warning on
      	      				indbad      = find(any(~isfinite(valn),2));
     	      				if ~isempty(indbad)
	             				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
              				end
      	      				exp0d.(nomc) = valn;
      				end
			else
				exp0d.(nomc) = NaN .* ones(size(cons.temps,1),size(exp0d.(nomc),2));
   			end
		end
	end
	if isfield(exp0d,'XDURt')
		exp0d.XDURt = cons.temps;
	end

        % pour le mode evolution
	noms = fieldnames(last_zs);
	tref = last_zs.temps;
	for l=1:length(noms)
		nomc = noms{l};
		val  = last_zs.(nomc);
		if length(val) > 1
			%valn  = interp1(tref,val,cons.temps,'nearest');	
			valn  = val((end-nbte+1):end);	
			valn(end + 1)     = val(end);
			last_zs.(nomc) = valn;
		end
	end
	last_zs.temps = cons.temps;
        if option.nb_nbi == 2
            last_zs.frnbi(end-1:end) = 1 + sqrt(-1);
        else
            last_zs.frnbi(end-1:end) = 1;
        end
    
	noms = fieldnames(last_profil);
	tref = last_profil.temps;
	for l=1:length(noms)
		nomc = noms{l};
		val  = last_profil.(nomc);
		if size(val,1) > 1
			%valn = interp1(tref,val,cons.temps,'nearest');	
			valn  = val((end-nbte+1):end,:);	
			valn(end + 1,:) = val(end,:);
			last_profil.(nomc) = valn;
		end
	end
	last_profil.temps  = cons.temps;
		
	% appel du 0d
	warning on
        %% option.short_dt = 'off';
        option.tswitch = Inf;
        if tol0d > 0
 		[zs,info,profli]= zerod(option,cons,geo,exp0d,tol0d,thyb,last_zs,last_profil);
        else
		[zs,info,profli]= zerod(option,cons,geo,exp0d,1e-3,thyb,last_zs,last_profil);
        end

elseif nargin >= 8

        % pour le mode evolution
	times = cons.temps;	
	if length(cons.temps) >= nbte
		temps = cons.temps(end-nbte+1:end);
		indsel  = (length(cons.temps) - nbte + 1):length(cons.temps);
	else
		temps = cons.temps;	
	        indsel  = 1:length(cons.temps);
	end		
	% indication de reduction
	fprintf('E %g ->',temps(end));		
	% modification des donnees
	noms = fieldnames(cons);
	for l=1:length(noms)
		nomc = noms{l};
		val  = cons.(nomc);
		valn = val(indsel);
		valn(end+1) = val(end);
		cons.(nomc) = valn;
	end
        % maximum tau_E = 30 s
	cons.temps(end) = min( cons.temps(end - 1) + 30, ...
                      max(cons.temps(end - 1) + last_zs.wth(end) ./last_zs.pin(end), 2 .* cons.temps(end - 1) - cons.temps(end-2)));

       %disp([cons.temps(end - 1) + last_zs.wth(end) ./last_zs.pin(end),2 .* cons.temps(end - 1) - cons.temps(end-2)]);
	noms = fieldnames(geo);
	for l=1:length(noms)
   		nomc = noms{l};
   		val  = geo.(nomc);
   		if ~isempty(val)
      			valn  = val(indsel);
 			valn(end+1) = val(end);
     			geo.(nomc) = valn;
   		end
	end

	noms = fieldnames(exp0d);
	tref = exp0d.temps;
	for l=1:length(noms)
   		nomc = noms{l};
   		val  = exp0d.(nomc);
  		if ~isempty(val) && (size(val,1) > 1)
			tloc = tref;
			indok = find(isfinite(tloc) & all(isfinite(val),2));
			val = val(indok,:);
			tloc = tloc(indok);
   			if length(tloc) >= 2;
      				if (size(val,1) == length(times(indok)))
					warning off
	         			valn  = interp1(tloc,val,cons.temps,'linear');		
					warning on
      	      				indbad      = find(any(~isfinite(valn),2));
     	      				if ~isempty(indbad)
	             				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
              				end
      	      				exp0d.(nomc) = valn;
      				end
			else
				exp0d.(nomc) = NaN .* ones(size(cons.temps,1),size(exp0d.(nomc),2));
   			end
		end
	end
	if isfield(exp0d,'XDURt')
		exp0d.XDURt = cons.temps;
	end

        % pour le mode evolution
	noms = fieldnames(last_zs);
	tref = last_zs.temps;
	for l=1:length(noms)
		nomc = noms{l};
		val  = last_zs.(nomc);
		if length(val) > 1
			%valn  = interp1(tref,val,cons.temps,'nearest');	
			valn  = val((end-nbte+2):end);	
			valn(end + 1)     = val(end);
			valn(end + 1)     = val(end);
			last_zs.(nomc) = valn;
		end
	end
	last_zs.temps = cons.temps;
        if option.nb_nbi == 2
            last_zs.frnbi(end-1:end) = 1 + sqrt(-1);
        else
            last_zs.frnbi(end-1:end) = 1;
        end
    
	noms = fieldnames(last_profil);
	tref = last_profil.temps;
	for l=1:length(noms)
		nomc = noms{l};
		val  = last_profil.(nomc);
		if size(val,1) > 1
			%valn = interp1(tref,val,cons.temps,'nearest');	
			valn  = val((end-nbte+2):end,:);	
			valn(end + 1,:) = val(end,:);
			valn(end + 1,:) = val(end,:);
			last_profil.(nomc) = valn;
		end
	end
	last_profil.temps  = cons.temps;		
	% appel du 0d
	warning on
	%[zs,info,profli]= zerod(option,cons,geo,exp0d,1e-2,thyb,last_zs,last_profil);
	% turn off tswitch option
        %%option.short_dt = 'off';
        option.tswitch = Inf;
        if tol0d > 0
 		[zs,info,profli]= zerod(option,cons,geo,exp0d,tol0d,thyb,last_zs,last_profil);
        else
		[zs,info,profli]= zerod(option,cons,geo,exp0d,1e-3,thyb,last_zs,last_profil);
        end
%     zs.pnbi'
%     zs.frnbi'
%     cons.pnbi'
else

	% temps pertinents des consignes
	k       = 31;
	lim     = 51;
	times   = cons.temps;
	t0     = zextraitnoeudlim(times,cons.ip,k,2,lim);
	t1     = zextraitnoeudlim(times,real(cons.iso) + imag(cons.iso),k,2,lim);
	t2     = zextraitnoeudlim(times,abs(cons.nbar),k,2,lim);
	t3     = zextraitnoeudlim(times,cons.picrh,k,2,lim);
	t4     = zextraitnoeudlim(times,cons.plh,k,2,lim);
	t5     = zextraitnoeudlim(times,real(cons.pnbi) + imag(cons.pnbi) + sqrt(real(cons.pnbi) + imag(cons.pnbi)),k,2,lim);
	t6     = zextraitnoeudlim(times,cons.pecrh,k,2,lim);
	t7     = zextraitnoeudlim(times,cons.hmore,k,2,lim);
	t8     = zextraitnoeudlim(times,real(cons.ftnbi) + imag(cons.ftnbi) + angle(cons.ftnbi),k,2,lim);
	t13     = zextraitnoeudlim(times,cons.xece,k,2,lim);
	
	% temps pertinents de la geometrie
	t9     = zextraitnoeudlim(times,geo.a,k,2,lim);
	t10    = zextraitnoeudlim(times,geo.R,k,2,lim);
	t11    = zextraitnoeudlim(times,geo.K,k,2,lim);
	t12    = zextraitnoeudlim(times,geo.d,k,2,lim);
	
	% vecteur temps pour cronos
	temps = union(t0,t1);
	temps = union(temps,t2);
	temps = union(temps,t3);
	temps = union(temps,t4);
	temps = union(temps,t5);
	temps = union(temps,t6);
	temps = union(temps,t7);
	temps = union(temps,t8);
	temps = union(temps,t9);
	temps = union(temps,t10);
	temps = union(temps,t11);
	temps = union(temps,t13);

	if isfield(option,'tswitch')
	      if isfinite(option.tswitch) && (option.tswitch >= min(cons.temps)) && (option.tswitch <= max(cons.temps))
                  ind   = find(times >= option.tswitch);
                  indsw = ind;
                  if ind < length(times)
		      indsw = cat(1,indsw,ind + 1);
                  end
                  if ind > 1
		      indsw = cat(1,indsw,ind - 1);
                  end
		  temps = union(temps,times(indsw));
	      end
	end

	if length(temps) < (length(times) /10)
		    tpl    = linspace(min(times),max(times),31)';
		    temps = union(temps,tpl);
 	end
	t14    = zextraitnoeudlim(times,cons.zeff,k,2,lim);
	temps = union(temps,t14);
	
	if ~isempty(thyb)
		indhyb = find((times >= min(thyb)) & (times <= max(thyb)));
		pas    = fix(length(indhyb)/31);
		if pas >= 1
			indhyb = indhyb(1:pas:end);
		end
		temps = union(temps,times(indhyb));
	end
	
	if length(temps) < 11
 		tpl    = linspace(min(times),max(times),11)';
		temps = union(temps,tpl);	
	end
	
	% securite
	indbad = find(diff(temps) >= ((max(temps) - min(temps)) / 10));
	nbl     = 1000;
	while(~isempty(indbad)) & (nbl > 0)
		temps =union(temps,temps(indbad) + ((max(temps) - min(temps)) / 10));
		indbad = find(diff(temps) >= ((max(temps) - min(temps)) / 10));
		nbl = nbl - 1;
	end
	
	% securite
	indbad = find(diff(temps) <= 1e-6);
	nbl     = 1000;
	while(~isempty(indbad)) & (nbl > 0)
		temps(indbad + 1) =[];
		indbad = find(diff(temps) <= 1e-6);
		nbl = nbl - 1;
	end

        % METIS fonctionne desormais pour tout pas de temps
        % on ne garde que les temps existants dans times pour eviter les artefacts
        indice_temps = unique(interp1(times,1:length(times),temps,'nearest','extrap'));
        temps = times(indice_temps);
        if length(temps) < 11
            temps = times;
        end
	
	% indication de reduction
	fprintf('Metis sample in fast mode : %d -> %d -> %d\n',length(times),length(temps),length(times));

	% modification des donnees
	noms = fieldnames(cons);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(cons,nomc);
		valn  = interp1(times,val,temps,'linear');
		indbad      = find(any(~isfinite(valn),2));
		if ~isempty(indbad)
			valn(indbad,:) = ones(length(indbad),1) * val(end,:);
		end
		cons = setfield(cons,nomc,valn);
	end
	cons.temps = temps;
	
	noms = fieldnames(geo);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(geo,nomc);
		if ~isempty(val)
			valn  = interp1(times,val,temps,'linear');
			indbad      = find(any(~isfinite(valn),2));
			if ~isempty(indbad)
				valn(indbad,:) = ones(length(indbad),1) * val(end,:);
			end
			geo = setfield(geo,nomc,valn);
		end
	end
	
	noms = fieldnames(exp0d);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(exp0d,nomc);
		if ~isempty(val)
			if (size(val,1) == length(times))
				valn  = interp1(times,val,temps,'linear');
				indbad      = find(any(~isfinite(valn),2));
				if ~isempty(indbad)
					valn(indbad,:) = ones(length(indbad),1) * val(end,:);
				end
				exp0d = setfield(exp0d,nomc,valn);
			end
		end
    end
    
	% appel du 0d
    option.short_dt = 'off';
	warning on
        if tol0d > 0
 		[zs,info,profli]= zerod(option,cons,geo,exp0d,tol0d,thyb);
        else
		[zs,info,profli]= zerod(option,cons,geo,exp0d,1e-2,thyb);
        end
	warning off
    
 	% modification des donnees de sortie
	noms = fieldnames(zs);
	for l=1:length(noms)
		nomc = noms{l};
		val  = getfield(zs,nomc);
		if length(val) == length(temps)
			val  = interp1(temps,val,times,'linear');
			zs = setfield(zs,nomc,val);
		end
	end
	profli.temps = temps;
	
	warning on


end

% handling external equilibrium
if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
    if ~isempty(equi_ext_full)
        setappdata(0,'METIS_EXTERNAL_CURDIF_EQUI',equi_ext_full);
    end
end
    




function  [xx,yy]= zextraitnoeudlim(x,y,nb,mode,lim)

if nargin < 5
 lim = inf;
end

[xx,yy]= zextraitnoeud(x,y,nb,0);
if mode == 1 
	return
elseif ~isfinite(lim)
	return
else
   while (length(xx) > lim) & (nb > 1)
	    if nb == 2
		 	nb = 1;
		 else
	    	nb = fix(nb / 2) + 1;
		 end
		 [xx,yy]= zextraitnoeud(x,y,nb,0);
   end

end
