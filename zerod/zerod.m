% calcul 0d a partir des donnees de cronos
% option = parametre du 0d
% temps  = data.gene.temps
% cons   = data.cons
% geo    = data.geo
% Vp     = volume plasma  (optionnel)
% Sp     = section poloidale du plasma  (optionnel)
% Sext   = surface externe du plasma  (optionnel)
function [zs,info,profli] = zerod(option,cons,geo,exp0d,tol0d,thyb,last_zs,last_profil)

% fonction autodeclarante
if nargin <= 1
    zs = zerod_param;
    return
end

% debut du calcul
clock_mem = clock;
cpu_mem   = cputime;

% initialisation de la structure profil
profli.xli = linspace(0,1,21);

if nargin < 6
    thyb = [];
end

if nargin < 5
    tol0d = 1e-3;
elseif isempty(tol0d)
    tol0d = 1e-3;
end
if isfield(option,'tol0d')
    if option.tol0d > 0
        tol0d = option.tol0d;
    end
end
if nargin < 4
    exp0d = [];
end

% securite
if ~isfield(exp0d,'vloop')
    exp0d.vloop = inf .* cons.temps;
end


% creation de la structure info
info = zero1t;
% compatibilite des parametres entre les differentes versions
[option,cons] = zerod_default_option(option,cons);

%  % securite logique
%  if option.matthews == 0
%  	% le rayonnement de coeur est calcule de maniere correcte et doit etre entierement deduit
%  	if option.fprad ~= 1
%  		disp('METIS info: key "fprad" must be set to 1 to be consistent with "matthews" key value 0');
%  	end
%  end
% securite option zeff
%  if (option.matthews == 1) & (option.zeff == 6)
%  	disp('METIS configuration: option.zeff set to 0 because Matthews law is already used to compute radiative power');
%  	option.zeff = 0;
%  end
if option.te_max > 1e5
    warning(sprintf('Maximum allowed temperature as been set to a larger value than 100 keV: check if others METIS setting are compatible with the new limit of %g keV',option.te_max/1e3))
end

% securite sur les impuretes
option.zimp = fix(option.zimp);
option.zmax = fix(option.zmax);
option.rimp = max(eps,min(1,option.rimp));

% external current diffusion and/or equilibrium
if isappdata(0,'EQUILIBRIUM_EXP') || isappdata(0,'CURDIFF_EXP')
    % generation on data from external data for METIS
    equi_ext = prepare_equi_ext(option,cons,geo,exp0d,profli);
    setappdata(0,'METIS_EXTERNAL_CURDIF_EQUI',equi_ext);
    % handling LCFS
    if isfield(equi_ext,'Rsepa') && isfield(equi_ext,'Zsepa')
        exp0d.Rsepa = equi_ext.Rsepa;
        exp0d.Zsepa = equi_ext.Zsepa;
        [profli.Rsepa,profli.Zsepa] = protect_sepa_z0(exp0d.Rsepa,exp0d.Zsepa,option.protect_sepa_z0);
        %profli.Rsepa = double(exp0d.Rsepa);
        %profli.Zsepa = double(exp0d.Zsepa);
    end    
    geo  = equi_ext.geo;
    cons = equi_ext.cons;
else
    % clear global data
    if isappdata(0,'METIS_EXTERNAL_CURDIF_EQUI')
        rmappdata(0,'METIS_EXTERNAL_CURDIF_EQUI');
    end
end


if nargin < 8
    
    % activation de la separatrice si disponible
    if isfield(exp0d,'Rsepa') & isfield(exp0d,'Zsepa')
        if all(isfinite(exp0d.Rsepa(:))) && all(isfinite(exp0d.Zsepa(:)))
            [profli.Rsepa,profli.Zsepa] = protect_sepa_z0(exp0d.Rsepa,exp0d.Zsepa,option.protect_sepa_z0);
            %profli.Rsepa = double(exp0d.Rsepa);
            %profli.Zsepa = double(exp0d.Zsepa);
        else
            fprintf('S');
        end
    end
    
    % activation XDUR si disponible
    profli.fxdurplh  = NaN .* ones(size(cons.temps)) * profli.xli;
    if isfield(exp0d,'XDURt') & isfield(exp0d,'XDURx') & isfield(exp0d,'XDURv')
        indok            = find(all(exp0d.XDURv>=0,2));
        if length(indok) >3
            fplh             = interp1(exp0d.XDURt(indok),exp0d.XDURv(indok,:),cons.temps,'nearest','extrap');
            vt = ones(size(cons.temps));
            if size(exp0d.XDURx,1) > 1
                xfplh            = min(1,max(0,interp1(exp0d.XDURt(indok),exp0d.XDURx(indok,:),cons.temps,'nearest','extrap')));
                profli.fxdurplh  = tsplinet(xfplh,fplh,vt * profli.xli);
            else
                profli.fxdurplh  = tsplinet(vt*exp0d.XDURx,fplh,vt * profli.xli);
            end
        end
    end
    
    
    % selon la precision demandee
    if (tol0d <= 1e-3)
        fprintf('initial run of ');
        %option.notimederivative = 1;
        option_loc = option;
        option_loc.tol0d = 2e-3;
        [zs,info,profli]= zerodfast(option_loc,cons,geo,exp0d,thyb);
        %option.notimederivative = 0;
        % reechantillonage des profils
        noms = fieldnames(profli);
        temps_pout = profli.temps;
        for k = 1:length(noms)
            nomc = noms{k};
            if size(profli.(nomc),1) == length(temps_pout)
                if all(~isfinite(profli.(nomc)(:)))
                    profli.(nomc) = NaN .* ones(size(cons.temps,1),size(profli.(nomc),2));
                elseif any(~isfinite(profli.(nomc)(:)))
                    profli.(nomc) = interp1(temps_pout,profli.(nomc),cons.temps,'nearest','extrap');
                else
                    profli.(nomc) = interp1(temps_pout,profli.(nomc),cons.temps,'linear','extrap');
                end
            end
        end
        
        % preserve les donnees de l'interpolation
        % activation de la separatrice si disponible
        if isfield(exp0d,'Rsepa') & isfield(exp0d,'Zsepa')
            if all(isfinite(exp0d.Rsepa(:))) && all(isfinite(exp0d.Zsepa(:)))
                [profli.Rsepa,profli.Zsepa] = protect_sepa_z0(exp0d.Rsepa,exp0d.Zsepa,option.protect_sepa_z0);
                %profli.Rsepa = double(exp0d.Rsepa);
                %profli.Zsepa = double(exp0d.Zsepa);
            else
                fprintf('S');
            end
        end
        
        % activation XDUR si disponible
        profli.fxdurplh  = NaN .* ones(size(cons.temps)) * profli.xli;
        if isfield(exp0d,'XDURt') & isfield(exp0d,'XDURx') & isfield(exp0d,'XDURv')
            indok            = find(all(exp0d.XDURv>=0,2));
            if length(indok) >3
                fplh             = interp1(exp0d.XDURt(indok),exp0d.XDURv(indok,:),cons.temps,'nearest','extrap');
                vt = ones(size(cons.temps));
		if size(exp0d.XDURx,1) > 1
		    xfplh            = min(1,max(0,interp1(exp0d.XDURt(indok),exp0d.XDURx(indok,:),cons.temps,'nearest','extrap')));
		    profli.fxdurplh  = tsplinet(xfplh,fplh,vt * profli.xli);
		else
		    profli.fxdurplh  = tsplinet(vt*exp0d.XDURx,fplh,vt * profli.xli);
		end
            end
        end

        
        
        option.init_geo = 1;
        option_eff = option;
        
        % initialisation de la bocucle de convergence
        nb         = 4;
        mem        = zs;
        proflimem  = profli;
        f          = 0.3;
        amorti     = 0.5;
        nb_plus3   = length(temps_pout);
        
        fprintf('start full run of METIS now ->\n');
        [hdlg,value] = z0dpatience;
        if ~isempty(value)
            set(hdlg,'userdata',[1-value value]);
        else
            set(hdlg,'userdata',[0.75 0.25]);
        end
    else
        % premier appel pour creation du jeux de donnees
        option_eff = option;
        option_eff.transitoire = 0;
        zs = zero1t(option_eff,cons,geo,[],0.9,profli);
        
        % initialisation de la bocucle de convergence
        nb      = 0;
        zs.dw      = inf;
        zs.dpfus   = inf;
        zs.dini    = inf;
        zs.diboot  = inf;
        zs.temps   = cons.temps;
        mem        = zs;
        if nargin >=8
            proflimem  = profli;
        else
            proflimem  = [];
        end
        f          = 1;
        amorti     = 0.5;
        nb_plus3   = length(cons.temps);
        
    end
    
else
    
    % premier appel pour creation du jeux de donnees
    option_eff = option;
    option_eff.transitoire = 1;
    zs  	 = last_zs;
    geo.vp   = zs.vp;
    geo.sp   = zs.sp;
    geo.sext = zs.sext;
    profli   = last_profil;
    
    
    % activation de la separatrice si disponible
    if isfield(exp0d,'Rsepa') & isfield(exp0d,'Zsepa')
        if all(isfinite(exp0d.Rsepa(:))) && all(isfinite(exp0d.Zsepa(:)))
            [profli.Rsepa,profli.Zsepa] = protect_sepa_z0(exp0d.Rsepa,exp0d.Zsepa,option.protect_sepa_z0);
            %profli.Rsepa = double(exp0d.Rsepa);
            %profli.Zsepa = double(exp0d.Zsepa);
        else
            fprintf('S');
        end
    elseif isfield(profli,'Rsepa') & isfield(profli,'Zsepa')
        if ~all(isfinite(profli.Rsepa(:))) || ~all(isfinite(profli.Zsepa(:)))
            fprintf('S');
            profli = rmfield(rmfield(profli,'Rsepa'),'Zsepa');
        end
    end
    
    % activation XDUR si disponible
    profli.fxdurplh  = NaN .* ones(size(cons.temps)) * profli.xli;
    if isfield(exp0d,'XDURt') & isfield(exp0d,'XDURx') & isfield(exp0d,'XDURv')
        indok            = find(all(exp0d.XDURv>=0,2));
        if length(indok) >3
            fplh             = interp1(exp0d.XDURt(indok),exp0d.XDURv(indok,:),cons.temps,'nearest','extrap');
            vt = ones(size(cons.temps));
            if size(exp0d.XDURx,1) > 1
                xfplh            = min(1,max(0,interp1(exp0d.XDURt(indok),exp0d.XDURx(indok,:),cons.temps,'nearest','extrap')));
                profli.fxdurplh  = tsplinet(xfplh,fplh,vt * profli.xli);
            else
                profli.fxdurplh  = tsplinet(vt*exp0d.XDURx,fplh,vt * profli.xli);
            end
         end
    end
    
    % initialisation de la bocucle de convergence
    nb      = 0;
    zs.dw      = inf;
    zs.dpfus   = inf;
    zs.dini    = inf;
    zs.diboot  = inf;
    zs.temps   = cons.temps;
    mem        = zs;
    if nargin >=8
        proflimem  = profli;
    else
        proflimem  = [];
    end
    f          = 1;
    amorti     = 0.5;
    nb_plus3   = length(cons.temps);
    
end


% securite sur iso
cons.iso = min(30,max(0,real(cons.iso))) + sqrt(-1) .* min(30,max(0,imag(cons.iso)));

% securite sur ip
if option.berror ~= 0
    cons.ip(2:end) = min(1e12,max(1,cons.ip(2:end)));
else
    cons.ip = min(1e12,max(1,cons.ip));
end
% securite sur nbar
cons.nbar = min(1e25,max(1e13,cons.nbar));
% securite sur a
geo.a = min(1e3,max(0.01,geo.a));
% securite sur R
geo.R = min(3e3,max(geo.a,geo.R));
% securite sur K
geo.K = min(5,max(0.1,geo.K));


option_eff.matthews    = 1;
if option.lhmode  == 1
    option_eff.lhmode = 2;
end


% si plotconv
if option.plotconv == 1
    option
end

if option.scaling == 4
    if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') && (option.exp_shape == 0)
	fprintf('METIS: energy fiting mode can''t be used in same time than Te and Ti are provided as external data; switching to standard mode.\n');
	option.scaling = 0;
    else
	fprintf('METIS: switching to energy fiting mode (using default scaling as starting point without ITB, MHD and dilution effect.\n');
    end
end
if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') && any(cons.hmore ~= 0) && (option.exp_shape == 0)
	fprintf('METIS: energy enhancement factor can''t be used in same time than Te and Ti are provided as external data; setting back  hmore reference to 1.\n');
        cons.hmore(:) = 1;
end

% pas de temps court 
nbptal = 3;
if option_eff.evolution == 1
    tinter = cons.temps;
else
    tinter = [];
end
nbt_loc  = length(cons.temps);
switch option.short_dt
  case 'full'
    noms_loc = fieldnames(cons);
    for k = 1:length(noms_loc)
      switch noms_loc{k}
        case 'temps'
          %nothing
        otherwise
            if length(cons.(noms_loc{k})) == nbt_loc
                cons.(noms_loc{k}) = z0sgolayfilt(cons.(noms_loc{k}),1,nbptal,[],1,tinter);
            end
        end
    end
    noms_loc = fieldnames(geo);
    for k = 1:length(noms_loc)
     switch noms_loc{k}
        case 'temps'
          %nothing
        otherwise
            if length(geo.(noms_loc{k})) == nbt_loc
                geo.(noms_loc{k}) = z0sgolayfilt(geo.(noms_loc{k}),1,nbptal,[],1,tinter);
            end
        end    
    end
    
    if isfield(profli,'Rsepa') && isfield(profli,'Zsepa') && ...
       ~isempty(profli.Rsepa)  && ~isempty(profli.Zsepa)
           profli.Rsepa = z0sgolayfilt( profli.Rsepa,1,nbptal,[],1,tinter);
           profli.Zsepa = z0sgolayfilt( profli.Zsepa,1,nbptal,[],1,tinter);       
    end
    
end
 


% pour le fit de vloop avec LH
etalhfit = 5e18;
etalhmem = 5e18;
feta  = 1e19;
scale_etalh = 1e19;
erreta = inf;
erretalast = inf;
errposmem = inf;
erretamem = inf;
change = 1;
count_stable = 0;
postab = [];
etatab = [];
fetatab = [];
plhmem = cons.plh;
nbarref = cons.nbar;


% flag d'initialisation de nep
if (option.ane >= 5) && (option.ane < 10)
    fnepini = 1;
else
    fnepini = 0;
end

% pre convergence for BgB and CDBM model if activated in standart mode
% during 3 first calls.
if option.kishape  == 0
   option_eff.kishape = 3;
end

% fit de wdia
switch option.gaz
    case {3,5,11}
        % with fusion power that can be dominating, more convergence loop
        % may be needed
        option.nbmax = option.nbmax + nb_plus3;
end
hmorefit   = cons.hmore;
dhmorefit  = 0;
modeh_nonconv = 1;
cristalisation = 0;
% boucle
nbtotbar = option.nbmax + double((option.vloop == 3) ||(option.vloop == 6))  .* 21 + (option.lhmode == 1) .* 31 + (option.scaling == 4) .* 11 + 1;
setappdata(0,'STOP_METIS',0);
fprintf('@:');
while  ((nb < (option.nbmax + double((option.vloop == 3) ||(option.vloop == 6)) .* 21 + (option.lhmode == 1) .* 31 + (option.scaling == 4) .* 11)) &  ...
        ((mean(zs.dw) > tol0d) | (mean(zs.dpfus) > tol0d) | ...
        (mean(zs.dini) > tol0d)| (mean(zs.diboot) > tol0d) | (nb < (7 + (option.lhmode == 1) .* 8)) | (change == 1) | (dhmorefit > tol0d) | modeh_nonconv))
    % compteur
    nb    = nb +1;
    %fprintf('.');
    change = 0;
    symbol ='.';
    
    % appel de zero1t
    if nb > 3
        option_eff = option;
    end
    if (option.lhmode  == 1) & (nb <= 10)
        change = 1;
        count_stable = 7;
    elseif (option.lhmode  == 1) & (nb > 10)
        option_eff.etalh  = etalhfit;
        if (abs(feta ./ scale_etalh) >= tol0d)
            change = 1;
        end
    end
    if (change == 0) && (option.lhmode  == 1)
        count_stable =  count_stable -1;
        if count_stable > 0
            change = 1;
        end
        %disp([etalhfit,mean(zs.etalh0)]);
    end
    
    % asservissement a tension par tour nulle avec LHCD
    if (option_eff.vloop == 3) & (nb > 3)
        dip        = (cons.ip - zs.ip);
        %dip        = f .* dip;
        %dip        = (dip + cat(1,dip(2:end),dip(end))) ./ 2;
        dip        = cat(1,dip(2:end),dip(end));
        cons.plh = max(plhmem ./ 10, (2 .* cons.plh + dip  ./ zs.etalh0 .* zs.nbar .* geo.R .* ...
            (plhmem > (max(plhmem)/100))) ./ 2);
        cons.plh = min(1.5 .* max(plhmem),cons.plh);
        %        figure(20);
        %        subplot(3,1,1)
        %        plot(cons.temps,cons.plh,cons.temps,plhmem);
        %        subplot(3,1,2)
        %        plot(cons.temps,cons.ip,cons.temps,zs.ip,cons.temps,mem.ip);
        %        subplot(3,1,3)
        %        plot(cons.temps,option.vref .* ones(size(cons.temps)),cons.temps,zs.vloop,cons.temps,mem.vloop);
        %        drawnow
    elseif (option_eff.vloop == 6) & (nb > 3)
        % asservissement a tension par tour nulle avec NBICD
        dip        = (cons.ip - zs.ip);
        dip        = cat(1,dip(2:end),dip(end));
        eta_nbi    = real(zs.inbicd) ./ max(1,real(zs.pnbi)) + sqrt(-1) .* imag(zs.inbicd) ./ max(1,imag(zs.pnbi));
        cons.pnbi  = cons.pnbi + 0.7 .* eta_nbi .* dip;
        %        figure(20);
        %        subplot(3,1,1)
        %        plot(cons.temps,real(cons.pnbi),'r',cons.temps,imag(cons.pnbi),'b',cons.temps,real(zs.pnbi),'m',cons.temps,imag(zs.pnbi),'c');
        %        subplot(3,1,2)
        %        plot(cons.temps,cons.ip,cons.temps,zs.ip,cons.temps,mem.ip,cons.temps,dip .* 10);
        %        subplot(3,1,3)
        %        plot(cons.temps,option.vref .* ones(size(cons.temps)),cons.temps,zs.vloop,cons.temps,mem.vloop);
        %        drawnow
    end
    
    % asservissement sur la fraction rayonnees
    if (option_eff.neasser == 2) & (nb > 3)
        epn  = (zs.ploss ./ zs.pin) - 0.2;   % seuil d'action a 80 % de puissance rayonnees
        % on declenche a 50 % de puissance rayonnee
        cons.nbar = min(nbarref,max(nbarref ./ 20,cons.nbar .* (1 + epn .* (epn < 0.3))));
    end
    
    if (nb < 3) & (option.evolution == 0)
        option_eff.frad = option.frad  .* nb./ 3;
    else
        option_eff.frad = option.frad;
    end
    
    if option.scaling == 4
        option_eff.scaling = 0;
        option_eff.dilution = 0;
        option_eff.sitb = 0;
        option_eff.smhd = 100;
        cons.hmore = hmorefit;
    end
    
    switch option.dwdt_method
        case 'working_point'
            option_eff.dwdt_method = 'none';
    end
    %   profli_in = profli;
    last = zs;
    [zs,profli] = zero1t(option_eff,cons,geo,mem,amorti,profli);
    
    %   save(sprintf('metis_loc_%d',nb),'zs','profli','option_eff','cons','geo','mem','amorti','profli_in');
    
    if option_eff.evolution == 1
        drawnow
    end
    
    % fin init_geo forcee
    option_eff.init_geo = 0;
    option.init_geo     = 0;
    
    
    % calcul de la variation
    %     zs.dw      = znonan(abs(zs.wth - mem.wth) ./ abs(zs.wth + mem.wth+1));
    %     zs.dpfus   = znonan(abs(zs.pin - mem.pin) ./  abs(zs.pin + mem.pin+1));
    %     zs.dini    = znonan(abs(zs.ini - mem.ini) ./  abs(zs.ini + mem.ini+1));
    %     zs.diboot  = znonan(abs(zs.iboot - mem.iboot) ./  abs(zs.iboot + mem.iboot+1));
    
    zs.dw      = znonan(abs(zs.wth - last.wth) ./ abs(mem.wth+1));
    zs.dpfus   = znonan(abs(zs.pin - last.pin) ./  abs(mem.pin+1));
    zs.dini    = znonan(abs(zs.ini - last.ini) ./  abs(mem.ini+1));
    zs.diboot  = znonan(abs(zs.iboot - last.iboot) ./  abs(mem.iboot+1));
    
    
    % convergence mode h
    if all(zs.modeh == last.modeh)
        modeh_nonconv	= 0;
    elseif cristalisation == 0
        modeh_nonconv	= 1;
    else
        modeh_nonconv	= 0;
    end
    
    %  figure(21)
    %  subplot(4,1,1)
    %  semilogy(zs.temps,zs.dw,'r',mem.temps,mem.dw,'b')
    %  set(gca,'ylim',[1e-3 1]);
    %  hold on
    %  subplot(4,1,2)
    %  semilogy(zs.temps,zs.dpfus,'r',mem.temps,mem.dpfus,'b')
    %  set(gca,'ylim',[1e-3 1]);
    %  hold on
    %  subplot(4,1,3)
    %  semilogy(zs.temps,zs.dini,'r',mem.temps,mem.dini,'b')
    %  set(gca,'ylim',[1e-3 1]);
    %  hold on
    %  subplot(4,1,4)
    %  semilogy(zs.temps,zs.diboot,'r',mem.temps,mem.diboot,'b')
    %  set(gca,'ylim',[1e-3 1]);
    %  hold on
    
    
    if option.evolution == 1
        zs.dw(:)      = zs.dw(end-1);
        zs.dpfus(:)   = zs.dpfus(end-1);
        zs.dini(:)    = zs.dini(end-1);
        zs.diboot(:)  = zs.diboot(end-1);
    end
    
    if option.plotconv == 1
        % figure convergence
        h = findobj(0,'type','figure','tag','err0d');
        if isempty(h)
            h = figure('tag','err0d');
        else
            figure(h)
        end
        clf
        semilogy(cons.temps,zs.dw,cons.temps,zs.dpfus,cons.temps,zs.dini, cons.temps, ...
            zs.diboot);
        title(sprintf('%g, #%d',f,nb))
        set(gca,'ylim',[1e-4 inf]);
        drawnow
    end
    
    %figure(19);clf;subplot(2,1,1);plot(cons.temps,cons.li,'b',cons.temps,zs.li,'r');subplot(2,1,2);plot(cons.temps,exp0d.vloop,'b',cons.temps,zs.vloop,'r');drawnow;
    
    % amortissement
    if (nb > fix(option.nbmax ./ 1.45)) && (change == 0)
        f = 0.3 .* f;
        symbol ='~';
        cristalisation = 1;
    elseif nb ==1
        f = 1;
        symbol ='i';
    elseif (nb < fix(option.nbmax ./ 2.5))
        f = 0.3;
        symbol ='.';
    elseif (option.lhmode  ~= 1)
        f = 0.9 .* f;
        symbol ='.';
        cristalisation = 1;
    end
    amorti = f .* (5/4);
    
    %     figure(187)
    %     subplot(3,1,1)
    %     plot(cons.temps,zs.wth,cons.temps,mem.wth);
    %     hold
    %     subplot(3,1,2)
    %     plot(cons.temps,zs.dwthdt,'b',cons.temps,mem.dwthdt,'r', ...
    %          cons.temps,z0dxdt(zs.wth,cons.temps),'c',cons.temps,z0dxdt(mem.wth,cons.temps),'m');
    %     hold
    %     subplot(3,1,3)
    %     plot(cons.temps,zs.ploss,'b',cons.temps,mem.ploss,'r')
    %     hold
    %     drawnow
    
    %    figure(21);clf
    %    subplot(3,1,1)
    %    semilogy(cons.temps,zs.telim,cons.temps,last.telim,cons.temps,mem.telim);
    %    subplot(3,1,2)
    %    semilogy(cons.temps,zs.nelim,cons.temps,last.nelim,cons.temps,mem.nelim);
    %    subplot(3,1,3)
    %    semilogy(cons.temps,zs.nelim .* zs.telim ./ zs.nebord ./ zs.tebord, ...
    %         cons.temps,last.nelim .* last.telim ./ last.nebord ./ last.tebord, ...
    %         cons.temps,mem.nelim .* mem.telim ./ mem.nebord ./ mem.tebord);
    %    drawnow
    %    %keyboard
    
    if option.evolution == 0
        mem            = zamor0(zs,mem,f);
        switch option.dwdt_method
            case 'working_point'
                mem = working_point(mem);
        end
        mem.modeh      = zs.modeh;
        mem.xpoint     = zs.xpoint;
        mem.asser      = zs.asser;
        mem.difcurconv = zs.difcurconv;
        mem.disrup     = zs.disrup;
        mem.disrup     = zs.disrup;
        
        mem.dw         = zs.dw;
        mem.dini       = zs.dini;
        mem.dpfus      = zs.dpfus;
        mem.diboot     = zs.diboot;
        
        switch option.sol_model
            case '2_points'
                %mem.telim = sqrt(zs.telim .* mem.telim);
                %mem.nelim = sqrt(zs.nelim .* mem.nelim);
                mem.telim = zs.xpoint .* zs.telim +  (~zs.xpoint) .* mem.telim;
                mem.nelim = zs.xpoint .* zs.nelim +  (~zs.xpoint) .* mem.nelim;
        end
        
    elseif nb == 1
        % pour le mode evolution
        noms = fieldnames(mem);
        for l=1:length(noms)
            nomc = noms{l};
            if length(mem.(nomc)) > 1
                mem.(nomc)(end-1:end) = zs.(nomc)(end-1:end);
            else
                mem.(nomc) = zs.(nomc);
            end
        end
    else
        memn            = zamor0(zs,mem,f);
        memn.modeh      = zs.modeh;
        memn.asser      = zs.asser;
        memn.difcurconv = zs.difcurconv;
        memn.disrup     = zs.disrup;
        memn.dw         = zs.dw;
        memn.dini       = zs.dini;
        memn.dpfus      = zs.dpfus;
        memn.diboot     = zs.diboot;
        
        switch option.sol_model
            case '2_points'
                %memn.telim = sqrt(zs.telim .* memn.telim);
                %memn.nelim = sqrt(zs.nelim .* memn.nelim);
                memn.telim = zs.xpoint .* zs.telim +  (~zs.xpoint) .* memn.telim;
                memn.nelim = zs.xpoint .* zs.nelim +  (~zs.xpoint) .* memn.nelim;
        end
        
        % pour le mode evolution
        noms = fieldnames(memn);
        for l=1:length(noms)
            nomc = noms{l};
            if length(mem.(nomc)) > 1
                mem.(nomc)(end-1:end) = memn.(nomc)(end-1:end);
            else
                mem.(nomc) = memn.(nomc);
            end
        end
        
        
    end
    
    % les profils
    if isempty(proflimem)
        proflimem = profli;
        
    elseif option.evolution == 0
        
        proflimem          = zamor0(profli,proflimem,f);
        % cas du profil de densite
        if all(mem.n0a >1) & (fnepini  ==1)
            proflimem.nep = profli.nep;
            fnepini = 0;
            
        end
        % cas du mode fit de l'efficacite LH
        if option.lhmode  == 1
            proflimem.jlh    = proflimem.jlh  .* ((mem.ilh./  max(eps,trapz(proflimem.xli,proflimem.spr .* proflimem.jlh,2))) * ones(size(proflimem.xli)));
        end
        switch option.dwdt_method
            case 'working_point'
                proflimem = working_point(proflimem);
        end
        
        % recopie
        profli             = proflimem;
    elseif nb == 1
        % pour le mode evolution
        noms = fieldnames(proflimem);
        for l=1:length(noms)
            nomc = noms{l};
            if size(proflimem.(nomc),1) > 1
                proflimem.(nomc)(end-1:end,:) = profli.(nomc)(end-1:end,:);
            else
                proflimem.(nomc) = profli.(nomc);
            end
        end
        % recopie
        profli             = proflimem;
        
    else
        
        proflimemn          = zamor0(profli,proflimem,f);
        % cas du profil de densite
        if all(mem.n0a >1) & (fnepini  ==1)
            proflimemn.nep = profli.nep;
            fnepini = 0;
        end
        % pour le mode evolution
        noms = fieldnames(proflimemn);
        for l=1:length(noms)
            nomc = noms{l};
            if size(proflimem.(nomc),1) > 1
                proflimem.(nomc)(end-1:end,:) = proflimemn.(nomc)(end-1:end,:);
            else
                proflimem.(nomc) = proflimemn.(nomc);
            end
        end
        % recopie
        profli             = proflimem;
    end
    
    %    % si pas de derivee spatiale pour la preconvergence en mode rapide
    %    if option.notimederivative  == 1
    %  	zs.dwbpdt(:) = 0;
    %  	zs.dwdt(:) = 0;
    %  	zs.dwthdt(:) = 0;
    %  	zs.dwmagtordt(:) = 0;
    %          disp('ici')
    %    end
    
    % recherche de etalh et pos
    if (option.lhmode  ~= 1)
        % on ne fait rien
    elseif nb < 11
        % preambule
        erretamem  = z0delta(exp0d.vloop,zs.vloop,cons.plh,cons.temps,thyb);
        best       = mem;
    elseif (option.lhmode  == 1)
        
        % recherche de l'efficacite
        erreta  = z0delta(exp0d.vloop,zs.vloop,cons.plh,cons.temps,thyb);
        errip   = z0delta(cons.ip,zs.ilh+ zs.iboot,cons.plh,cons.temps,thyb);
        if (abs(feta ./ scale_etalh) < tol0d)
            % on arrete
            chpos = 1;
        elseif (errip < 0)
            % il faut reduire
            etalhfit = max(1e17,etalhmem - feta);
            mem      = best;
            erreta = erretamem;
            % on diminue le pas
            feta = abs(feta) ./ sqrt(2);
            symbol ='#';
        elseif (abs(erreta) >= abs(erretamem))
            % le nouvel essai est moins bon que le precedent
            % on diminue le pas
            feta = abs(feta) ./ sqrt(2);
            % le signe est change si l'erreur change de signe
            %if (sign(erreta) ~= sign(erretamem))
            %	feta = - abs(feta) ./ sqrt(2);
            %else
            %	feta = abs(feta) ./ sqrt(2);
            %end
            % retour a l'etat precedent
            etalhfit = etalhmem;
            erreta = erretamem;
            mem    = best;
            etalhfit = min(3e19,max(-3e19.* (option.etalh <=0),etalhfit  - sign(erreta)  .* feta));
        else
            % C'EST MEILLEUR
            % on conserve l'essai
            erretamem = erreta;
            etalhmem  = etalhfit;
            best      = mem;
            %if (errip < 0)
            %	etalhfit = min(3e19,max(-3e19.* (option.etalh <=0),etalhfit  - abs(feta)));
            %	feta = abs(feta) ./ sqrt(2);
            %else
            %	etalhfit = min(3e19,max(-3e19.* (option.etalh <=0),etalhfit  - sign(erreta) .* feta));
            %end
            etalhfit = min(3e19,max(-3e19.* (option.etalh <=0),etalhfit  - sign(erreta) .* feta));
            if abs(etalhfit) == 3e19
                feta = abs(feta) ./ sqrt(2);
            end
            if  sign(erreta) < 0
                symbol ='^';
            else
                symbol ='v';
            end
        end
    end
    
    etatab = cat(1,etatab,etalhmem+ sqrt(-1) .* etalhfit);
    fetatab = cat(1,fetatab,feta);
    
    if option.scaling == 4
        hmorefitmem   = hmorefit;
        switch option.machine
            case 'TS'
                hmorefit   = medfilt1(min(2,max(0.3,hmorefit .* ( 1 + 0.3 .* (max(eps,exp0d.w) - zs.wdia) ./  ...
                    (max(eps,exp0d.w) + mem.wdia)))),5);
            otherwise
                hmorefit   = medfilt1(min(2,max(0.3,hmorefit .* ( 1 + 0.3 .* (medfilt1(max(eps,exp0d.w),3) - zs.w) ./  ...
                    (max(eps,exp0d.w) + mem.w)))),5);
        end
        dhmorefit = abs(hmorefit - hmorefitmem); % la normalisation est 1
    end
    
    % pas de temps court   
    switch option.short_dt
        case 'full'
            noms_loc = fieldnames(mem);
            for k = 1:length(noms_loc)
                switch noms_loc{k}
                    case 'temps'
                        %nothing
                    otherwise
                        if length(mem.(noms_loc{k})) == nbt_loc
                            mem.(noms_loc{k}) = z0sgolayfilt(mem.(noms_loc{k}),1,nbptal,[],1,tinter);
                        end
                end
            end
            noms_loc = fieldnames(profli);
            for k = 1:length(noms_loc)
                switch noms_loc{k}
                    case {'temps','Rsepa','Zsepa'}
                        %nothing
                    otherwise
                        if size(profli.(noms_loc{k}),1) == nbt_loc
                            profli.(noms_loc{k}) = z0sgolayfilt(profli.(noms_loc{k}),1,nbptal,[],1,tinter);
                        end
                end
            end   
        case 'on'
            noms_loc = fieldnames(mem);
            for k = 1:length(noms_loc)
                switch noms_loc{k}
                    case {'pohm','vloop','vmes','drmdt','irun', ...
                          'dwmagtordt','dwbpdt'}
                        if length(mem.(noms_loc{k})) == nbt_loc
                            mem.(noms_loc{k}) = z0sgolayfilt(mem.(noms_loc{k}),1,nbptal,[],1,tinter);
                        end
                end
            end
            noms_loc = fieldnames(profli);
            for k = 1:length(noms_loc)
                switch noms_loc{k}
                    case {'jrun','pohm','dpsidt','epar','ej','tep','tip','qei'}
                        if size(profli.(noms_loc{k}),1) == nbt_loc
                            profli.(noms_loc{k}) = z0sgolayfilt(profli.(noms_loc{k}),1,nbptal,[],1,tinter);
                        end
                end
            end 
            mem.te0  = interp1(profli.temps,profli.tep(:,1),cons.temps,'linear','extrap');
            tem_loc  = trapz(profli.xli,profli.tep .* profli.vpr,2) ./ ...
                          max(eps,trapz(profli.xli,profli.vpr,2));
            tim_loc  = trapz(profli.xli,profli.tip .* profli.vpr,2) ./ ...
                          max(eps,trapz(profli.xli,profli.vpr,2));
            mem.tem  = interp1(profli.temps,tem_loc,cons.temps,'linear','extrap');
            tite_loc = max(1/10,min(10,tim_loc ./ tem_loc));
            mem.tite = interp1(profli.temps,tite_loc,cons.temps,'linear','extrap');

    end
    
    fprintf(symbol);
    
    if option.evolution ~= 0
        drawnow
    else
        z0dpatience(nb/nbtotbar);
    end
    
end

if option.evolution == 1
    if zs.disrup(end-1) > 0
        fprintf('D');
    end
end


% fin du clacul
fprintf(' in %g s (cpu = %g s)\n',etime(clock,clock_mem),cputime - cpu_mem);

% sortie coherente avec les profils
zs = mem;
%  zcompstruct(zs,mem);
%  zplotstruct(zs,mem,'o');
%  disp('====================')
%  zcompstruct(profli,proflimem);

if option.lhmode == 1;
    h = findobj(0,'type','figure','tag','z0dasser');
    if isempty(h)
        h=figure('tag','z0dasser');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',1,'color',[1 1 1])
    plot(1:length(etatab),real(etatab),'r',1:length(etatab),real(etatab),'.b', ...
        1:length(etatab),imag(etatab),'g+',1:length(fetatab),fetatab,'oc')
    title(sprintf('LH current : %g A/w/m^2\n',etalhmem));
    xlabel('# of call');
    ylabel('LH current drive efficiency (A/w/m^2)');
    
end

if option.evolution == 1
    if (zs.dw(end-1) > 0.05) | (zs.dini(end-1) > 0.05)| (zs.dpfus(end-1) > 0.05)| (zs.diboot(end-1) > 0.05)
        disp('Warning : bad convergence : check it !');
    end
elseif any(zs.dw(:) > 0.05) | any(zs.dini(:) > 0.05)| any(zs.dpfus(:) > 0.05)| any(zs.diboot(:) > 0.05)
    disp('Warning : bad convergence at some time step : check it !');
end

if option.evolution == 1
    % rien
elseif any(zs.disrup(:) > 0)
    ind = min(find(zs.disrup~=0));
    fprintf('possible shot disruption @ %g s (%d): radiative power exceeds input power; you must check the result\n',cons.temps(ind),ind);
end

%zs.nbar  = cons.nbar;
zs.nb    = nb;

% sauvegarde de hmorefit dans hitb en cas de fit sur Wdia
if option.scaling == 4
    zs.hitb = hmorefit;
end

% calcul de la rotation toroidale
%zs.wrad = z0rot(zs,option,cons,geo);
%  [zs.wrad,void_snbi,void_sicrh,void_sfus,void_sripth, ...
%           void_sriplh,void_sripicrh,void_sturb,void_fact,profli] = ...
%                                          z0rot3(zs,option,cons,geo,profli);


% si besoin plot de LH
if isfield(profli,'fxdurplh') & (option.dlh == 0) & (option.wlh == 0)
    zs.efficiency = zs.etalh0;
elseif option.wlh == 0
    zs.efficiency = zs.etalh0;
else
    % rien
end

% recherche de donnees non valides
% recherche de imag et nan
noms = fieldnames(zs);
noms(strmatch('stf',noms,'exact')) = [];
stf  = 0;
for kzw = 1:length(noms)
    nomc = noms{kzw};
    %varc = getfield(zs,nomc);
    varc = zs.(nomc);
    if any(~isfinite(varc(:)))
        fprintf('NaN in %s\n',nomc);
        stf = stf + 1;
    end
    if (isempty(findstr(nomc,'idn')) &&  isempty(findstr(nomc,'nbi'))) || (option.nb_nbi == 1)
        if any(imag(varc(:)))
            fprintf('Imag in %s\n',nomc);
            stf = stf + sqrt(-1);
            if option.evolution == 1
                zs.(nomc) = real(zs.(nomc));
            end
        end
    end
end
zs.stf     =  stf;

profli.temps = cons.temps;

% % recompute precise xii and xie after wave form relaxation algorithme
% ee   = 0.1602176462e-18;
% rhomax = profli.rmx(:,end);
% ve     = ones(size(profli.xli));
% % calcul des coefficients de transport associes
% tepd1       = pdederive(profli.xli,profli.tep,0,2,2,1);
% tipd1       = pdederive(profli.xli,profli.tip,0,2,2,1);
% % attention vpr = V' * rhomax !
% xiea         = - profli.qe  ./ max(eps, profli.vpr) ./ (ee .* tepd1) ./ ...
%                   max(1e13,profli.nep) ./ profli.grho2 .* (rhomax * ve) .^ 2;
% xiia         = - profli.qi ./ max(eps,profli.vpr) ./ (ee .* tipd1) ./  ...
%                  max(1e13,profli.nip) ./ profli.grho2.* (rhomax * ve) .^ 2;
% %  prolongation
% xiea         = pchip(profli.xli(2:end),real(xiea(:,2:end)),profli.xli);
% xiia         = pchip(profli.xli(2:end),real(xiia(:,2:end)),profli.xli);
% % limitation
% profli.xie  = xiea;
% profli.xii  = xiia;


% si besoin plot de LH
if isfield(profli,'fxdurplh') & (option.dlh == 0)
    % rien
elseif option.wlh == 0
    % rien
elseif option.evolution == 0
    vt = ones(size(cons.temps));
    switch option.gaz
        case 1
            agaz    = 1 .* vt;
            zgaz    = 1 .*vt;
        case 2
            agaz    = 2 .* vt;
            zgaz    = 1 .* vt;
        case 3
            agaz    = 2.5 .* vt;
            zgaz    = 1 .* vt;
        case 5
            % same amount of D and He3 ?
            agaz    = 2.5 .* vt;
            zgaz    = 1.5 .* vt;
        case 11
            % less than 20% of boron
            agaz    = 1 .* vt;
            zgaz    = 1 .* vt;
        otherwise
            agaz    = 4 .* vt;
            zgaz    = 2 .* vt;
    end
    % directivite
    switch option.lhmode
        case 0
            directivity = 0.8;
        case 3
            directivity = abs(option.etalh);
        case 4
            directivity = abs(option.etalh);
        otherwise
            directivity = 0.8;
    end
    % facteur de propagation LHCD
    qcyl          =  5 .* geo.a .^ 2 .* geo.b0 ./ (zs.ip./1e6) ./ geo.R.* ( 1 + (geo.a ./ geo.R) .^ 2);
    qbord = profli.qjli(:,end);
    switch option.upshiftmode
        case {'newmodel','newmodel + tail'}
            upshift      = option.fupshift .* vt;
        otherwise
            upshift      = option.fupshift .* max(eps,qbord ./ qcyl - 1);
    end

    z0lhacc2lobes(option.freqlh.*1e9,option.npar0,option.wlh,agaz,zgaz,cons.temps,profli.xli,profli.nep,profli.tep,...
        profli.qjli,profli.Raxe,profli.rmx,profli.spr,profli.vpr,geo.b0,cons.plh, ...
        option.xlh.*vt,option.dlh.*vt,1, ...
        directivity,0.5*vt,upshift,1,option.upshiftmode,option.npar_neg);
    
end


% calcul de la puissance electrique
rfan     = (3.56e6 + 14.03e6) ./ 3.56e6 ;
pout     = cons.plh + cons.pecrh + cons.picrh + real(cons.pnbi) + imag(cons.pnbi) + zs.pohm + option.mul_blanket .* rfan .* zs.pfus;
pelectot = option.carnot .* pout;
zs.pelec = pelectot .* (1 - option.aux)  - 1 ./ option.effinj .* (cons.plh + cons.pecrh + cons.picrh +  ...
    real(cons.pnbi) + imag(cons.pnbi));

% recalcul de betaptot pour eviter des propagation d'erreur
% betaptot n'est pas utilise dans le coeur de METIS
zs.betaptot  = zs.betap .* zs.w ./ zs.wth;

% fonction d'amortissement;
function out = zamor0(e1,e2,f)

out =[];

noms = fieldnames(e1);
for k = 1:length(noms)
    nomc = noms{k};
    %v1   = getfield(e1,nomc);
    %v2   = getfield(e2,nomc);
    v1   = e1.(nomc);
    v2   = e2.(nomc);
    s    = f .* v1 + (1-f) .* v2;
    
    % pour debuggage
    %if any(imag(s))
    %  fprintf('imag in %s\n',nomc);
    %  keyboard
    %end
    ind  = find(~isfinite(s));
    if ~isempty(ind)
        s(ind) = v1(ind);
    end
    %out = setfield(out,nomc,s);
    out.(nomc) = s;
end

% prend la derniere valeur et la recopie
function out = working_point(e1)

out =[];

noms = fieldnames(e1);
for k = 1:length(noms)
    nomc = noms{k};
    %v1   = getfield(e1,nomc);
    v1   = e1.(nomc);
    switch nomc
        case {'t','time','temps'}
            s = v1;
        otherwise
            if size(v1,1) > 1
                if all(isfinite(v1(end,:)))
                    s = ones(size(v1,1),1) * v1(end,:);
                else
                    s = v1;
                end
            else
                s = v1;
            end
            
            ind  = find(~isfinite(s));
            if ~isempty(ind)
                s(ind) = v1(ind);
            end
    end
    %out = setfield(out,nomc,s);
    out.(nomc) = s;
end


function x = znonan(x)

if all(~isfinite(x(:)))
    r = 0;
else
    r   = mean(x(isfinite(x)));
    if isempty(r)
        r = 0;
    end
end
x(~isfinite(x)) = r;


function e = z0delta(exp0d,sim,p,t,inter);

if nargin < 5
    inter =[];
end

if ~isempty(inter)
    ind = find((t>= min(inter)) & (t <= max(inter)));
    exp0d = exp0d(ind);
    sim   = sim(ind);
    p   	= p(ind);
    t  	 = t(ind);
end


e = (exp0d - sim) .* p;

indnan = find(~isfinite(e));
e(indnan) = [];
t(indnan) = [];
if isempty(e)
    e = inf;
else
    e = trapz(t,e);
end



function out = z0sgolayfilt(in,order,nbptal,w,d,tinter)

if isempty(in)
  out  = in;
  return
end

if nargin < 4
    w =[];
end
if nargin < 5
    d= 1;
end
if nargin < 6
    tinter = [];
end
%%fprintf(' %d',length(tinter));

if iscomplex(in)
    out =             z0sgolayfilt(real(in),order,nbptal,w,d,tinter) + ...
          sqrt(-1) .* z0sgolayfilt(imag(in),order,nbptal,w,d,tinter);
elseif all(in >=0)
    if ~isempty(tinter)
      out = max(min(in(:)),antialisasing_evolution(in,tinter));
    else
      out = max(min(in(:)),sgolayfilt(in,order,nbptal,w,d));
    end
elseif all(in <= 0)
    if ~isempty(tinter)
      out = min(max(in(:)),antialisasing_evolution(in,tinter));
    else
      out = min(max(in(:)),sgolayfilt(in,order,nbptal,w,d));  
    end 
else
    if ~isempty(tinter)
      out = antialisasing_evolution(in,tinter);
    else
      out = sgolayfilt(in,order,nbptal,w,d);
    end
end

function out = antialisasing_evolution(in,tinter)

persistent tinter_mem
persistent mat_inv

out = in;
tinter = tinter - tinter(1);
if isempty(tinter_mem) || any(tinter ~= tinter_mem) || isempty(mat_inv)
  %mat = cat(2,tinter .^ 2,tinter,ones(size(tinter)));
  mat = cat(2,tinter,ones(size(tinter)));
  mat_inv = pinv(mat);
  tinter_mem = tinter;
end
for k=1:size(in,2)
    %pp_alt = polyfit(tinter,in(:,k),2)
    pp = (mat_inv * in(:,k))';
    %out(:,k) = pp(1) .* (tinter .^ 2) + pp(2) .* tinter + pp(3);
    out(:,k) = pp(1) .* tinter + pp(2);
end

%figure(21);clf;
%plot(tinter,out,'ob',tinter,in,'r.')
%keyboard
%pause(0.1);
