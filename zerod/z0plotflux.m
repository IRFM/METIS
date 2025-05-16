%  memo E. Nardon, 11/02/2011 +
% ref : M. Nassi, Fusio Technology 24 (1993), p 50- et references.
% ref : P. Helander Collisional transport ... , p 162 - 165
% Calcule le flux consomme en fonction du temps.
% On utilise l'approximation Lsol (self du solenoide) = Msol,p (mutuelle entre le solenoide et le plasma) qui serait valable si le solenoide etait infini.
% breakdown : utilisation uniquement des donnees de METIS + coefficient de couplage
% le model a la bonne reponse en densite, tension par tour initiale et puissance de chauffage EC
%couplage = coefficient de transfert de l'energy entre le CS et le plasma en tenant compte de la disspation dans les coques, au claquage le temps de diffusion du courant dans le plasma est petit
available_flux = Inf;
if isempty(post.z0dinput.machine)
    post.z0dinput.machine = 'unknown';
end
switch upper(post.z0dinput.machine)
    case 'JET'
        couplage = 0.2;%(fer + coques resitives)
        available_flux = 34;
    case {'TS','WEST'}
        couplage = 0.2; % proche de 0.2 d'apres les mesures (fer + coques resitives )
        available_flux = 9.8  - (-7.82);
    case 'ITER'
        couplage = 0.15; % coques faible resistance
        available_flux = 227;
    case {'JT60-SA','JT-60SA'}
        couplage = 0.15; % coques faible resistance
        available_flux = 41;
    otherwise
        if exist('couplage','var')
            % rien
        else
            % couplage = 1/3; % bon tokamak (cf MHD & FUSION Magnets de RJ Thome)
            if ~isempty(findstr('DEMO',upper(post.z0dinput.machine))) || ~isempty(findstr('REACT',upper(post.z0dinput.machine)))
                couplage = 0.15;
            else
                couplage = 1/3; % bon tokamak (cf MHD & FUSION Magnets de RJ Thome), valeur conservatrice.
            end
            if ~exist('NOPLUTFLUX_HELIOS','var')
                prompt={sprintf('Coupling coefficient between CS/PF -> plasma for breakdown / burn-through\nJET/TS = 0.2, ITER = 0.15, optimized tokamak = 0.3')};
                name='Coupling coefficient';
                numlines=1;
                defaultanswer={sprintf('%g',couplage)};
                answer=inputdlg(prompt,name,numlines,defaultanswer);
                if ~isempty(answer)
                    couplage = min(1,max(eps,str2num(answer{1})));
                end
            end
        end
end
if isfield(post.z0dinput.option,'available_flux') && ~isfinite(post.z0dinput.option.available_flux) && isfinite(available_flux)
    post.z0dinput.option.available_flux = available_flux;
end
mu0 = 4*pi*10^(-7);
% model simple (conservation de l'energie)
eioniz = post.z0dinput.option.eioniz; % molecule  -> ions + electrons (par ion) + acceleration des electrons
if eioniz == 0
    eioniz = z0eioniz_div(post.zerod.tem(1),post.zerod.nem(1));
end
ploss_break  =  max(post.profil0d.qe(1,:),[],2) + max(post.profil0d.qi(1,:),[],2);
%ploss_break = max(post.zerod.ploss(1),post.zerod.pin(1) - post.zerod.prad(1));
if ploss_break < (post.zerod.ip(1) ^ 2 .* post.zerod.RR(1) ./ 1000)
    disp('unable to compute ploss at start of plasma');
    dt_break = NaN;
    breakdown_flux = 0;
else
    dt_break =(post.zerod.vp(1) .* post.zerod.nem(1) .* 1.602176462e-19 .* eioniz + post.zerod.w(1)) ./ max(1,ploss_break);
    breakdown_flux = post.zerod.vloop(1) .* dt_break ./ couplage;
end
% dans ce cas le calcul METIS inclu le breakdown
if post.z0dinput.option.berror ~= 0
    if isfield(post.zerod,'eddy_current') && isfield(post.zerod,'flux_edge_cor')
        [void_taue,void_f,ieddy,void_nec,void_tref,void_prf,flux_bord_cor,temps_ed] = z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons, ...
            post.zerod, post.zerod.taue,post.zerod.eddy_current,post.zerod.flux_edge_cor);
    else
        [void_taue,void_f,ieddy,void_nec,void_tref,void_prf,flux_bord_cor,temps_ed] = z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons, ...
            post.zerod, post.zerod.taue);
    end
    breakdown_flux = - max(0,interp1(post.z0dinput.cons.temps,flux_bord_cor,temps_ed,'linear','extrap'));
    dt_break = post.z0dinput.cons.temps(1) - temps_ed;
end

t    = post.zerod.temps;
geo  = post.z0dinput.geo;
exp0d = post.z0dinput.exp0d;
if isfield(exp0d,'edgeflux')
    exp0d.flux = exp0d.edgeflux;
end
Raxe    = (geo.R + post.zerod.d0);
Ka   = sqrt(post.zerod.sp ./ pi) ./ geo.a;
% Self externe du plasma (formule approchee pour un profil de courant plat) :
L_plasma = mu0 * geo.R .* (log(8 *  geo.R ./ ((post.zerod.sp/pi).^0.5)) - 2);
[Lext,MBv] = z0lextmutbv(geo.R,geo.a,Ka);
% papier de reference :
% Hirshman and Neilson Phys. of Fluids 29 (3) 1986, p 792
% S. Y. Mak and K. young in Phys. Educ. 21 (1986) p111-
% L_plasma ets plus precise que Lext !
% pour MBV on n'a pas le choix
L_plasma_dot = z0dxdt(L_plasma,t);
ip_dot = z0dxdt(post.zerod.ip,t);

% calcul du flux plasma
%flux_ext_dot  = 0.5 * L_plasma_dot.*post.zerod.ip + L_plasma.*ip_dot;
%flux_external_alt = cumtrapz(t ,flux_ext_dot) + L_plasma(1) .* post.zerod.ip(1);
% ameliore la precision ( a partir de d/dt(1/2 L*I*I) = dflux/dt * I)
flux_external     = L_plasma .* post.zerod.ip - cumtrapz(t ,0.5 * L_plasma_dot.*post.zerod.ip);
%flux_external_test  = L_plasma(1) .* post.zerod.ip(1) + cumtrapz(t,L_plasma .* ip_dot)  + cumtrapz(t ,0.5 * L_plasma_dot.*post.zerod.ip);
%figure;plot(t,flux_external,'r',t,flux_external_test,'b');
%keyboard

% flux manquant dans la self du a li (il faut reprendre la normalisation de metis!)
L_i = mu0 * geo.R .* post.zerod.li / 2;
L_i_dot  = z0dxdt(L_i,t);
%flux_plasma_li = cumtrapz(t , 0.5 * L_i_dot .* post.zerod.ip + L_i .* ip_dot) + L_i(1) .* post.zerod.ip(1) ;
% ameliore la precision
flux_plasma_li = L_i .* post.zerod.ip - cumtrapz(t , 0.5 * L_i_dot .* post.zerod.ip);
flux_plasma_li = max(0,flux_plasma_li); % probleme de presision dans la descente courant
%
flux_plasma = flux_external + flux_plasma_li;

% flux_plamsa_int
flux_plasma_int = interp1(post.profil0d.temps,post.profil0d.psi(:,1) - post.profil0d.psi(:,end),t,'pchip','extrap') .* 2 .* pi;

% flux au bord du plasma
flux_psi_bord    = - interp1(post.profil0d.temps,post.profil0d.psi(:,end) - post.profil0d.psi(1,1),t,'pchip','extrap') .* 2 .* pi;
flux_psi_centre  = - interp1(post.profil0d.temps,post.profil0d.psi(:,1) - post.profil0d.psi(1,1),t,'pchip','extrap') .* 2 .* pi;

% calcul de la diff�rence via l'inductance interne (papier Ejima). definition differente de li
%delta_psi = L_i .* post.zerod.ip + 2/3 .* cumtrapz(t ,L_i .* ip_dot);
%delta_psi = delta_psi - delta_psi(1) + flux_psi_bord(1) - flux_psi_centre(1);
%figure;plot(t,delta_psi,'r',t,flux_psi_bord - flux_psi_centre,'b',t,flux_psi_bord - flux_psi_centre-delta_psi,'k', ...
%            t,cumtrapz(t ,L_i .* ip_dot),'m',t,cumtrapz(t ,L_i_dot .* post.zerod.ip),'c');drawnow
%keyboard

%flux du champ vertical
% calcul Wesson p 120-123.
bv =  mu0 .* post.zerod.ip ./ (4 .* pi .* geo.R) .* (log(8 .* geo.R ./((post.zerod.sp/pi).^0.5)) + post.zerod.betaptot + post.zerod.li ./ 2 - 3/2);
% ref : M. Nassi
flux_bv =  (1 - MBv) .* pi .* geo.R .^ 2 .* bv;

% flux inductif
flux_induc = cumtrapz(t ,post.zerod.vloop) + flux_plasma_int(1);

% ref :Ohmic flux consumption during initial operation of the NSTX spherical torus
% J.E. Menard et al, NF 2001 41 1197
% recalcul de Pohm sans les securite pour la convergence de METIS
pohm     = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.ej,2);
pohm     = interp1(post.profil0d.temps,pohm,t,'pchip',0);
indbad   = find(abs(pohm) > abs(pi .*  post.zerod.pohm));
if ~isempty(indbad)
    pohm(indbad) =  post.zerod.pohm(indbad);
end
flux_res = cumtrapz(t , pohm ./ post.zerod.ip) + breakdown_flux;
%flux_res_alt = cumtrapz(t , post.zerod.pohm ./ post.zerod.ip) + breakdown_flux;

% Flux totale du poloidal
% Ref : M. Nassi, Fusion Technology vol 24 (1993) p 50-64
% flux total (avec champ vertical)
% c'est le flux qu'il faut contre balancer avec le systeme poloidale (flux CS  + flux champ vertical)
flux_tot     = flux_plasma   + flux_res + flux_bv;
flux_tot_alt = flux_psi_bord + flux_external  + flux_bv + breakdown_flux;

% for verification
if 0
    vloop = post.zerod.vloop;
    vloop(find(~isfinite(vloop))) = sqrt(-1);
    flux_tot_nardon = cumtrapz(t ,vloop) + cumtrapz(t,L_plasma .* ip_dot)  + cumtrapz(t ,L_plasma_dot.*post.zerod.ip)+ breakdown_flux;
    flux_tot_nardon = flux_tot_nardon - flux_tot_nardon(1) + flux_tot(1);
    flux_tot_nardon_alt = cumtrapz(t,(L_plasma + L_i) .* ip_dot)  + 0.5 .* cumtrapz(t ,(L_plasma_dot + L_i_dot) .* post.zerod.ip) + flux_res;
    flux_tot_nardon_alt = flux_tot_nardon_alt - flux_tot_nardon_alt(1) + flux_tot(1);
    figure;plot(t,flux_tot,'r',t,flux_tot_alt,'b',t,flux_tot_nardon,'m',t,flux_tot_nardon_alt,'c', ...
        t,flux_tot - flux_tot_nardon,t,flux_plasma + flux_res,'k',t,flux_bv,'g');
    keyboard;
end
% calcul de flux CS
flux_cs     = flux_tot     -  pi .* geo.R .^ 2 .* bv;
flux_cs_alt = flux_tot_alt -  pi .* geo.R .^ 2 .* bv;

% flux resistif alternatif
flux_res_ri = flux_psi_bord - flux_plasma_li + breakdown_flux;

% precalcul : flux mesurer a la surface du plasma
vloop = post.zerod.vloop;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = cumtrapz(t,vloop,1);
flux(find(imag(vloop))) = NaN;
flux_loop = real(flux) + breakdown_flux + flux_plasma_int(1);

% precalcul : flux mesurer par une boucle fixe
vloop = post.zerod.vmes;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = cumtrapz(t,vloop,1);
flux(find(imag(vloop))) = NaN;
flux_meas = real(flux) + breakdown_flux + flux_plasma_int(1);


% flux resistif de l'axe
% le plus simple est d'utiliser le flux interne
flux_res_axe = flux_psi_centre ;
% offset correction
flux_res_axe  = flux_res_axe - flux_res_axe(1) + breakdown_flux;


% corection flux vertical (car pas precis)
%  fdot_res_axe   = z0dxdt(flux_res_axe,t);
%  flux_res_axe   = max(0,flux_res_axe(1)) + cumtrapz(t,max(0,fdot_res_axe));
%  flux_bv_cor    = flux_bv_ref - flux_res_axe;
%  MBv_cor        = 1 - flux_bv_cor ./ ( pi .* geo.R .^ 2 .* bv);

% correction offset flux_cons (phase non dercrite par la simulation)
%flux_cons       = flux_cons + flux_plasma_li(1) + flux_res_axe(1);
%flux_meas       = flux_meas + flux_plasma_li(1) + breakdown_flux + flux_res_axe(1);
%flux_loop       = flux_loop + flux_plasma_li(1) + breakdown_flux + flux_res_axe(1);
%flux_psi_bord   = flux_psi_bord  + flux_plasma_int(1) + breakdown_flux + flux_res_axe(1);

% flux resitif (par soustraction)
%flux_res_ri =  flux_tot - flux_external - flux_plasma_li  + (flux_res_axe - breakdown_flux) - flux_plasma_int(1);

% calcul de Cejima
dd = abs(post.zerod.ip - max(post.zerod.ip)) ./ max(post.zerod.ip);
ind_ip_max = find(dd  <= 1e-2,1);
Cejima_db = flux_res(ind_ip_max) ./ mu0 ./ post.zerod.ip(ind_ip_max) ./ geo.R(ind_ip_max);
Cejima = (flux_res(ind_ip_max) - breakdown_flux) ./ mu0 ./ post.zerod.ip(ind_ip_max) ./ geo.R(ind_ip_max);


if ~exist('NOPLUTFLUX_HELIOS','var')
    
    % calcul de la duree du choc
    if isfield(post.z0dinput.option,'available_flux') && isfinite(post.z0dinput.option.available_flux)
        ind_duration     = find(flux_tot >= post.z0dinput.option.available_flux,1);
        ind_duration_alt = find(flux_tot_alt >= post.z0dinput.option.available_flux,1);
        duration = t(ind_duration);
        duration_alt = t(ind_duration_alt);
        % most probable
        [nbin,ipbin] = hist(post.zerod.ip,max(3,length(post.zerod.ip)/10));
        ipflatop = sum((nbin == max(nbin)) .* ipbin) ./ max(1, sum((nbin == max(nbin))));
        ind_flattop      = find(post.zerod.ip >= ipflatop,1);
        time_flattop     = t(ind_flattop);
        duration_flattop = duration - time_flattop;
        duration_flattop_alt = duration_alt - time_flattop;
        flux_flow        = flux_tot(ind_duration) - flux_tot(ind_flattop);
        flux_flow_alt    = flux_tot_alt(ind_duration_alt) - flux_tot_alt(ind_flattop);
        
        fprintf('----------------------------------------------------------------------\n');
        try
            fprintf('%s@%d :\n',post.z0dinput.option.machine,post.z0dinput.option.shot);
        catch
            fprintf('%s@%d :\n',post.z0dinput.machine,post.z0dinput.shot);
        end
        fprintf('shot duration between %g and %g (s)\n',min(duration,duration_alt),max(duration,duration_alt));
        fprintf('maximum flat-top duration between %g and %g (s)\n', ...
            min(duration_flattop,duration_flattop_alt),max(duration_flattop,duration_flattop_alt));
        fprintf('maximum flat-top duration between %g and %g (h)\n', ...
            min(duration_flattop,duration_flattop_alt)/3600,max(duration_flattop,duration_flattop_alt)/3600);
        fprintf('flux consumption between %g and %g (Wb/s)\n', ...
            min(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt), ...
            max(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt));
        fprintf('----------------------------------------------------------------------\n');
        
        
        
        
    end
    
    fullscreen = get(0,'ScreenSize');
    h = findobj(0,'type','figure','tag','z0plotflux');
    if isempty(h)
        h=figure('tag','z0plotflux');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen)
    
    hold on
    %plot(t,flux_cons,'b',t,flux_tot,'r',t,flux_plasma_li + flux_external,'m',t,flux_res,'k',t,flux_res_ri,':k',t,flux_res_axe,'k-.');
    plot(t,flux_cs,'b',t,flux_tot,'r',t,flux_plasma_li + flux_external,'m',t,flux_plasma_li,'c', t,flux_external,'c-.',...
         t,flux_res,'k',t,flux_res_axe,'k-.');
    leg_exp = 1;
    if post.z0dinput.mode_exp == 0
        fluxexp = -post.profil0d.psi(:,end) .* 2 .* pi;
        flux0 = mean(fluxexp(isfinite(fluxexp))) - mean(flux(isfinite(flux_loop)));
        plot(t,flux_loop,'g',t,flux_meas,'g-.',post.profil0d.temps,fluxexp - flux0 ,'c');
    else
        if isfield(exp0d,'flux')
            flux0 = -mean(flux_loop(isfinite(flux_loop))) / 2 / pi - mean(real(exp0d.flux(isfinite(exp0d.flux))));
        else
            flux0 = 0;
            exp0d.flux = NaN .* t;
        end
        if all(exp0d.flux == 0) ||all(~isfinite(exp0d.flux))
            leg_exp = 0;
        end
        plot(t,flux_loop,'g',t,flux_meas,'g-.',t,-(real(exp0d.flux) + flux0) .* 2 .* pi,'c');
    end
    grid
    xlabel('Time (s)' )
    ylabel('flux (Wb)' )
    if isfinite(dt_break)
        if dt_break >= 0
            bdcsf = flux_cs(1);
        else
            bdcsf = interp1(post.z0dinput.cons.temps,flux_cs,temps_ed,'linear','extrap');
        end
    else
        bdcsf = NaN;
    end
    if isfinite(bdcsf)
        try
            title(sprintf('METIS : %s@%d/FLUX (breakdown CS flux = %g Wb, C_{Ejima} = %g, C_{Ejima + bd} = %g)',post.z0dinput.option.machine,post.z0dinput.option.shot,bdcsf,Cejima,Cejima_db));
        catch
            title(sprintf('METIS : %s@%d/FLUX (breakdown CS flux = %g Wb, C_{Ejima} = %g, C_{Ejima + bd} = %g)',post.z0dinput.machine,post.z0dinput.shot,bdcsf,Cejima,Cejima_db));
        end
    else
        try
            title(sprintf('METIS : %s@%d/FLUX (breakdown failed, C_{Ejima} = %g)',post.z0dinput.option.machine,post.z0dinput.option.shot,Cejima));
        catch
            title(sprintf('METIS : %s@%d/FLUX (breakdown failed, C_{Ejima} = %g)',post.z0dinput.machine,post.z0dinput.shot,Cejima));
        end
    end
    if leg_exp
        legend('Total CS Flux','Total flux (CS + vertical field)','Total plasma flux (L_p * I_p)', ...
            'Internal flux (L_i * I_p)', 'External flux (L_{ext} * I_p)', ...
            'Resitive consumed flux (int(t,Pohm / Ip))','Resistive flux @ magnetic axis',...
            'Computed flux loop (plasma surface)','Computed flux loop (fixed loop)','Measured flux loop','location','best')
    else
        legend('Total CS Flux','Total flux (CS + vertical field)','Total plasma flux (L_p * I_p)', ...
            'Internal flux (L_i * I_p)', 'External flux (L_{ext} * I_p)', ...
            'Resitive consumed flux (int(t,Pohm / Ip))','Resistive flux @ magnetic axis', ...
            'Computed flux loop (plasma surface)','Computed flux loop (fixed loop)','location','best')
    end
    
    if verLessThan('matlab', 'R2015b')
        xx=[t(:);t(end:-1:1)];
        yy=[max(flux_cs,flux_cs_alt);min(flux_cs(end:-1:1),flux_cs_alt(end:-1:1))];
        hl=patch(xx,yy,[0.7,0.7,1],'erasemode','xor','edgecolor','none');
        yy=[flux_induc + breakdown_flux;flux_psi_bord(end:-1:1) + breakdown_flux];
        hl=patch(xx,yy,[0.7,1,0.7],'erasemode','xor','edgecolor','none');
        yy=[max(flux_res,flux_res_ri);min(flux_res(end:-1:1),flux_res_ri(end:-1:1))];
        hl=patch(xx,yy,[0.7,0.7,0.7],'erasemode','xor','edgecolor','none');
        yy=[max(flux_tot,flux_tot_alt);min(flux_tot(end:-1:1),flux_tot_alt(end:-1:1))];
        hl=patch(xx,yy,[1,0.7,0.7],'erasemode','xor','edgecolor','none');
    else
        xx=[t(:);t(end:-1:1)];
        yy=[max(flux_cs,flux_cs_alt);min(flux_cs(end:-1:1),flux_cs_alt(end:-1:1))];
        hl=patch(xx,yy,[0.7,0.7,1],'edgecolor','none');
        alpha(hl,0.7);
        yy=[flux_induc + breakdown_flux;flux_psi_bord(end:-1:1) + breakdown_flux];
        hl=patch(xx,yy,[0.7,1,0.7],'edgecolor','none');
        alpha(hl,0.7);
        yy=[max(flux_res,flux_res_ri);min(flux_res(end:-1:1),flux_res_ri(end:-1:1))];
        hl=patch(xx,yy,[0.7,0.7,0.7],'edgecolor','none');
        alpha(hl,0.7);
        yy=[max(flux_tot,flux_tot_alt);min(flux_tot(end:-1:1),flux_tot_alt(end:-1:1))];
        hl=patch(xx,yy,[1,0.7,0.7],'edgecolor','none');
        alpha(hl,0.7);
    end
    if dt_break < 0
        plot([t(1)-dt_break,t(1)],[breakdown_flux,breakdown_flux],'--g',[t(1)-dt_break,t(1)],[breakdown_flux,breakdown_flux],'ok');
    else
        plot([t(1)-dt_break,t(1)],[breakdown_flux,breakdown_flux+ flux_plasma_int(1)],'--g',[t(1)-dt_break,t(1)],[breakdown_flux,breakdown_flux+ flux_plasma_int(1)],'ok');
    end
    plot(t(ind_ip_max),flux_tot(ind_ip_max),'om');
    
    if isfield(post.z0dinput.option,'available_flux') && isfinite(post.z0dinput.option.available_flux)
        plot(min(t),post.z0dinput.option.available_flux,'r*');
        if ~isempty(duration) && ~isempty(duration_alt)
            plot([min(duration,duration_alt),max(duration,duration_alt)],[0,0],'*r');
        elseif ~isempty(duration)
            plot(duration,0,'*r');
        end
        
        ylim = get(gca,'ylim');
        xlim = get(gca,'xlim');
        xtxt = min(xlim) + 0.03 .* (max(xlim) - min(xlim));
        ytxt = max(ylim) - 0.05 .* (max(ylim) - min(ylim));
        text(xtxt,ytxt,sprintf('shot duration between %g and %g (s)\n',min(duration,duration_alt),max(duration,duration_alt)),'fontsize',12,'fontweight','bold','fontname','times');
        ytxt = max(ylim) - 0.1 .* (max(ylim) - min(ylim));
        text(xtxt,ytxt,sprintf('maximum flat-top duration between %g and %g (s)\n', ...
            min(duration_flattop,duration_flattop_alt),max(duration_flattop,duration_flattop_alt)),'fontsize',12,'fontweight','bold','fontname','times');
        ytxt = max(ylim) - 0.15 .* (max(ylim) - min(ylim));
        text(xtxt,ytxt,sprintf('maximum flat-top duration between %g and %g (h)\n', ...
            min(duration_flattop,duration_flattop_alt)/3600,max(duration_flattop,duration_flattop_alt)/3600),'fontsize',12,'fontweight','bold','fontname','times');
        ytxt = max(ylim) - 0.2 .* (max(ylim) - min(ylim));
        text(xtxt,ytxt,sprintf('flux consumption between %g and %g (Wb/s)\n', ...
            min(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt), ...
            max(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt)),'fontsize',12,'fontweight','bold','fontname','times');
        
    end
    set(h,'Papertype','A4','paperOrientation','landscape','PaperPositionMode','auto')
    
    edition2
    
end


