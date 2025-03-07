%  memo E. Nardon, 11/02/2011 +
% ref : M. Nassi, Fusio Technology 24 (1993), p 50- et references.
% ref : P. Helander Collisional transport ... , p 162 - 165
% Calcule le flux consomme en fonction du temps.
% On utilise l'approximation Lsol (self du solenoide) = Msol,p (mutuelle entre le solenoide et le plasma) qui serait valable si le solenoide etait infini.

% breakdown : utilisation uniquement des donnees de METIS + coefficient de couplage
% le model a la bonne reponse en densite, tension par tour initiale et puissance de chauffage EC
%couplage = coefficient de transfert de l'energy entre le CS et le plasma en tenant compte de la disspation dans les coques, au claquage le temps de diffusion du courant dans le plasma est petit

function [flux_tot, flux_cs, flux_plasma_tot, flux_res] = z0fluxconsum(z0dinput, zerod, profil0d)

post.z0dinput = z0dinput;
post.zerod = zerod;
post.profil0d = profil0d;

switch upper(post.z0dinput.machine)
    case 'JET'
        couplage = 0.2;%(fer + coques resitives)
    case 'TS'
        couplage = 0.2; % proche de 0.2 d'apres les mesures (fer + coques resitives )
    case 'ITER'
        couplage = 0.15; % coques faible resistance
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
%             if ~exist('NOPLUTFLUX_HELIOS','var')
%                 prompt={sprintf('Coupling coefficient between CS/PF -> plasma for breakdown / burn-through\nJET/TS = 0.2, ITER = 0.15, optimized tokamak = 0.3')};
%                 name='Coupling coefficient';
%                 numlines=1;
%                 defaultanswer={sprintf('%g',couplage)};
%                 answer=inputdlg(prompt,name,numlines,defaultanswer);
%                 if ~isempty(answer)
%                     couplage = min(1,max(eps,str2num(answer{1})));
%                 end
%             end
        end
end
mu0 = 4*pi*10^(-7);
% model simple (conservation de l'energie)
eioniz = post.z0dinput.option.eioniz; % molecule  -> ions + electrons (par ion) + acceleration des electrons
if eioniz == 0
    eioniz = z0eioniz_div(post.zerod.tem(1),post.zerod.nem(1));
end
ploss_break = post.zerod.ploss(1);
if ploss_break < (post.zerod.ip(1)/1000)
    ploss_break = post.zerod.ploss(2);
end
dt_break =(post.zerod.vp(1) .* post.zerod.nem(1) .* 1.602176462e-19 .* eioniz + post.zerod.w(1)) ./ max(1,ploss_break);
breakdown_flux = post.zerod.vloop(1) .* dt_break ./ couplage;

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
%ip_dot = z0dxdt(post.zerod.ip,t);

% calcul du flux plasma
%flux_ext_dot  = 0.5 * L_plasma_dot.*post.zerod.ip + L_plasma.*ip_dot;
%flux_external_alt = cumtrapz(t ,flux_ext_dot) + L_plasma(1) .* post.zerod.ip(1);
% ameliore la precision
flux_external     = L_plasma .* post.zerod.ip - cumtrapz(t ,0.5 * L_plasma_dot.*post.zerod.ip);

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
flux_res = cumtrapz(t , pohm ./ post.zerod.ip) + breakdown_flux;
%flux_res_alt = cumtrapz(t , post.zerod.pohm ./ post.zerod.ip) + breakdown_flux;

% Flux totale du poloidal
% Ref : M. Massi, Fusion Technology vol 24 (1993) p 50-64
% flux total (avec champ vertical)
% c'est le flux qu'il faut contre balancer avec le systeme poloidale (flux CS  + flux champ vertical)
flux_tot     = flux_plasma   + flux_res + flux_bv;
flux_tot_alt = flux_psi_bord + flux_external  + flux_bv + breakdown_flux;

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
Cejima = flux_res(ind_ip_max) ./ mu0 ./ post.zerod.ip(ind_ip_max) ./ geo.R(ind_ip_max);

flux_plasma_tot = flux_plasma_li + flux_external;

% if ~exist('NOPLUTFLUX_HELIOS','var')

%     % calcul de la duree du choc
%     if isfield(post.z0dinput.option,'available_flux')
%         ind_duration     = find(flux_tot >= post.z0dinput.option.available_flux,1);
%         ind_duration_alt = find(flux_tot_alt >= post.z0dinput.option.available_flux,1);
%         duration = t(ind_duration);
%         duration_alt = t(ind_duration_alt);
%         % most probable
%         [nbin,ipbin] = hist(post.zerod.ip,max(3,length(post.zerod.ip)/10));
%         ipflatop = sum((nbin == max(nbin)) .* ipbin) ./ max(1, sum((nbin == max(nbin))));
%         ind_flattop      = find(post.zerod.ip >= ipflatop,1);
%         time_flattop     = t(ind_flattop);
%         duration_flattop = duration - time_flattop;
%         duration_flattop_alt = duration_alt - time_flattop;
%         flux_flow        = flux_tot(ind_duration) - flux_tot(ind_flattop);
%         flux_flow_alt    = flux_tot_alt(ind_duration_alt) - flux_tot_alt(ind_flattop);

%         fprintf('----------------------------------------------------------------------\n');
%         try
%             fprintf('%s@%d :\n',post.z0dinput.option.machine,post.z0dinput.option.shot);
%         catch
%             fprintf('%s@%d :\n',post.z0dinput.machine,post.z0dinput.shot);
%         end
%         fprintf('shot duration between %g and %g (s)\n',min(duration,duration_alt),max(duration,duration_alt));
%         fprintf('maximum flat-top duration between %g and %g (s)\n', ...
%             min(duration_flattop,duration_flattop_alt),max(duration_flattop,duration_flattop_alt));
%         fprintf('maximum flat-top duration between %g and %g (h)\n', ...
%             min(duration_flattop,duration_flattop_alt)/3600,max(duration_flattop,duration_flattop_alt)/3600);
%         fprintf('flux consumption between %g and %g (Wb/s)\n', ...
%             min(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt), ...
%             max(flux_flow./duration_flattop,flux_flow_alt ./ duration_flattop_alt));
%         fprintf('----------------------------------------------------------------------\n');

%     end


% end

