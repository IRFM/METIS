% script d'estimation des flux de matiere
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3
 
% compatibilite ascendante
if ~isfield(post.zerod,'nbar_nat') || all(~isfinite(post.zerod.nbar_nat))
     ulh = 0.25;
     fh  = post.zerod.modeh;
     post.zerod.nbar_nat   = min(post.zerod.negr,max(1e13, 1e20 .* (post.z0dinput.geo.b0 ./ post.zerod.q95 ./ post.z0dinput.geo.R) .^ 0.6  .*  (ulh + (1 - ulh) .* fh)));
     pellet_fraction = post.zerod.frac_pellet .* (post.zerod.frac_pellet > 1e-2);
else
     pellet_fraction = post.zerod.frac_pellet;
end

% choix du coefficiet de recyclage
if exist('NOINTER','var')
	    eta_p = 0.5;
	    filter_width = 11;
else
    prompt={'pellet fuelling efficiency','filter width (number of points)'};
    name='METIS gas bilan';
    numlines=1;
    defaultanswer={'0.5','11'};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answer)
	    eta_p = min(0.95,max(1e-3,str2num(answer{1})));
	    filter_width = fix(str2num(answer{2}));
    else
	    eta_p = 0.5;
	    filter_width = 11;
    end
end
Recycling = post.z0dinput.option.Recycling;
eta_g_default = 0.1;

% temps
mat.temps        = post.zerod.temps;
% cette formule n'est pas precise numeriquement
% flux de matiere sortant
mat.output      = interp1(post.profil0d.temps,post.profil0d.ge(:,end) .* post.profil0d.grho2(:,end) .* post.profil0d.vpr_tor(:,end), ...
			   post.zerod.temps,'pchip','extrap') ./ (4.41e-4 .* phys.avo);
% bilan contenu du plasma
ntot              = trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.nep,2) ./ (4.41e-4 .* phys.avo);
mat.contents      = interp1(post.profil0d.temps,ntot,post.zerod.temps,'pchip','extrap');
mat.plasmanetflux = z0dxdt(mat.contents,post.zerod.temps); % flux net qui sort du plasam

% source de recylcage dans le plasma (a partir de la source de neutre qui entre dans le plasam)
mat.s0_in       = post.zerod.n0a ./ (4.41e-4 .* phys.avo);  % recyclage  + gas puff

% effect a recycling fraction (at LCFS not in the divertor)
% sorte de minimum pour le flux de gaz qui circule dans la SOL
switch post.z0dinput.option.configuration
case {0,1}
    mat.s0       = mat.s0_in  ./ max(eps,post.z0dinput.option.fn0a); 
case {2,3}
    mat.s0       = mat.s0_in  ./ max(eps,post.zerod.xpoint .* post.z0dinput.option.fn0a_div + (~post.zerod.xpoint) .* post.z0dinput.option.fn0a);
otherwise
    mat.s0       = mat.s0_in  ./ max(eps,post.z0dinput.option.fn0a_div); 
end
% s0_star est la source qui sort du plasma
% computation of suggested fn0a
fn0a_div = mat.output(post.zerod.xpoint ~= 0) ./ mat.s0(post.zerod.xpoint ~= 0);
fn0a_div = mean(fn0a_div(isfinite(fn0a_div) & (fn0a_div>0)));
fn0a     = mat.output(post.zerod.xpoint == 0) ./ mat.s0(post.zerod.xpoint == 0);
fn0a     = mean(fn0a(isfinite(fn0a) & (fn0a>0)));
disp('Suggested setting for next METIS run:');
disp('(must be between 0 and 1 to be valid)');
fn0a
fn0a_div

if filter_width > 1
	mat.s0          = medfilt1(mat.s0,filter_width);
	mat.s0_in       = medfilt1(mat.s0_in,filter_width);
	mat.output      = medfilt1(mat.output ,filter_width);
end
mat.bilan_s0    = cumtrapz(post.zerod.temps,mat.s0);
mat.bilan_s0_in = cumtrapz(post.zerod.temps,mat.s0_in);
% security
output_mem       = mat.output;
mat.output       = max(mat.output,1.1 .* mat.s0);
mat.bilan_output = cumtrapz(post.zerod.temps,mat.output);

% source due au glacon (dans le plasma)
mat.pellet_star   = interp1(post.profil0d.temps,trapz(post.profil0d.xli,post.profil0d.spellet .* post.profil0d.vpr,2), ...
				post.zerod.temps,'pchip','extrap')  ./ (4.41e-4 .* phys.avo);
mat.pellet_star(~isfinite(mat.pellet_star)) = 0;
if (filter_width > 1) &&  (post.z0dinput.option.pif < 1)
	mat.pellet_star          = medfilt1(mat.pellet_star,filter_width);
end
mat.bilan_pellet_star  = max(eps,cumtrapz(post.zerod.temps,mat.pellet_star));
% flux de  matiere du a l'injection de gla�on
mat.pellet = mat.pellet_star ./ eta_p;
mat.bilan_pellet = mat.bilan_pellet_star ./ eta_p;

% source due a l'injection de neutres
mat.snbi_star        = real(post.zerod.pnbi) ./ post.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo) + ...
                       imag(post.zerod.pnbi) ./ post.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo);
if filter_width > 1
	mat.snbi_star          = medfilt1(mat.snbi_star,filter_width);
end
mat.bilan_nbi_star   = max(eps,cumtrapz(post.zerod.temps,mat.snbi_star));
% flux de matiere du a l'injection de neutres
mat.snbi             = real(post.zerod.pnbi)  ./ post.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(post.zerod.frnbi)) + ...
                       imag(post.zerod.pnbi)  ./ post.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(post.zerod.frnbi));
if filter_width > 1
	mat.snbi          = medfilt1(mat.snbi,filter_width);
end
mat.perte_snbi       = real(post.zerod.pnbi)  ./ post.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(post.zerod.frnbi)) .* (1 - real(post.zerod.frnbi)) + ...
                       imag(post.zerod.pnbi)  ./ post.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(post.zerod.frnbi)) .* (1 - imag(post.zerod.frnbi));                      
if filter_width > 1
	mat.perte_snbi          = medfilt1(mat.perte_snbi,filter_width);
end
mat.bilan_nbi   = max(eps,cumtrapz(post.zerod.temps,mat.snbi));

% calcul gaz puff (a partir du scaling pour la densite naturelle)
fact_nbar = post.zerod.nem ./ max(1,post.zerod.nbar); 
if (post.z0dinput.option.tauhemul < 0) && (post.z0dinput.option.Recycling < 1)
        % modele qui prend en compte le confinement reel et le recyclage dans le divertor
	% ref Stangeby section 6.7 
        % ref originale : D. Reiter et al, PPCF vol 33 (1991) p 1579-1600
	tau_ref    = post.zerod.tauhe - post.z0dinput.option.Recycling ./ (1 - post.z0dinput.option.Recycling) .* post.zerod.taup;
else
	tau_ref    = post.zerod.tauhe;
end
% gaz fuelling effciency
mat.eta_g_low  =  post.zerod.taup ./ tau_ref;
% source correspondant a la difference entre la densit� et la densit� naturelle du plasma
mat.dfuelling_bilan_dt  = fact_nbar .* (post.zerod.nbar - post.zerod.nbar_nat) ./ tau_ref .* post.zerod.vp ./ (4.41e-4 .* phys.avo) - ...
                          (mat.snbi_star + mat.pellet_star);
if filter_width > 1
	mat.dfuelling_bilan_dt     = medfilt1(mat.dfuelling_bilan_dt,filter_width);
end
mat.fuelling_bilan         = cumtrapz(post.zerod.temps,mat.dfuelling_bilan_dt);
% source de gaz puff dans le plasma si efficacite de 1
mat.gas_puff          =  mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt > 0);
mat.wall_pumping      = -mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt < 0);
if ~isfinite(fn0a)
  fn0a = post.z0dinput.option.fn0a;
end
if ~isfinite(fn0a_div)    
  fn0a_div = post.z0dinput.option.fn0a_div;
end
switch post.z0dinput.option.configuration
case {0,1}
    mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a); 
case {2,3}
    mat.gas_puff      = mat.gas_puff ./ max(eps,post.zerod.xpoint .* fn0a_div + (~post.zerod.xpoint) .* fn0a);
otherwise
    mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a_div); 
end

% source totale dans le plasma
mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
% bilan = input - losses = 0 + changement de densite du plasma
mat.bilan = mat.input - mat.plasmanetflux;
% correction de gaz puff
mat.gas_puff = mat.gas_puff - 1.1 .* mat.bilan .* (mat.bilan < 0);
%mise a jour  bilan = input - losses = 0 + changement de densite du plasma
mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
mat.bilan = mat.input - mat.plasmanetflux;
% pumping
mat.pumping = max(mat.bilan,mat.wall_pumping);
mat.bilan_pumping = cumtrapz(post.zerod.temps,mat.pumping);
% correction to gas puff
puff_correc  = (mat.pumping - mat.bilan);
mat.gas_puff = mat.gas_puff + puff_correc .* (puff_correc> 0);
mat.bilan_gas_puff = cumtrapz(post.zerod.temps,mat.gas_puff);
mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
% bilan final
mat.bilan = mat.input - mat.plasmanetflux  -mat.pumping;
mat.bilan_input   = cumtrapz(post.zerod.temps,mat.input);


% flux to divertor/limiter
mat.flux_divlim_total = mat.output + mat.gas_puff  + (1 - eta_p) .* mat.pellet + mat.perte_snbi - mat.s0_in;
mat.recycle           = mat.flux_divlim_total - mat.pumping;
mat.Reff              = mat.recycle ./ mat.flux_divlim_total;
mat.bilan_recycle     = cumtrapz(post.zerod.temps,mat.recycle);
mat.eta_g_high        = mat.s0_in ./ max(eps,mat.output + mat.gas_puff + mat.recycle);

%
tag=sprintf('z0plotgaz');
h = findobj(0,'type','figure','tag',tag);
if isempty(h)
	h=figure('tag',tag);
else
	figure(h);
end   
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultaxeslinewidth',3,'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',12,'Position',[400 40 1024 960],...
'PaperPositionMode','auto')
k=3;
subplot(k,1,1)
plot(mat.temps,mat.recycle ./ 10,mat.temps,mat.gas_puff,mat.temps,mat.pellet,...
     mat.temps,mat.snbi,mat.temps,mat.pumping,mat.temps,mat.wall_pumping,'-.');
title(sprintf('Zerod : %s@%d/gaz bilan',post.z0dinput.machine,post.z0dinput.shot));
legend('Recycle / 10','Gas puff','Pellet','NBI','Pumping','Wall_{pumping} (possiblility, if not saturated)','Location','NorthWest','orientation','horizontal')
%legend(gca,'boxoff')
ylabel('Pa.m^3/s (electrons)')
z0loglin(gca);

subplot(k,1,2)
semilogy(mat.temps,mat.bilan_output,mat.temps,mat.bilan_gas_puff + mat.bilan_pellet + mat.bilan_nbi,mat.temps,mat.contents);
legend('LCFS flux integral','Fuelling','Plasma contents','Location','NorthWest','orientation','horizontal')
ylabel('Pa.m^3 (electrons)')
set(gca,'ylim',[min(mat.contents),Inf]);
z0loglin(gca);

subplot(k,1,3)
plot(mat.temps,mat.Reff,'r',mat.temps,max(mat.eta_g_low,mat.eta_g_high),'b',mat.temps,min(mat.eta_g_low,mat.eta_g_high),'c');
legend('Recycling coefficient','Gas puff efficiency (high estimation)','Gas puff efficiency (low estimation)')
set(gca,'ylim',[0 1.2]);
z0loglin(gca);
setappdata(h,'Already_resize',true);
joint_axes(h,k);


% not interresting for other mixture
if post.z0dinput.option.gaz == 3
    [burn_fraction,averaged_burn_fraction,source_T_plasma,source_T_pellet,source_T_gaz,source_T_nbi, ...
                                    source_leakage,source_recycling]=z0burn_fraction(post,mat.eta_g_high,eta_p,0);
    fprintf('Averaged burn fraction = %g %%\n',averaged_burn_fraction .* 100);

    kiloT  = cumtrapz(post.z0dinput.cons.temps,(source_T_pellet + source_T_gaz + source_T_nbi) .* 4.41e-4 .* phys.avo .* 3.016 .* phys.ua) +  ...
                    post.zerod.nTm(1) .*  post.zerod.vp(1) .* 3.016 .* phys.ua;
                    
    gramme_seconde  = (source_T_pellet + source_T_gaz + source_T_nbi) .* 4.41e-4 .* phys.avo .* 3.016 .* phys.ua .* 1e3;
    averaged_gramme_seconde = max(kiloT) ./ (max(post.z0dinput.cons.temps) - min(post.z0dinput.cons.temps)) .* 1000;
    
    % sources - (internal plasma consumption + T plasma storage)
    divertor_flux   = - (source_T_plasma - source_leakage + source_recycling) + source_T_pellet + source_T_gaz + source_T_nbi - source_leakage + source_recycling;
    
    tag=sprintf('burn_fraction');
    h = findobj(0,'type','figure','tag',tag);
    if isempty(h)
	    h=figure('tag',tag);
    else
	    figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultaxeslinewidth',3,'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',12,'Position',[400 40 1024 960],...
    'PaperPositionMode','auto')
    
    k = 3;
    subplot(k,1,1)     
    plot(post.z0dinput.cons.temps,max(0,min(1,burn_fraction)) .* 100,post.z0dinput.cons.temps,max(0,gramme_seconde));
    legend(sprintf('Burn Fraction (%% , averaged burn fraction = %g %%)',averaged_burn_fraction .* 100), ...
           sprintf('Input source (g/s, averaged consumption = %g g/s)',averaged_gramme_seconde));
    xlabel('time (s)');
    title(sprintf('METIS : %s@%d / Burn Fraction (into core plasma, shot comsumption = %g kg)',post.z0dinput.machine,post.z0dinput.shot, max(kiloT)));
    set(gca,'yscale','log');
    if isfinite(averaged_burn_fraction) && (averaged_burn_fraction > 0)
	set(gca,'ylim',[0 ceil(averaged_burn_fraction .* 100) .* 3]);
    end
    z0loglin(gca);

    subplot(k,1,2)     
    plot(post.z0dinput.cons.temps,source_T_pellet,post.z0dinput.cons.temps,source_T_gaz,post.z0dinput.cons.temps,source_T_nbi);
    xlabel('time (s)');
    ylabel('Tritium input sources (take into account efficiency, Pa m^3 / s)');
    legend('Pellet','Gaz puff','NBI');
    z0loglin(gca);
     
    subplot(k,1,3)     
    plot(post.z0dinput.cons.temps,source_T_plasma,post.z0dinput.cons.temps,source_leakage,post.z0dinput.cons.temps,source_recycling,post.z0dinput.cons.temps,divertor_flux);
    xlabel('time (s)');
    ylabel('Tritium plasma bilan (inside core plasma, Pa m^3 / s)');
    legend('Input','Losses (plasma to divertor)','Recycling (direct from SOL)','divertor flux');
    z0loglin(gca);

    joint_axes(h,k);

end
