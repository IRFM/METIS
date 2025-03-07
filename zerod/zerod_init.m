% cette fonction prpeare les donnees pour metis a partir de differente source de donnees
function z0dinput = zerod_init(mode_exp,shot,gaz,temps)


% mode :
%  0 -> cronos
%  1 -> TS
%  11 -> TS pre shot
%  2 -> JET
%  -1 -> pas de donnees  (simulation vide)
%  -2 -> pour metis-evolution


if nargin < 1
	mode_exp = 0;
end
if nargin < 2
   shot = [];
end
if nargin < 3
   gaz = [];
end
if nargin < 4
   temps = [];
end

% suppression des donnees experimentales de remplacement
switch mode_exp
case {-1,-2}
	% juste un jeux de donnees vide; pas d'incoherence possible avec les nouvelles: l'effacement n'est pas utile.
otherwise
	rm_cs4m;
end

% parametre et info des parametres 
info = zerod;
z0dinput.option        = info.valeur;
z0dinput.info          = info.info;
% langue
langue                 =  'anglais';
z0dinput.langue        =  'anglais';
% variable de sorties
z0dinput.zsinfo        = zero1t;
z0dinput.profinfo      = z0dprofinfo;
z0dinput.mode_exp      = mode_exp;


% selon le mode
switch mode_exp
case 0
    % acces au donnees de CRONOS
    z0dinput = zerod_init_cronos(mode_exp,shot,gaz,temps,z0dinput);

case 1
    % acces au donnees de Tore Supra
    z0dinput = zerod_init_ts(mode_exp,shot,gaz,temps,z0dinput);

case 11
    % acces au donnees de Tore Supra avant choc (preparation)
    z0dinput = zerod_init_ts_preshot(mode_exp,shot,gaz,temps,z0dinput);

case 12
    % acces au donnees de Tore Supra
    z0dinput = zerod_init_west(mode_exp,shot,gaz,temps,z0dinput);

case 13
    % acces au donnees de Tore Supra avant choc (preparation)
    z0dinput = zerod_init_west_preshot(mode_exp,shot,gaz,temps,z0dinput);

case 2
    % acces au donnees de JET via MSD+
    z0dinput = zerod_init_jet(mode_exp,shot,gaz,temps,z0dinput);

case 3
    % acces au donnees de DIIID via MSD+
    z0dinput = zerod_init_diiid(mode_exp,shot,gaz,temps,z0dinput);

case 4
    % acces au donnees de ASDEX-U
    z0dinput = zerod_init_asdex_u(mode_exp,shot,gaz,temps,z0dinput);

case 5
    % acces au donnees de COMPASS
    z0dinput = zerod_init_compass(mode_exp,shot,gaz,temps,z0dinput);
    
case 6
    % acces au donnees de EAST
    z0dinput = zerod_init_east(mode_exp,shot,gaz,temps,z0dinput);

case 7
    % acces au donnees de EAST
    z0dinput = zerod_init_tcv(mode_exp,shot,gaz,temps,z0dinput);

case 117 
    % lecture des donnees dans l'UAL
    z0dinput = zerod_init_ual(mode_exp,shot,gaz,temps,z0dinput);

case 118 
    % Read UAL data for IMAS
    z0dinput = zerod_init_imas(mode_exp,shot,gaz,temps,z0dinput);

case 119 
    % Read data from code system results.
    z0dinput = zerod_init_code_system(mode_exp,shot,gaz,temps,z0dinput);
    
case 201 
    % Read data from code system results.
    z0dinput = zerod_init_st40(mode_exp,shot,gaz,temps,z0dinput);

case -1
    % creation d'un jeux de donnees vide pour machine quelconque
    z0dinput = zerod_init_empty(mode_exp,shot,gaz,temps,z0dinput);

case -3
    % creation d'un jeux de donnees vide pour machine quelconque
    z0dinput = zerod_init_empty(mode_exp,shot,gaz,temps,z0dinput,'nodialog');

case -2 
    % creation d'un jeux de donnees pour metis_evolution
    z0dinput = zerod_init_evolution(mode_exp,shot,gaz,temps,z0dinput);

otherwise
	error('ZEROD_INIT : unknow mode_exp');
end
if isempty(z0dinput)
	disp('intialisation cancelled');
	z0dinput=[];
	return
end

% mise a jour de la structure experimentale vide
noms = fieldnames(z0dinput.zsinfo);
if ~isfield(z0dinput,'exp0d')
    z0dinput.exp0d=[];
end
exp0d  = z0dinput.exp0d;
if isfield(exp0d,'temps')
	texp = exp0d.temps;
else
	texp = z0dinput.cons.temps;
	exp0d.temps = texp;
end
nbt  = length(texp);
%
vtnan = NaN .* ones(nbt,1);

for k = 1:length(noms)
	nomc = noms{k};
	if isfield(exp0d,nomc)
		var = getfield(exp0d,nomc);
		if length(var) ~= nbt
			disp('dimension mismatch')
			var = mean(var(isfinite(var))) .* ones(nbt,1);
	  		exp0d = setfield(exp0d,nomc,var);
		else
			% si donnnees non valides
			fnan = imag(var);
			var  = real(var);
			ind  = find(fnan~=0 & var == 0);
			if  ~isempty(ind)
				var(ind) = NaN;
			end 
	  		exp0d  = setfield(exp0d,nomc,var);
		end
	else
	  	exp0d = setfield(exp0d,nomc,vtnan);
	end
end

% donnees experimentale
z0dinput.exp0d = exp0d;


if ~isfield(z0dinput.cons,'xece')
      z0dinput.cons.xece = zeros(size(z0dinput.cons.temps));
end

% mise a jour de cons.iso
if ~isfield(z0dinput.cons,'iso')
   z0dinput.cons.iso = zeros(size(z0dinput.cons.temps)); 
elseif length(z0dinput.cons.iso) == 1
   z0dinput.cons.iso = z0dinput.cons.iso .* ones(size(z0dinput.cons.temps)); 
end
% consigne d'injection de tritium par nbi (fraction de la puissance)
%z0dinput.cons.ftnbi = min(1,z0dinput.cons.iso .* 0.5);
if ~isfield(z0dinput.cons,'ftnbi')
  z0dinput.cons.ftnbi  = 0 .* z0dinput.cons.iso;
end

% securite mise en forme et NaN
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
   nomc = noms{k};
   val = getfield(z0dinput.cons,nomc);
   val(~isfinite(val)) = 0;
   z0dinput.cons = setfield(z0dinput.cons,nomc,val(:));
end
noms = fieldnames(z0dinput.geo);
for k=1:length(noms)
   nomc = noms{k};
   val = getfield(z0dinput.geo,nomc);
   val(~isfinite(val)) = 0;
   z0dinput.geo = setfield(z0dinput.geo,nomc,val(:));
end

% securite sur le zeff
% It is OK for option 5 and 11
if z0dinput.option.gaz == 4
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(2.2,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
else
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(1.1,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
end


% gestion auto du ripple
if strcmp(z0dinput.machine,'TS')
   z0dinput.option.rip = 1;
else
   z0dinput.option.rip = 0;

end

% parametre de nom de la machine
z0dinput.option.machine = z0dinput.machine;
z0dinput.option.shot = z0dinput.shot;


%securite largeur LH
if ~isfinite(z0dinput.option.dlh)
	z0dinput.option.dlh = 0.2;
	z0dinput.option.xlh = 0.2;
end

if ~isfinite(z0dinput.option.npar0)
	z0dinput.option.npar = 2;
end

% securite sur le sens de ECRH
% par dedaut perpendiculaire
if ~isfinite(z0dinput.option.sens)
	z0dinput.option.sens = 0;
end


% securite geo
z0dinput.geo.a = max(z0dinput.geo.a,1e-2);
z0dinput.geo.R = max(z0dinput.geo.R,3e-2);
z0dinput.geo.K = max(z0dinput.geo.K,0.1);
z0dinput.geo.b0 = max(z0dinput.geo.b0,1e-4);


% securite sur einj
z0dinput.option.einj  = max(1,z0dinput.option.einj);
if isfield(z0dinput.option,'einj2')
    z0dinput.option.einj2 = max(1,z0dinput.option.einj2);
end

% add possible missing field on profile description dependning on function
% creating the simulation
if ~isfield(z0dinput,'profinfo')
    z0dinput.profinfo      = z0dprofinfo;
end

% transfert dans le workspace
if nargout == 0
   zassignin('base','z0dinput',z0dinput);
end
disp('==> Data ready !');

txt = 'Metis : Fast tokamak simulator';
if isfield(z0dinput,'run')
    txt = sprintf('%s (%s@%d for run = %d and occ = %s)',txt,z0dinput.machine,z0dinput.shot,z0dinput.run,z0dinput.occ);
else
    txt = sprintf('%s (%s@%d)',txt,z0dinput.machine,z0dinput.shot);
end
setappdata(0,'METIS_INTERFACE_TITLE',txt);
if isappdata(0,'METIS_FILENAME');
    rmappdata(0,'METIS_FILENAME');
end

[hfig,h] = zuiformhandle('zeroda');
if isempty(hfig)
  return
end

set(hfig,'name',txt);


function yy = zechan(x,y,xx,mode)

if isempty(x) | isempty(y)
	yy = sqrt(-1) .* ones(size(xx));
else
	yy =interp10d(x,y(:,end),xx,mode);
end


function yi = interp10d(x,y,xi,methode)

if (size(x,1)>1) & (size(x,2)==1)
	indnok = 0;
	nb     = 100;
	while (~isempty(indnok)) & (nb >0)
		indnok = find(diff(x)<=0);
		if ~isempty(indnok)
			x(indnok) = [];
			y(indnok,:) = [];
	
		end
		nb = nb - 1;
	end
end
yi = interp1(x,y,xi,methode);



% fonction d'interpolation pour consigene TS
function  yi = interp1c(x,y,xi)

	[x,ind] = sort(x);
	y       = y(ind);
	while (any(diff(x) <=0))
		indbad = find(diff(x) <=0);
		x(indbad+1) = x(indbad+1) + 1e-6;
	end
	yi = interp1(x,y,xi,'linear');


% calcul de de xece en fonction du R resonance
% c'est un calcul rapide approche
function [xece,nh,apol] = r2xece(r0,z0,d0,a,rant,zant,phi_pol,phi_tor,b0,freq_ghz)

fce = 28 .* b0;
rmin = r0 - a;
rmax = r0 + a;
rc1  = 1 .* r0 .* fce ./ freq_ghz;   
rc2  = 2 .* r0 .* fce ./ freq_ghz;   
rc3  = 3 .* r0 .* fce ./ freq_ghz;   
% harmonique
inh1 = ((rc1 >= rmin) .*  (rc1 <= rmax)); 
inh2 = ((rc2 >= rmin) .*  (rc2 <= rmax)); 
inh3 = ((rc3 >= rmin) .*  (rc3 <= rmax)); 
nh   =  inh1 + (1 - inh1) .* 2 .* inh2 + (1 - inh1) .* (1 - inh2) .* 3 .*  inh3;
% position du maximum
rc   = nh .* r0 .* fce ./ freq_ghz;
rs   = rant ./ sqrt(2) ./ max(eps,abs(sin(phi_tor))) .* sqrt( 1 - sqrt(1 - 4 .* rc .^ 2 ./ rant .^ 2 .* sin(phi_tor) .^ 2));
rs(abs(sin(phi_tor)) < 0.01) = rc(abs(sin(phi_tor)) < 0.01); 
rabs = medfilt1((rc + rs) ./ 2,3);
% 1ere estimation
rho = sqrt((rabs - r0) .^ 2 +  (zant + (rant - rabs) .* tan(phi_pol) - z0) .^ 2);
%convergence
for k=1:10
	d    = d0 .* (1 - (rho ./a) .^ 2);
	rhos = sqrt((rabs - (r0 + d)) .^ 2 +  (zant + (rant - rabs) .* tan(phi_pol) - z0) .^ 2);
	rho  = 0.6 .* rho + 0.4 .* rhos; 
end
xece = medfilt1(rho ./ a,3);
apol = medfilt1(acos(min(1,max(-1,(rabs - (r0 + d)) ./ rho))) / pi .* 180,3);


% etalonnage des consigne sure le choc de reference
function cons = etalon(mesure_exp,cons_exp,cons)

% cas 0
if all(cons == 0)
	fprintf('1 * consigne + 0   (consigne nulle)\n');
	return
end

% 1 recherche des point valides
delta = abs(mesure_exp - cons_exp) ./ max(cons);
indok = find((delta <= 0.5) & (cons_exp > eps) & (mesure_exp > eps));
if length(indok) < 7
	fprintf('1 * consigne + 0 (donnees insuffisantes)\n');
	return
end

% etalonnage lineaire
mm = mean(cons_exp(indok));
ss = std(cons_exp(indok));
if ss ==0
	ss = 1;
end
xx = (cons_exp(indok) - mm) ./ ss;
warning off
pp = polyfit(xx,mesure_exp(indok),1);
warning on
if any(~isfinite(pp))
	fprintf('1 * consigne + 0 (donnees insuffisantes)\n');
	return
end
if pp(1) < 0
	fprintf('1 * consigne + 0 (pente negative)\n');
	return
end
% nouvelle consigne corrigee
xx = (cons - mm - min(cons)) ./ ss;
cons = max(0,polyval(pp,xx) .* (cons > 0)) + min(cons);

fprintf('%g * consigne + %g\n',pp(1) ./ ss , pp(2) - pp(1) .* mm ./ ss)



%============================================================
%============================================================
% debut de la section pour les differentes source de donnees
%============================================================
%============================================================

%============================================================
% acces au donnees de CRONOS
%============================================================
function z0dinput = zerod_init_cronos(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
% recherche des infos
z0dinput.cons.temps    = evalin('base','data.gene.temps');
% choix de la composition
compo = evalin('base','param.compo');
if compo.z(1) == 1
  if compo.a(1) == 1
      z0dinput.option.gaz  = 1;
  elseif compo.a(1) == 2;
      if any( (compo.a == 3) & (compo.z == 1))
	z0dinput.option.gaz  = 3;
      else
	z0dinput.option.gaz  = 2;
      end
  else
      z0dinput.option.gaz  = 3;
  end
else
    z0dinput.option.gaz  = 4;
end  

% zmax
z0dinput.option.zmax = max(compo.z);
z0dinput.option.zimp = max(compo.z(z0dinput.option.zmax~=compo.z));
z0dinput.option.rimp = 0.1;
z0dinput.option.zeff = 0;

% coefficient de reflexion des parois
cimp = evalin('base','param.cons.impur');
if isfield(cimp,'rimp')
  z0dinput.option.rimp = cimp.rimp;
end
if isfield(cimp,'reflex')
  z0dinput.option.rw = cimp.reflex;
end
% rapport isotopique
if isfield(cimp,'cmin1')
  z0dinput.option.cmin = cimp.cmin1;
end


if z0dinput.option.gaz  == 3
  nhnd = evalin('base','data.cons.nhnd');
  nTnD = imag(nhnd);
  nTnDm = mean(nTnD(isfinite(nTnD)));
  if isempty(nTnD) | (nTnDm == 0)
      if isfield(cimp,'cmin1') & isfield(cimp,'cmin2')
	    indT = find((compo.a == 3) & (compo.z == 1));
	    indD = find((compo.a == 2) & (compo.z == 1));
	    if isempty(indT)
		    z0dinput.cons.iso = 0
	    elseif indT == 1
		    if indD == 2
			    z0dinput.cons.iso = 1./ cimp.cmin1;
		    else
			    z0dinput.cons.iso = 1./ cimp.cmin2;
		    end
	    elseif indT == 2
		    z0dinput.cons.iso = cimp.cmin1;
	    elseif indT == 3
		    z0dinput.cons.iso = cimp.cmin2;
	    elseif (indT == 4) | (indT == 5)
		    impur = evalin('base','data.impur.impur');
		    nDnT  = squeeze(impur(:,1,indT)) ./ squeeze(impur(:,1,indD));
		    nDnT(~isfinite(nDnT)) = mean(nDnT(isfinite(nDnT)));
		    if ~isempty(nDnT)
			    z0dinput.cons.iso = nDnT;
		    else  
			    disp('rapport nTnD indisponible');
			    z0dinput.cons.iso = 0.01;
		    end
	    end
	else 
	    z0dinput.cons.iso =1;
	end
  else
      z0dinput.cons.iso = nTnD;
  end
else 
	    z0dinput.cons.iso = 0;
end
% mode FCI
cfci = evalin('base','param.cons.fci');
if isfield(cfci,'mode')
  if iscell(cfci.mode)
	modefci = cfci.mode{1};
  else
	modefci =  cfci.mode;
  end
  switch   deblank(modefci)
  case 'HMIN_2T'
      z0dinput.option.fwcd  = 0;
      z0dinput.option.mino  = 'T';
  case 'HMIN_H'
      z0dinput.option.fwcd  = 0;
      z0dinput.option.mino  = 'H';
  case 'HMIN_He3'
      z0dinput.option.fwcd  = 0;
      z0dinput.option.mino  = 'He3';
  case 'HMIN_He'
      z0dinput.option.fwcd  = 0;
      z0dinput.option.mino  = 'He4';
  case 'FWCD'
      z0dinput.option.fwcd = 1;
  case 'FWEH'
      z0dinput.option.fwcd = 2;
  otherwise
      z0dinput.option.fwcd = 0;
  end
else
  z0dinput.option.fwcd = 0;
end
if isfield(cfci,'frequence')
  z0dinput.option.freq = max(1,cfci.frequence(1));
else
  z0dinput.option.freq = 57;
end
% mode d'utilisation du zerod
modev = evalin('base','data.mode.cons.psi');
if mean(modev(isfinite(modev)))  > 0.5
  z0dinput.option.vloop = 1;
else
  z0dinput.option.vloop = 0;
end
z0dinput.option.ohm = 1;

% reglage de fce
cfce = evalin('base','param.cons.fce');
pfce = abs(evalin('base','data.cons.fce'));
if isfield(cfce,'synergie')
  z0dinput.option.synergie = cfce.synergie(1);
else
  z0dinput.option.synergie = 0;
end
if isfield(cfce,'angle_tor')
  sens = sum(pfce .* (ones(size(pfce,1),1) *sin(cfce.angle_tor./180.*pi)),2); 
  if any(isfinite(sens))
    sens = mean(sens(isfinite(sens)));
  else
    sens = 0;
  end
  if ~isempty(sens)
      if isfinite(sens)
	    z0dinput.option.sens = sign(sens);
      else
	    z0dinput.option.sens = 0;	 
      end
  else
      z0dinput.option.sens = 0;
  end
else
  z0dinput.option.sens = 0;
end
pfce = abs(evalin('base','data.source.fce.el'));
indok = find(all(isfinite(pfce),2));
if isempty(indok)
  z0dinput.cons.xece = zeros(size(z0dinput.cons.temps));
else
  indon = find(any(pfce(indok,:) >0,2));
  if isempty(indon)
      z0dinput.cons.xece = zeros(size(z0dinput.cons.temps));
  else
      pfce(~isfinite(pfce)) = 0;
      mask = (pfce == (max(pfce,[],2)*ones(1,size(pfce,2))));
      xece = sum(mask .* (ones(size(pfce,1),1) * evalin('base','param.gene.x')),2) ./ max(1,sum(mask,2));
      z0dinput.cons.xece = max(0,min(1,xece));
  end
end

% angle ece ?
phys   = evalin('base','param.phys');
geoin  =  evalin('base','data.geo');
if isfield(cfce,'modpolar') & isfield(cfce,'freq_ghz')
    if cfce.modpolar(1) == 2
	    bk       = 2 * pi * phys.me ./ 1 ./phys.e .* cfce.freq_ghz(1) .* 1e9;
    else 
	    bk       = 2 * pi * phys.me ./ 2 ./phys.e .* cfce.freq_ghz(1) .* 1e9;
    end
    b0  = mean(geoin.b0(isfinite(geoin.b0)));
    r0  = mean(geoin.b0(isfinite(geoin.r0)));
    a   = mean(geoin.b0(isfinite(geoin.a)));
    rk  = bk ./ b0 .* r0;
    if (rk - r0) < -(a/3)
	z0dinput.option.angle_ece  = 180;
    elseif (rk - r0) > (a/3)
	z0dinput.option.angle_ece  = 0;
    else
	z0dinput.option.angle_ece  = 90;
    end
else
  z0dinput.option.angle_ece  = 90;
end


% injection de neutres
machine = evalin('base','param.from.machine');
switch machine
case 'TS'
  z0dinput.option.angle_nbi = 90;
  z0dinput.option.modeh     = 0;
  z0dinput.option.rip        = 1;
  z0dinput.option.etalh     = 0.8;
  z0dinput.option.lhmode    = 3;
  z0dinput.option.xlh        = 0.2;
  z0dinput.option.dlh        = 0.3;
  z0dinput.option.wlh        = 32 .* 11e-3; 
  z0dinput.option.rtang       = 1.47;
  z0dinput.option.nphi = 30;
case 'JET'
  z0dinput.option.angle_nbi = 60;
  z0dinput.option.modeh     = 1;
  z0dinput.option.etalh     = 0.8;
  z0dinput.option.lhmode    = 3;
  z0dinput.option.xlh        = 0.4;
  z0dinput.option.dlh        = 0.3;
  z0dinput.option.nphi     =  15;
  z0dinput.option.wlh        = 32 .* 11e-3; 
case 'ITER'
  z0dinput.option.angle_nbi = 75;  
  z0dinput.option.modeh     = 1;
  z0dinput.option.etalh     = 0.8;
  z0dinput.option.lhmode    = 3;
  z0dinput.option.xlh        = 0.7;
  z0dinput.option.dlh        = 0.4;
  z0dinput.option.nphi     =  12;
  z0dinput.option.wlh        = 24 .* 22e-3 + 11e-3; % a verifier !!!

case 'ASDEX'
  z0dinput.option.angle_nbi = 60;
  z0dinput.option.modeh     = 1;
otherwise
  z0dinput.option.angle_nbi = 0;
  z0dinput.option.modeh     = 1;
end   

z0dinput.option.qdds      = 0.95;
z0dinput.option.kidds     = 3;

% lecture de etalh
z0dinput.cons.temps    = evalin('base','data.gene.temps');
z0dinput.exp0d       = [];
hybcons = evalin('base','param.cons.hyb');
xdur = evalin('base','data.prof.xdur');
if isfield(hybcons,'efficacite')
  etalh = hybcons.efficacite;
  if isfinite(etalh(1))
      z0dinput.option.etalh     =  etalh(1) ;
      z0dinput.option.lhmode    =  2;
  end
end
if all(isfinite(xdur(:))) & ~all(xdur(:) == 0)
    if isfield(hybcons,'xdur') 
	    if  hybcons.xdur == 1
		    g3s  = sum(xdur,1);
		    d3   = linspace(0,1,length(g3s));
		    %xlh  = trapz(d3,(1-d3) .* d3.* g3s,2) ./trapz(d3,(1-d3) .* g3s,2);
		    xlh  = d3(max(find(g3s ==max(g3s))));
		    dlh  = sqrt(trapz(d3,d3 .^ 2 .* g3s,2) ./trapz(d3,g3s,2) - xlh.^ 2);
		    z0dinput.option.xlh = xlh;
		    z0dinput.option.dlh = dlh;
		    z0dinput.option.wlh  = 0;
	    else
		    z0dinput.option.xlh =  hybcons.centre(1);
		    z0dinput.option.dlh =  hybcons.largeur(1); % attention pb de definition
	    end
    end
		    
    z0dinput.exp0d.XDURt = z0dinput.cons.temps;
    z0dinput.exp0d.XDURx = evalin('base','param.gene.x');
    z0dinput.exp0d.XDURv = xdur;		

elseif isfield(hybcons,'xdur')
  z0dinput.option.xlh =  hybcons.centre(1);
  z0dinput.option.dlh =  hybcons.largeur(1);
end

consin                 = evalin('base','data.cons');
z0dinput.cons.ip       = consin.ip;
z0dinput.cons.flux     = consin.flux;
ne                     = evalin('base','data.prof.ne');
x                      = linspace(0,1,size(ne,2));
z0dinput.cons.nbar     = trapz(x,ne,2);
z0dinput.cons.picrh    = sum(abs(consin.fci),2);
z0dinput.cons.plh      = sum(abs(consin.hyb),2);
try
    warning off
    npar0                  = mean(trapz(z0dinput.cons.temps,sum(abs(consin.hyb) .* angle(consin.hyb),2),1) ./ ...
			      max(1,trapz(z0dinput.cons.temps,sum(abs(consin.hyb),2),1)));
    z0dinput.option.npar0  =npar0;
    warning on
catch
    z0dinput.option.npar0  = 2;
end
if ~isfinite(z0dinput.option.npar0)
	    z0dinput.option.npar0 = 2;
end
z0dinput.cons.pnbi     = sum(real(consin.idn),2);
cidn = evalin('base','param.cons.idn');
if isfield(cidn,'energie')
  energie = cidn.energie;
  indok   = find(isfinite(energie) & all(isfinite(consin.idn),1));
  energie = ones(size(consin.idn,1),1) * energie;
  einj    = trapz(z0dinput.cons.temps,sum(energie(:,indok) .* consin.idn(:,indok),2),1) ./  ...
	    max(1,trapz(z0dinput.cons.temps,sum(consin.idn(:,indok),2),1));
  if isfinite(einj)
	z0dinput.option.einj = einj;
  else
	z0dinput.option.einj = max(cidn.energie(isfinite(cidn.energie)));
  end
  if isempty(z0dinput.option.einj)
      z0dinput.option.einj = 1e5; % valeur de repli
  end
  if z0dinput.option.einj < 1e4
	z0dinput.option.einj = 1e4;
  end
else
  z0dinput.option.einj = 1e5; % valeur de repli
end
if isfield(cidn,'charge') & isfield(cidn,'masse')
    indok   = find(isfinite(cidn.charge) & isfinite(cidn.masse) & all(isfinite(consin.idn),1));
    charge  = cidn.charge(indok);
    masse   = cidn.masse(indok);
    pvnbi   = consin.idn(:,indok);
    pnbi    = max(1,sum(real(pvnbi),2));
    vt      = ones(size(z0dinput.cons.temps));
    switch z0dinput.option.gaz
    case 2
	    ind1 = find((charge == 1) & (masse == 1));
	    ind2 = find((charge == 1) & (masse == 2));
	    if ~isempty(ind1) & ~isempty(ind2);
		p1 =  sum(pvnbi(:,ind1),2) ./ pnbi;
		p2 =  sum(pvnbi(:,ind2),2) ./ pnbi;
		z0dinput.cons.ftnbi = p1 ./ max(1,p1 + p2);
	    elseif isempty(ind1)  & ~isempty(ind2)
		z0dinput.cons.ftnbi = 0 .* vt;
	    elseif ~isempty(ind1) & isempty(ind2)
		z0dinput.cons.ftnbi = vt;
	    else
		z0dinput.cons.ftnbi = 0 .* vt;
	    end
    case 3
	    ind3 = find((charge == 1) & (masse == 3));
	    ind2 = find((charge == 1) & (masse == 2));
	    if ~isempty(ind3) & ~isempty(ind2);
		p3 =  sum(pvnbi(:,ind3),2) ./ pnbi;
		p2 =  sum(pvnbi(:,ind2),2) ./ pnbi;
		z0dinput.cons.ftnbi = p3 ./ max(1,p3 + p2);
	    elseif isempty(ind3)  & ~isempty(ind2)
		z0dinput.cons.ftnbi = 0 .* vt;
	    elseif ~isempty(ind3) & isempty(ind2)
		z0dinput.cons.ftnbi = vt;
	    else
		z0dinput.cons.ftnbi = 0 .* vt;
	    end
    otherwise
	    z0dinput.cons.ftnbi = 0 .* z0dinput.cons.temps;
    end
else
	  z0dinput.cons.ftnbi = 0 .* z0dinput.cons.temps;
end

z0dinput.cons.pecrh    = sum(abs(consin.fce),2);
z0dinput.cons.hmore    = ones(size(z0dinput.cons.temps));
z0dinput.cons.zeff     = evalin('base','data.gene.zeffm');
z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = consin.zeffm(~isfinite(z0dinput.cons.zeff));
z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = 3;

li                     = evalin('base','data.gene.li');
librute                = li;
li(~isfinite(li))     = 1;
z0dinput.option.li    = li(1);
% 1- transformation des structure Cronos -> 0D
z0dinput.geo.a     = geoin.a;
z0dinput.geo.R     = geoin.r0;
z0dinput.geo.K     = geoin.e1;
z0dinput.geo.d     = (abs(geoin.trh1) + abs(geoin.trb1)) ./ 2;
z0dinput.geo.b0    = geoin.b0;
z0dinput.geo.z0    = geoin.z0;
z0dinput.geo.vp    = evalin('base','data.gene.volume');
z0dinput.geo.sp    = evalin('base','data.gene.surface');
z0dinput.geo.sext  = evalin('base','data.equi.vpr(:,end)');

% si la separatrice est utilise
if mean((geoin.mode == 2)) > 0.5
    sepa.R = double(geoin.R);
    sepa.Z = double(geoin.Z);
    % calcul des moments
    % la courbe doit etre fermee
    if (sepa.R(1,1) ~= sepa.R(1,end)) | (sepa.Z(1,1) ~= sepa.Z(1,end))
	    sepa.R(:,end+1) = sepa.R(:,1);
	    sepa.Z(:,end+1) = sepa.Z(:,1);
    end

    % calcul des moments
    % centre pour angle d'integration
    rc = mean(sepa.R,2);
    zc = mean(sepa.Z,2);
    vc = ones(1,size(sepa.R,2));
    uc = unwrap(angle((sepa.R-rc*vc) + sqrt(-1) .* (sepa.Z  -zc*vc)));
    uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
    uc(:,1)   = uc(:,end) + 2 .* pi;
    xu    = linspace(0,1,length(vc));
    %dudx  = pdederive(xu,uc,2,2,2,1);
    %dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
    %dudx(:,end) = dudx(:,1);
    dRdx  = pdederive(xu,sepa.R,2,2,2,1);
    dZdx  = pdederive(xu,sepa.Z,2,2,2,1);
    % calcul de R0 et Z0
    maskrmax  = (sepa.R == (max(sepa.R,[],2) * vc));
    geo.z0        = sum(sepa.Z .* maskrmax,2) ./ sum(maskrmax,2);
    % recalcul des parametres sur le vecteur final
    rmin  = min(sepa.R,[],2);
    rmax  = max(sepa.R,[],2);
    geo.a = 0.5 .* (rmax - rmin);
    geo.R = 0.5 .* (rmax + rmin);
    zmin  = min(sepa.Z,[],2);
    zmax  = max(sepa.Z,[],2);
    %geo.z0      =(zmax + zmin) ./ 2;
    geo.K    = abs(trapz(xu,sepa.Z .*  dRdx,2) ./ pi ./ geo.a .^ 2);
    %geo.K  = (zmax -zmin) ./ 2 ./ geo.a;
    rzmax = geo.R;
    rzmin = geo.R;
    for k = 1:size(sepa.Z,1)
	    rzmax(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmax(k))));
	    rzmin(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmin(k))));
    end
    uu   =  angle(rzmax - geo.R + sqrt(-1) .* (zmax - geo.z0));
    ul   =  angle(rzmin - geo.R + sqrt(-1) .* (zmin - geo.z0));
    tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
    tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
    tm   =  (tl + tu) ./ 2;
    geo.d = sin(tm);


    z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
    z0dinput.geo.z0      = geo.z0;     % centre geometricque du plasma en Z (m)
    z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
    z0dinput.geo.K       = geo.K;     % elongation (b/a)
    z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
    z0dinput.geo.b0      = geoin.b0; % champ toroidal a R = 6.2 m (T)


    z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
    z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)




    t  = asin(max(0,min(1,geo.d)));
    u  = linspace(0,2.*pi,201);
    vu = ones(size(u));
    Rtest  = geo.R *vu + (geo.a * vu) .* cos(ones(size(geo.R,1),1) * u + t * sin(u));
    Ztest = (geo.a .* geo.K) * sin(u);

    h = findobj(0,'type','figure','tag','z0geosepa');
    if isempty(h)
    h=figure('tag','z0geosepa');
    else
    figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,1:length(geo.R),Rtest,Ztest+geo.z0 * vu,'color','r','marker','o','linestyle','none');
    zplotprof(gca,1:length(geo.R),sepa.R,sepa.Z,'color','b','marker','none','linestyle','-');
    axis('square')
    axis('equal')

end

% pour le fit de l'energie
z0dinput.exp0d.w = evalin('base','data.gene.wdia');

z0dinput.machine   = evalin('base','param.from.machine');
z0dinput.shot      = evalin('base','param.from.shot.num');

z0dinput.option.signe = evalin('base','param.gene.signe.ip .* param.gene.signe.b0');


%============================================================
% acces au donnees de Tore Supra
%============================================================
function z0dinput = zerod_init_ts(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
% cas donnees TS
z0dinput.exp0d         = [];
% 1 - numero du choc
if isempty(gaz) | isempty(shot)

    try 
	numshot = fix(evalin('base','param.from.shot.num'));
    catch
	  numshot = 32299;
    end 

    switch langue 
    case 'francais'
	    prompt={'Numero du choc :','Charge du gaz principal :'};
	    def={sprintf('%d',numshot),'1'};
	    dlgTitle='Lecture des donnees de Tore Supra';
    otherwise
	    prompt={'shot number :','charge of main gas:' };
	    def={sprintf('%d',numshot),'1'};
	    dlgTitle='Access to Tore Supra data';
    end
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
            z0dinput = [];
	    return
    end
    shot  = str2num(answer{1});
    gaz  = str2num(answer{2});

end

% lecture des donnees
[ip,tip] = tsbase(shot,'sipmes');
if isempty(ip)
	disp('No plasma')
        z0dinput = [];
	return
end
if gaz < 0
  tt = tip(tip >=0);
else
  tt = tip(tip >=0 & ip > 3e-2);
end
if isempty(tt)
	disp('No plasma current');
        z0dinput = [];
    	return
end

if ~isempty(temps)
	% on prend le vecteur donnee en entree
elseif gaz < 0
    gaz = abs(gaz);
    temps = tt;
    temps(find(diff(temps)<=0)) =[];
else
    dt    = max(mean(diff(tt)),0.1);
    temps = (min(tt):dt:max(tt))';

    if length(temps) > 1001
	temps = linspace(min(tt),max(tt),1001)';
    end	
end
%
% probleme tip
%	
indtip = find(diff(tip) <= 0);
if ~isempty(indtip)
    tip(indtip) = [];
    ip(indtip) = [];
end
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = max(1,interp10d(tip,ip,temps,'nearest') .* 1e6);
ip   = z0dinput.cons.ip;

% flux au bord
[fluxbord,tfluxbord] = tsbase(shot,'gfluxnoy%4');
if isempty(fluxbord)
    [fluxbord,tfluxbord] = tsbase(shot,'gfluxnoy%3');
    z0dinput.cons.flux    = interp10d(tfluxbord,fluxbord,temps,'nearest')./ 2 ./ pi;
else
    z0dinput.cons.flux    = interp10d(tfluxbord,fluxbord,temps,'nearest') ./ 2 ./ pi;
end

[gplasma,tplasma]       = tsbase(shot,'gplasma');
if ~isempty(gplasma)
    tplasma                = tplasma(:,1);
    %
    % probleme tplasma
    %	
    indtplasma = find(diff(tplasma) <= 0);
    if ~isempty(indtplasma)
	tplasma(indtplasma) = [];
	gplasma(indtplasma,:) = [];
    end

    z0dinput.geo.a         = interp10d(tplasma,gplasma(:,3),temps,'nearest');
    z0dinput.geo.R         = interp10d(tplasma,gplasma(:,1),temps,'nearest');
    z0dinput.geo.z0         = interp10d(tplasma,gplasma(:,2),temps,'nearest');
    e                      = interp10d(tplasma,gplasma(:,4),temps,'nearest');
    z0dinput.geo.K         = 1 + max(-0.25,min(e,0.25));
    z0dinput.geo.d         = zeros(size(temps));
else
    [R0,tplasma]       = tsbase(shot,'srmaj');
    z0dinput.geo.R     = interp10d(tplasma,R0,temps,'nearest');
    [a,tplasma]       = tsbase(shot,'samin');
    z0dinput.geo.a     = interp10d(tplasma,a,temps,'nearest');
    [z0,tplasma]       = tsbase(shot,'szpos');
    z0dinput.geo.z0     = interp10d(tplasma,z0,temps,'nearest');
    [e1,te1,ce1]        = tsbase(shot,'sellip');
    if ~isempty(e1)
	  z0dinput.geo.K     = interp10d(te1,e1,temps,'nearest');
    else
	  z0dinput.geo.K     = ones(size(temps));
    end
    z0dinput.geo.d    = zeros(size(temps));
end

% separatrice
[grho,tgrho]   = tsbase(shot,'grho');
if ~isempty(grho)
	tgrho          = tgrho(:,1);
	indnok   = find(all(grho < 0,2));
	grho(indnok,:) = [];
	tgrho(indnok) = [];
	indnok   = find(any(grho < 0.1,2));
	grho(indnok,:) = [];
	tgrho(indnok) = [];
end
if ~isempty(grho)

	grho     = interp1(tgrho,grho,temps,'nearest','extrap');
	% attente d'acces base
	zaxe           = 0;
	raxe           = 2.42;

	alpha          = (0:15:345) ./ 180 .* pi;
	vt             = ones(size(grho,1),1);
	rr             = raxe + grho .* cos(vt*alpha);
	zz             = zaxe + grho .* sin(vt*alpha);
	rr(:,end+1)    = rr(:,1);
	zz(:,end+1)    = zz(:,1);
	alpha(end+1)   = 2*pi;
	teta           = linspace(0,2*pi,201);
	Rext           = tsplinet(vt*alpha,rr,vt*teta);
	Zext           = tsplinet(vt*alpha,zz,vt*teta);
	Rext(:,end)    = Rext(:,1);
	Zext(:,end)    = Zext(:,1);

	% recalcul des parametres pour verification
	ve    = ones(1,size(Rext,2));
	rmin  = min(Rext,[],2);
	rmax  = max(Rext,[],2);
	ra    = max(1,0.5 .* (rmin + rmax));
	a     = max(0.01,0.5 .* (rmax - rmin));
	zmin  = min(Zext,[],2);
	zmax  = max(Zext,[],2);
	za    = (zmin + zmax) ./ 2;
	b     = 0.5 .* (zmax - zmin);
	k     = max(0.5,b ./ a);
	mask1 = (Zext == (max(Zext,[],2)*ve));
	mask2 = (Zext == (min(Zext,[],2)*ve));

	rzmax = max(Rext .* mask1,[],2);
	rzmin = max(Rext .* mask2,[],2);
	cl    = ra - rzmin;
	cu    = ra - rzmax;
	d     = (cl+cu) ./2 ./ a;
	z0dinput.exp0d.Rsepa = Rext;
	z0dinput.exp0d.Zsepa = Zext - za * ones(1,size(Zext,2));
	z0dinput.geo.a       = a;
	z0dinput.geo.R       = ra;
	z0dinput.geo.z0      = za;
	z0dinput.geo.K       = k;
	z0dinput.geo.d       = d;

end

[itor,titor]           = tsbase(shot,'sitor');
if ~isempty(itor)
    ind        	       = find(titor >=0);
    itor                   = itor(ind);
    titor                  = titor(ind);
    [b,a]                  = butter(11,0.1);
    itor                   = filtfilt(b,a,itor);
else
	itor=tsmat(shot,'EXP=T=S;GENERAL;ITOR')*ones(size(temps));
	titor=temps;
end
rb0                    = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi;
rb0                    = interp10d(titor,rb0,temps,'nearest');
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
itor                   = interp10d(titor,itor,temps,'nearest');

[piqne,tne]            = tsbase(shot,'spiq');
[nl,tnl]               = tsbase(shot,'gnl');
if isempty(nl)
    [nl,tnl]           = tsbase(shot,'gtrnl');	    
end
[nli,tnli]             = tsbase(shot,'snli');
nli                    = nli .* 1e19;

if isempty(nl)
    disp('No density measurment available')
    z0dinput = [];
    return
elseif all(nl(:) <= 3e17)
    disp('No density measurment available')
    z0dinput = [];
    return	
else
    tnl=tnl(:,1);
    ind = find(diff(tnl) <=0);
    while(~isempty(ind) & ~isempty(tnl))
      tnl(ind) = [];
      nl(ind,:)=[];
      ind = find(diff(tnl) <=0);
    end
end
nl = nl - mean(nl((tnl <0) & (tnl > -2)));
nbar                   = interp10d(tnl,max(nl,[],2),temps,'nearest')./ 2 ./ z0dinput.geo.a;
% securite anti Nan
nbarm                  = mean(nbar(isfinite(nbar)));
nbar(~isfinite(nbar))  = nbarm;
nbar                   = max(1e17,nbar);

if ~isempty(nli)
	% suppresion de l'offset
	nli = nli - mean(nli((tnli <0) & (tnli > -2)));
	nbari          = interp10d(tnli,nli,temps,'nearest')./ 2 ./ z0dinput.geo.a;
	nbarm                  = mean(nbari(isfinite(nbari)));
	nbari(~isfinite(nbari))  = nbarm;
	nbari                   = max(1e17,nbari);

	indfit = find((abs(nbari - nbar)./nbari) < 0.1);
	if length(indfit) >11		
		pp     = polyfit(nbari(indfit)./1e19,nbar(indfit)./1e19,1);
		nbarm   = polyval(pp,nbari./1e19) .* 1e19;
		bb      = min(1,max(0,temps - 1));
		nbarm   = bb .* nbarm + (1- bb) .* nbar; 
		%figure(23);clf;plot(temps,nbar,'b',temps,nbari,'r',temps,nbarm,'go');
		nbar    = nbarm;
	end
end
z0dinput.cons.nbar     = nbar;
	
[nemoy,tne]            = tsbase(shot,'snmoy');
if ~isempty(nemoy)
    nemoy                   = interp10d(tne,nemoy,temps,'nearest');
    piqne                  = max(1.1,min(5,interp10d(tne,piqne,temps,'nearest')));
    ane                    = piqne - 1;
    ne0                    = nemoy .* piqne;
else
    [nemoy,tne]            = tsbase(shot,'strnmoy');
    if ~isempty(nemoy)
	nemoy                   = interp10d(tne,nemoy,temps,'nearest');
	piqne                  = 1.5 .* ones(size(itor));
	ane                    = piqne - 1;
	ne0                    = nbar .*  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi));
    else
	nemoy = nbar;
	piqne = 1.5 .* ones(size(itor));
	ane                    = piqne - 1;
	ne0                    = nbar .*  (2 .* gamma(ane + 1.5) ./ gamma(ane + 1 ) ./ sqrt(pi));
    end
end
	
	
  [gbilan,tbilan]         = tsbase(shot,'gbilan');
  if ~isempty(gbilan)
      tbilan                 = tbilan(:,1);
      %
      % probleme tbilan
      %	
      indtbilan = find(diff(tbilan) <= 0);
      if ~isempty(indtbilan)
	tbilan(indtbilan) = [];
	gbilan(indtbilan,:) = [];
      end
      z0dinput.cons.picrh    =  max(0,interp10d(tbilan,gbilan(:,3),temps,'nearest')) .* 1e6;
      z0dinput.cons.plh      = max(0,interp10d(tbilan,gbilan(:,2),temps,'nearest')) .* 1e6;
      z0dinput.cons.pnbi     = zeros(size(temps)); 
      pecrh                  = zecrh(shot,temps);
      z0dinput.cons.pecrh    = max(0,sum(pecrh,2).*1e6);
      z0dinput.cons.pecrh    = z0dinput.cons.pecrh .* (z0dinput.cons.pecrh >=1e4);
      z0dinput.cons.hmore    = ones(size(temps));
      ptot                   =  max(0,interp10d(tbilan,gbilan(:,5),temps,'nearest')) .* 1e6;
      pohm                   =  interp10d(tbilan,gbilan(:,1),temps,'nearest') .* 1e6;
else
      [pfci,tpfci]=tsbase(shot,'spuiss');
      if isempty(pfci)
	  [pfci,tpfci]=tsbase(shot,'gptra');
	  if isempty(pfci)
		  [pfci,tpfci]=tsbase(shot,'gpfci');
	  end
	  if ~isempty(pfci)
		  tpfci = tpfci(:,1);
		  pfci  = sum(pfci .* (pfci >=0),2);
	  end
      end
      if ~isempty(pfci)
	  z0dinput.cons.picrh    =  max(0,interp10d(tpfci,pfci,temps,'nearest')) .* 1e6;
      else
	  z0dinput.cons.picrh    =  zeros(size(temps));
      end
      [phyb,tphyb]=tsbase(shot,'GPHYB%3');
      if isempty(phyb) 
	  [phyb,tphyb]=tsbase(shot,'GHYB%3');
      end
      if ~isempty(phyb) 
	  z0dinput.cons.plh      = max(0,interp10d(tphyb,phyb,temps,'nearest')) .* 1e6;
      else
	  z0dinput.cons.plh      = zeros(size(temps));
      end
      z0dinput.cons.pnbi     = zeros(size(temps)); 
      z0dinput.cons.pecrh    = zeros(size(temps));
      z0dinput.cons.hmore    = ones(size(temps));
      [pohm,tpohm]=tsbase(shot,'spohm');
      if ~isempty(phyb) 
	  pohm      = max(0,interp10d(tpohm,pohm,temps,'nearest')) .* 1e6;
      else
	  pohm      = zeros(size(temps));
      end
      ptot                   =  pohm + z0dinput.cons.picrh + z0dinput.cons.plh;
end

if shot >= 44000

	[zeffm,tzeffm]         = tsbase(shot,'szfbrmtan');
	if ~isempty(zeffm)
		zeffm          = interp10d(tzeffm,zeffm,temps,'nearest');
	else
		[zeffm,tzeffm]         = tsbase(shot,'szfbrm');
		if ~isempty(zeffm)
			zeffm          = interp10d(tzeffm,zeffm,temps,'nearest');
		else
			zeffm          = 2.3 .* ones(size(temps));
		end
	end
        
elseif shot > 36599
      [zeffm,tzeffm]         = tsbase(shot,'szfbrm');
      if ~isempty(zeffm)
	zeffm          = interp10d(tzeffm,zeffm,temps,'nearest');
	[zeffm1,tzeffm1]         = tsbase(shot,'szfbrmtan');
	if ~isempty(zeffm1)
	    zeffm1          = interp10d(tzeffm1,zeffm1,temps,'nearest');
	    indnofci = find((z0dinput.cons.picrh<= 1e6) & isfinite(zeffm) & isfinite(zeffm1));
	    if length(indnofci) > 5
		  rapz = max(.3,min(3,mean(zeffm1(indnofci) ./ zeffm(indnofci))))
		  zeffm = zeffm .* rapz;
	    end
	end
      end
else
	  [zeffm,tzeffm]         = tsbase(shot,'szfbrm');
	  if ~isempty(zeffm)
		    zeffm          = interp10d(tzeffm,zeffm,temps,'nearest');
	  end
end
if isempty(zeffm)
	z0dinput.cons.zeff    = zeffscaling(z0dinput.cons.nbar./ 1e19,ptot ./ 1e6,ip ./ 1e6, ...
			z0dinput.geo.a,z0dinput.geo.R,itor,gaz);
else
	if gaz == 2
		z0dinput.cons.zeff    = max(2,min(16,zeffm));
	else
		z0dinput.cons.zeff    = max(1,min(7,zeffm));
	end

end                  

d0                    = 0.1 .* ones(size(temps));
[gqbeli,tgqbeli]      = tsbase(shot,'gqbeli');
if ~isempty(gqbeli)
    tgqbeli               = tgqbeli(:,1);
    %
    % probleme tqbeli
    %	
    indtgqbeli = find(diff(tgqbeli) <= 0);
    if ~isempty(indtgqbeli)
      tgqbeli(indtgqbeli) = [];
      gqbeli(indtgqbeli,:) = [];
    end
    li                    = interp10d(tgqbeli,abs(gqbeli(:,3)),temps,'nearest');
    librute               = li;
    li(~isfinite(li))     = 1;
    li                    = min(10,max(0.1,li));
    li                    = medfilt1(li,7);
    z0dinput.option.li    = li(7);
    
    % d0 util pour la suite
    d0                    = interp10d(tgqbeli,abs(gqbeli(:,5)),temps,'nearest');
else
    % utilsiation de tequila
    [jpiq,tjpiq]         = tsbase(shot,'SJPIQ');
    if isempty(jpiq)
	z0dinput.option.li = 1;
	z0dinput.option.limode = 1;
	librute               = ones(size(temps));

    else
	lipiq                 = log(1.65 + 0.89 .* jpiq);
	li                    = interp10d(tjpiq,lipiq,temps,'nearest');
	librute               = li;
	li(~isfinite(li))     = 1;
	li                    = min(10,max(0.1,li));
	li                    = medfilt1(li,7);
	z0dinput.option.li      = li(7);
    end

end
% donnee calculee dans le zerod
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];
	
	
% les parametres
[pant,tpant]          = tsbase(shot,'gpuifci');
if ~isempty(pant)
  indant              = find(max(pant(:,1:3),[],1) > 0.3);
  if isempty(indant)
    indant            = 1;
  end
  frequence           = tsbase(shot,'sfreqfci');
  if ~isempty(frequence)
    z0dinput.option.freq = max(1,mean(frequence(indant)));
  end
end
z0dinput.option.mino     = 'H';
geom.a                   = z0dinput.geo.a;
geom.r0                  = z0dinput.geo.R;
geom.b0                  = z0dinput.geo.b0;
geom.d0                  = d0;
pos                      = geom.r0 + geom.a + 0.02;
[scenar,scenstr,Rres]    = scenarzerodfci(geom,z0dinput.option.freq,pos,z0dinput.option.mino);
if isfinite(Rres)
  z0dinput.option.fwcd = 0;
else
  z0dinput.option.fwcd = 2;
end
z0dinput.option.angle_nbi = 90;
z0dinput.option.modeh     = 0;
z0dinput.option.rip        = 1;
z0dinput.option.nphi = 30;
z0dinput.option.etalh     = 0.8;
z0dinput.option.lhmode    = 3;
z0dinput.option.xlh        = 0.2;
z0dinput.option.dlh        = 0.3;
z0dinput.option.zmax       = 8;
z0dinput.option.zimp       = 6;
z0dinput.option.rimp       = 0.3;
z0dinput.option.zeff       = 0;
z0dinput.option.ane        = 0;
z0dinput.option.qdds      = 0.95;
z0dinput.option.kidds     = 3;
z0dinput.option.xiioxie   = 0.5;


if shot <= 28353
	z0dinput.option.signe = 1;
else
	z0dinput.option.signe = 1;	
end

% rapport iso
[rgaz,rdcx] = lit_nhnd_ts(shot);
if ~isempty(rdcx)
	z0dinput.option.cmin = min(0.3,max(0.01,rdcx ./ (1 - rdcx)));
	z0dinput.option.mino = 'H';
elseif ~isempty(rgaz)
	z0dinput.option.cmin = min(0.3,max(0.01,rgaz ./ (1 - rgaz)));	
	z0dinput.option.mino = 'H';
end
% on tente le temps reel 'GPROFILEVSPX'
[g3,t3]=tsbase(shot,'GPROFILEVSPX');
if ~isempty(g3)
	[valid,tvalid]=tsbase(shot,'GPROFCHKVSPX%1'); 
	if ~isempty(valid)
		ind = find(valid >= 0.999);
		fprintf('SPX : %d invalid measurement on %d\n',length(ind),size(g3,1));
		if size(valid,1) == size(g3,1)
			g3(ind,:) = [];
			t3(ind,:) = [];	
		else
			g3 = [];
			t3 = [];
		end
	end
end
if isempty(g3) || (size(g3,1) < 7) || all(g3(:) == 0)
	  % recherche des xdur
	  for k = 0:9
		[g3,t3,d3,c3]=tsbase(shot+k/10,'GHXR3');
		if ~isempty(g3)
			break
		end
	  end
else
		t3 = t3(:,1);
		d3 = linspace(0,1,20);
		g3 = g3(:,1:20);	
end

% si presence hybride
phyb  = [];
tphyb =[];
if fix(shot)<20000   
	[phyb,tphyb,void,cert] = tsbase(fix(shot),'ghyb');
else
	[phyb,tphyb,void,cert] = tsbase(fix(shot),'gphyb');
	if size(tphyb,2) == 1
		tphyb = tphyb * ones(1,size(phyb,2));
	end
end
if (fix(shot)<20000) & isempty(phyb) 
	[phyb,tphyb,void,cert] = tsbase(fix(shot),'gphyb');
	if size(tphyb,2) == 1
		tphyb = tphyb * ones(1,size(phyb,2));
	end
end 
if ~isempty(phyb)
[npar01,npar02]= litphashyb(shot,tphyb,phyb);
npar0          =  npar01 .* phyb(:,1) +  npar02 .* phyb(:,2);
phybt          =  phyb(:,1) +  phyb(:,2);
indok          =  find(isfinite(npar0));
if isempty(indok)
    npar0 = 1.8;
else
    npar0 = sum(npar0(indok)) ./ sum(phybt(indok)); 
end
else
npar0 = 1.8;
end
if ~isempty(g3) & ~all(g3(:) == 0)
    
	indbad = find(diff(t3) <=0);
	while (~isempty(indbad) & ~isempty(t3))
		g3(indbad,:) = [];
		t3(indbad)   = [];
		indbad = find(diff(t3) <=0);			
	end


	XDUR.v = g3;
	XDUR.t = t3;
	XDUR.x = d3;

	g3s  = sum(g3,1);
	%xlh  = trapz(d3,(1-d3) .* d3.* g3s,2) ./trapz(d3,(1-d3) .* g3s,2);
	xlh  = d3(max(find(g3s ==max(g3s))));
	dlh  = sqrt(trapz(d3,d3 .^ 2 .* g3s,2) ./max(eps,trapz(d3,g3s,2)) - xlh.^ 2);
	z0dinput.option.xlh = xlh;
	z0dinput.option.dlh = dlh;
	z0dinput.option.lhmode = 3;
	z0dinput.option.npar0 = npar0;
	z0dinput.option.etalh  = min(1,max(0.1,2.01 - 0.63 .* npar0));
	z0dinput.option.wlh = 0;
	z0dinput.exp0d.XDURt = t3;
	z0dinput.exp0d.XDURx = d3;
	z0dinput.exp0d.XDURv = g3;
	
    

else
	z0dinput.option.npar0 =npar0;
	z0dinput.option.wlh = 11e-3 .* 32;
	z0dinput.option.lhmode = 3;
end

% lecture des donnees ecrh directement dans le diagnostic
[xika1,tbad]=tsbase(shot,'sika1');  		% gyrotron A1 cathode current
[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
if ~isempty(xika1)  & ~isempty(xHAUTTOR)
	[prxA1,tA1]=tsbase(shot,'spia1');   		% gyrotron A1 power
	[prxA2,tA2]=tsbase(shot,'spia2');   		% gyrotron A2 power
	[sonde,tsonde]=tsbase(shot,'sonderf');   	% RF probe on reflectometer
	%[xHAUTTOR,tHAUTTOR]=tsbase(shot,'SHAUTTOR');		% Toroidal injection anlge top mirror (A1) hysteresis corrected
	[xHAUTPOL,tHAUTPOL]=tsbase(shot,'SHAUTPOL');		% Poloidal injection anlge top mirror (A1) hysteresis corrected
	[xMILTOR,tMILTOR]=tsbase(shot,'SMILTOR');		% Toroidal injection anlge central mirror (A2) hysteresis corrected
	[xMILPOL,tMILPOL]=tsbase(shot,'SMILPOL');		% Poloidal injection anlge central mirror (A2) hysteresis corrected
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% powers in kW
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	factA1=1; 		% Calibration factor to be applied to  A1
	factA2=1; 		% Calibration factor to be applied to  A2
	PA1=max(0,1e6*prxA1/factA1);  % A1 power in W
	PA2=max(0,1e6*prxA2/factA2);	% A2 power in W
	Ptot=PA1+PA2;		% total power
	
	% position
	rA2  = 3.5300;
	rA1  = 3.5300;
	zA2  = 0;
	zA1  = 0.2000;

	% puissance
	z0dinput.cons.pecrh      = max(0,interp10d(tbad,Ptot,temps,'nearest'));
	z0dinput.cons.pecrh      = z0dinput.cons.pecrh .* (z0dinput.cons.pecrh >=1e4);
	% position 
	phi_tor_a1 = interp10d(tbad,xHAUTTOR,temps,'nearest')./180*pi;	
	phi_tor_a2 = interp10d(tbad,xMILTOR,temps,'nearest')./180*pi;	
	phi_pol_a1 = interp10d(tbad,xHAUTPOL,temps,'nearest')./180*pi;	
	phi_pol_a2 = interp10d(tbad,xMILPOL,temps,'nearest')./180*pi;	
	
	% calcul de la position
	[xece1,nh1,apol1] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,d0,z0dinput.geo.a,rA1,zA1,phi_pol_a1,phi_tor_a1,z0dinput.geo.b0,118);
	[xece2,nh2,apol2] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,d0,z0dinput.geo.a,rA2,zA2,phi_pol_a2,phi_tor_a2,z0dinput.geo.b0,118);
	%[xece,nh,apol] = r2xece(r0,z0,d0,a,rant,zant,phi_pol,phi_tor,b0,freq_ghz)
	
	% moyenne
	pa1 = interp10d(tbad,PA1,temps,'nearest');	
	pa2 = interp10d(tbad,PA2,temps,'nearest');		
	z0dinput.cons.xece   = max(0,min(1,(xece1 .* pa1 + xece2 .* pa2) ./ max(1,pa1 + pa2)));
	z0dinput.option.sens = sign(trapz(temps,phi_tor_a1 .* pa1 + phi_tor_a2 .* pa2)) ;
	z0dinput.option.angle_ece = trapz(temps,abs(apol1) .* pa1 + abs(apol2) .* pa2) ./ max(1,trapz(temps,pa1 +pa2));
end


if gaz == 1
	z0dinput.option.gaz = 2;
else
	z0dinput.option.gaz = 4;
end
z0dinput.option.zmax  = 8;

% read gaz puff reference
[sdeb,tdeb]           = tsbase(shot,'sdeb');
if ~isempty(tdeb)
  if z0dinput.option.gaz  == 4
	sdeb = (4.41e-4 .* 6.02214199e23) .* sdeb;   % conversion to electron/s from  Pa.m^3/s
  else 
	sdeb = 2 .* (4.41e-4 .* 6.02214199e23) .* sdeb;   % conversion to electron/s from  Pa.m^3/s  
  end
  sdeb = sgolayfilt(sdeb,1,51);
  sdeb(sdeb < 0) = 0;
  z0dinput.cons.nbar  = z0dinput.cons.nbar + sqrt(-1) .* interp10d(tdeb,sdeb,temps,'nearest');  
end
% prefill pressure 
[gptore,tptore]       = tsbase(shot,'gptore');
if ~isempty(gptore)
  pressure   = gptore(:,1) .* 10 .^ gptore(:,2);
  indpress = find((tptore(:,1) < 0) & (tptore(:,1) >= -0.2));
  z0dinput.option.p_prefill = mean(pressure(indpress));
end
% breakdown data
if shot > 36599
   % ciel
   z0dinput.option.L_eddy = 1.1327e-05;
   z0dinput.option.R_eddy = 0.0095;
   z0dinput.option.berror = 0; % not a measure !
end

% donnees idn 
[tcronos,pnbi_cronos,beam_energy,frac,amass] = read_nbi_power(shot,temps);
z0dinput.cons.pnbi = pnbi_cronos;
z0dinput.option.rtang = 1.47;
z0dinput.option.angle_nbi = 90;
z0dinput.option.einj = 1e3 .* beam_energy .* (frac(1) + frac(2) / 2 + frac(3) / 3);
if amass == 2
    z0dinput.cons.ftnbi = zeros(size(temps));
else
    z0dinput.cons.ftnbi = ones(size(temps));
end 

% donnees experimentales
z0dinput.exp0d.temps = temps;
z0dinput.exp0d.pin   = z0dinput.cons.picrh + z0dinput.cons.plh + z0dinput.cons.pecrh + pohm;
		      
z0dinput.exp0d.ploss = ptot;
if ~isempty(zeffm)
	z0dinput.exp0d.zeff  = z0dinput.cons.zeff;
end
z0dinput.exp0d.ane  = ane;
z0dinput.exp0d.nem  = nemoy;
z0dinput.exp0d.ne0  = ne0;
if~isempty(gbilan)
	  z0dinput.exp0d.w     = interp10d(tbilan,gbilan(:,6),temps,'nearest') .* 1e6;
else
	  [wdia,tmag] =   tsbase(shot,'swdia'); 
	  z0dinput.exp0d.w     = interp10d(tmag,wdia,temps,'nearest').* 1e6;
end
z0dinput.exp0d.dwdt  = pdederive(temps,z0dinput.exp0d.w,2,2,1,1); 
if ~isempty(gbilan)
	  z0dinput.exp0d.taue  = interp10d(tbilan,gbilan(:,7),temps,'nearest');
else
	  z0dinput.exp0d.taue  = NaN .* ones(size(temps));
end
z0dinput.exp0d.pw    = ptot;

[gshr,tsh]         = tsbase(shot,'gshr');
[gshte,tsh]      = tsbase(shot,'gshtenv');
if isempty(gshte)
	[gshte,tsh]      = tsbase(shot,'gshte');
end
if isempty(gshr)
	[gshte,tsh]      = tsbase(shot,'gtrshte');
	[gshr,tsh]      = tsbase(shot,'gtrshr');
	if ~isempty(tsh)
		tsh             = tsh(:,1);
		gshte           = gshte ./ 1000;
		gshr           = gshr ./ 1000;
	end
end
if ~isempty(gshte)
      indtok             = find(tsh >2 & tsh < (max(tsh) -1));
      indok              = find(all(gshte(indtok,:)>=0,1));
      gshr               = mean(gshr(:,indok),1);
      gshte              = gshte(:,indok);
      indok              = find(gshr > 2.35 & gshr <= 2.45);
      if isempty(indok)
	      d     = abs(gshr - 2.4);
	      indok = max(find(d == min(d)));
      end
      te0                = mean(gshte(:,indok),2) .* 1e3;
      tshc = tsh;
      te0c = te0;
      indr = find(diff(tsh)<=0);
      if ~isempty(indr)
	      tshc(indr) = [];
	      te0c(indr) = [];
      end
      z0dinput.exp0d.te0   = interp10d(tshc,te0c,temps,'nearest');
      [ti,tti]           = tsbase(shot,'stibrag');
      if ~isempty(ti)
	      [teb,tti]           = tsbase(shot,'stebrag');
	      z0dinput.exp0d.tite   = interp10d(tti,ti./teb,temps,'nearest');
      end
      nom=tsbase(shot,'SCIONBRAG');
      [deplw,tw]=tsbase(shot,'SDEPLW');
      if ~isempty(deplw)
	    if strcmp(lower(deblank(nom)),'fer')
		  fact =3;
	    elseif  strcmp(lower(deblank(nom)),'chrome')
		  fact = 2.5;
	    else
		  fact =NaN;
	    end
	    % Ip sens inverse trigo & deplw >0 sens trigo , choc >28424
	    if shot >28424
	      z0dinput.exp0d.wrad =  - interp10d(tw,deplw .* 1e5,temps,'nearest') .* fact .* 1e3 ./ z0dinput.geo.R; 
	    else
	      z0dinput.exp0d.wrad =  - interp10d(tw,deplw .* 1e5,temps,'nearest') .* fact .* 1e3 ./ z0dinput.geo.R; 
	    end
	    z0dinput.exp0d.wrad =  z0dinput.exp0d.wrad  ./ (1 +  z0dinput.exp0d.ane);  % valeur moyenne approximative
	end
else
      [tethom,tthom,zthom]      = tsbase(shot,'GTETHOM');
      if isempty(tethom)
	      [temic,tmic] = tsbase(shot,'gmictenv');
	      if isempty(temic)
		      [temic,tmic] = tsbase(shot,'gmicte');
		      if isempty(temic)
			      disp('No temperature measurement available')
		      end
	      end
	      if ~isempty(temic)
		      [rmic,tmic,voies,cert] = tsbase(shot,'gmicrnv');
		      indtok             = find(tmic >2 & tmic < (max(tsh) -1));
		      indok              = find(all(temic(indtok,:)>=0,1));
		      rmic               = mean(rmic(:,indok),1);
		      temic              = gshte(:,indok);
		      indok              = find(rmic > 2.35 & rmic <= 2.45);
		      if isempty(indok)
			      d     = abs(rmic - 2.4);
			      indok = max(find(d == min(d)));
		      end
		      te0                = mean(temic(:,indok),2) .* 1e3;
		      tshc = tmic;
		      te0c = te0;
		      indr = find(diff(tsh)<=0);
		      if ~isempty(indr)
			      tshc(indr) = [];
			      te0c(indr) = [];
		      end
		      z0dinput.exp0d.te0   = interp10d(tshc,te0c,temps,'nearest');
	      else
		      [tefab,tfab] = tsbase(shot,'gtefab');
		      if isempty(tefab)
			      disp('No temperature measurement available')
		      else
			      z0dinput.exp0d.te0   = interp10d(tfab,max(tefab .* (tefab >0),[],2),temps,'nearest').* 1e3;	
		      end

		      
	      end
      else

	      indok              = find(abs(zthom) == min(abs(zthom)));
	      te0                = mean(tethom(:,indok),2) .* 1e3;
	      z0dinput.exp0d.te0   = interp10d(tthom,te0,temps,'nearest');
      end
end
	
z0dinput.exp0d.pohm  = pohm;
%        z0dinput.exp0d.vloop  = interp10d(tbilan,gbilan(:,9),temps,'nearest') .* 1e6;
if size(gbilan,2) < 9
	[vl,tvl]           = tsbase(shot,'svsur');
	if isempty(vl)
		[vl,tvl]           = tsbase(shot,'sv217');
	end
	if isempty(vl)
		[vl,tvl]           = tsbase(shot,'sv235');
	end
%
% probleme tvl
%	
	indtvl = find(diff(tvl) <= 0);
	if ~isempty(indtvl)
		tvl(indtvl) = [];
		vl(indtvl) = [];
	end

	z0dinput.exp0d.vloop = interp10d(tvl,medfilt1(vl,11),temps,'nearest');
else
	z0dinput.exp0d.vloop = interp10d(tbilan,gbilan(:,9),temps,'nearest');
end
if ~isempty(gqbeli)
    z0dinput.exp0d.qa    = interp10d(tgqbeli,gqbeli(:,1),temps,'nearest'); 
    z0dinput.exp0d.betap = interp10d(tgqbeli,gqbeli(:,2),temps,'nearest'); 
else
    [qa,tqa]=tsbase(shot,'sqpsi');
    if ~isempty(qa)
	  z0dinput.exp0d.qa    = interp10d(tqa,qa,temps,'nearest'); 
    end
    [beta,tb] =tsbase(shot,'SDIAM'); 
    if ~isempty(beta)
	  z0dinput.exp0d.betap    = interp10d(tb,beta,temps,'nearest'); 
    end

end
z0dinput.exp0d.ip    = z0dinput.cons.ip;

[prad,tprad]       = tsbase(shot,'sprad');
if isempty(prad)
	[prad,tprad]       = tsbase(shot,'gbilan%11');        
end
if ~isempty(prad)
	z0dinput.exp0d.prad  = interp10d(tprad,prad,temps,'nearest') .* 1e6;
end
z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
z0dinput.exp0d.li    = librute;
z0dinput.exp0d.picrh = z0dinput.cons.picrh;
z0dinput.exp0d.plh   = z0dinput.cons.plh;
z0dinput.exp0d.pecrh = z0dinput.cons.pecrh;

z0dinput.exp0d.edgeflux    = z0dinput.cons.flux;
z0dinput.machine     = 'TS';
z0dinput.shot        = shot;
      
z0dinput.option.available_flux = 10.3 - (-10.1);  % usual 9.8  - (-7.82)

%============================================================
% acces au donnees de Tore Supra avant choc (preparation)
%============================================================
function z0dinput = zerod_init_ts_preshot(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
% cas donnees TS preparation
z0dinput.exp0d         = [];

% 1 - numero du choc
if isempty(gaz) & isempty(shot)

    numshot = 0;
    refshot = 0;
    try 
	refshot = evalin('base','param.from.shot.num');
    catch
	try
		refshot = evalin('base','post.z0dinput.shot');
	catch
		refshot = 0;		
	end
    end 

    switch langue 
    case 'francais'
	    prompt={'Numero du choc :','Choc de reference :'};
	    def={sprintf('%d',numshot),sprintf('%d',refshot)};
	    dlgTitle='Lecture des donnees de Tore Supra';
    otherwise
	    prompt={'shot number :','reference shot number :'};
	    def={sprintf('%d',numshot),sprintf('%d',refshot)};
	    dlgTitle='Access to Tore Supra data';
    end
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
            z0dinput = [];
	    return
    end
    shot     = str2num(answer{1});
    refshot  = str2num(answer{2});
else
      refshot = gaz;	
end
%function checkup
% Check-up du pilote avant choc

% *** lecture des consignes principales
% APOLO
Ip   = tsmat(shot,'APOLO;Ip_Ref;VALEURS');
if isempty(Ip)
	disp('No plasma') 
        z0dinput = [];
	return
end
R    = tsmat(shot,'APOLO;Pos_R;VALEURS');
Z    = tsmat(shot,'APOLO;Pos_Z;VALEURS');
El    = tsmat(shot,'APOLO;Pos_El;VALEURS');
Tr    = tsmat(shot,'APOLO;Pos_T;VALEURS');
Idiv = tsmat(shot,'APOLO;I_Div;VALEURS');
Flag_I0 = tsmat(shot,'APOLO;+VME_POL;Flag_I0');
Flag_Lang = tsmat(shot,'APOLO;+VME_POL;Flag_Lang');

% AGAZ
Nl   = tsmat(shot,'AGAZ;Nl;VALEURS');
D2   = tsmat(shot,'AGAZ;Prop_D2;VALEURS');
if ~isempty(D2)     
	D2 = max(D2(:,2));
end
He   = tsmat(shot,'AGAZ;Prop_He;VALEURS');     
if ~isempty(He)     
	He = max(He(:,2));
end
if isempty(D2) & isempty(He)
	D2 = 1;
	He = 0;
end	
if isempty(He)
	He = 0;
end	
if isempty(D2)
	D2 = 0;
end

% AHYB
phybtot  = tsmat(shot,'AHYB;P_Hyb_Tot;VALEURS');
phybc1   = tsmat(shot,'AHYB;P_Cp1;VALEURS');
phybc2   = tsmat(shot,'AHYB;P_Cp2;VALEURS');
phyb_dis = tsmat(shot,'AHYB;+VME_HYB;Dispatch_P');
ahybtot  = tsmat(shot,'AHYB;Ph_Globale;VALEURS');
ahybc1   = tsmat(shot,'AHYB;Ph_Cp1;VALEURS');
ahybc2   = tsmat(shot,'AHYB;Ph_Cp2;VALEURS');

% AFCI
pfci1 = tsmat(shot,'AFCI;Consigne_Q1;VALEURS');
pfci4 = tsmat(shot,'AFCI;Consigne_Q2;VALEURS');
pfci5 = tsmat(shot,'AFCI;Consigne_Q5;VALEURS');
pfcip  = tsmat(shot,'AFCI;+VME_FCI;Puissance');

% AFCE
pfce1 = tsmat(shot,'AFCE;Power_A1;VALEURS');
pfce2 = tsmat(shot,'AFCE;Power_A2;VALEURS');
%
onoffa1 = tsmat(shot,'DFCE;PlasmaParamA;EnableA1');
onoffa2 = tsmat(shot,'DFCE;PlasmaParamA;EnableA2');
%
angvara1 = tsmat(shot,'DFCE;PlasmaParamA;EnableRTTop');
angvara2 = tsmat(shot,'DFCE;PlasmaParamA;EnableRTMid');
%
duty1 = tsmat(shot,'AFCE;Mod_Duty_A1;VALEURS');
duty2 = tsmat(shot,'AFCE;Mod_Duty_A2;VALEURS');
freq1 = tsmat(shot,'AFCE;Mod_Freq_A1;VALEURS');
freq2 = tsmat(shot,'AFCE;Mod_Freq_A2;VALEURS');
pol1 = tsmat(shot,'AFCE;PolAng_Haut;VALEURS');
pol2 = tsmat(shot,'AFCE;PolAng_Mil;VALEURS');
tor1 = tsmat(shot,'AFCE;TorAng_Haut;VALEURS');
tor2 = tsmat(shot,'AFCE;TorAng_Mil;VALEURS');	
%

% *** lecture donnees controle commande
Itor = tsmat(shot,'EXP=T=S;GENERAL;ITOR');
if Itor == 0
	disp('No plasma')
        z0dinput = [];
	return
end
Lim  = tsmat(shot,'EXP=T=S;LIMITEURS;LIMITEUR');
hyb  = tsmat(shot,'EXP=T=S;ANTENNES;POSHYB');
fci  = tsmat(shot,'EXP=T=S;ANTENNES;POSCI');
chochyb  = tsmat(shot,'EXP=T=S;CHAUFFAGES;CHOCHYB');
chocfci  = tsmat(shot,'EXP=T=S;CHAUFFAGES;CHOCFCI');
chocfce  = tsmat(shot,'EXP=T=S;CHAUFFAGES;CHOCFCE');
chocidn  = tsmat(shot,'EXP=T=S;CHAUFFAGES;CHOCIDN');

% creation de la basetemps commune
dtu = Inf;
t = cat(1,R(:,1),Z(:,1),Ip(:,1),Nl(:,1),El(:,1),Tr(:,1));
if ~isempty(pfci1)
	t = cat(1,t,pfci1(:,1));
end
if ~isempty(pfci4)
	t = cat(1,t,pfci4(:,1));
end
if ~isempty(pfci5)
	t = cat(1,t,pfci5(:,1));
end
if ~isempty(pfce1)
	if onoffa1 == 1
		t = cat(1,t,pfce1(:,1));
	end
end
if ~isempty(duty1)
	if onoffa1 == 1
		t = cat(1,t,duty1(:,1));
	end
end
if ~isempty(freq1)
	if onoffa1 == 1
		t = cat(1,t,freq1(:,1));
		dtu = min(dtu,min(duty1(:,2))./max(freq1(:,2)) ./ 5);
	end
end
if ~isempty(pol1)
	if onoffa1 == 1
		t = cat(1,t,pol1(:,1));
	end
end
if ~isempty(tor1)
	if onoffa1 == 1
		t = cat(1,t,tor1(:,1));
	end
end
if ~isempty(pfce2)
	if onoffa2 == 1
		t = cat(1,t,pfce2(:,1));
	end
end
if ~isempty(duty2)
	if onoffa2 == 1
		t = cat(1,t,duty2(:,1));
	end
end
if ~isempty(freq2)
	if onoffa2 == 1
		t = cat(1,t,freq2(:,1));
		dtu = min(dtu,min(duty2(:,2))./max(freq2(:,2)) ./ 5);
	end
end
if ~isempty(pol2)
	if onoffa2 == 1
		t = cat(1,t,pol2(:,1));
	end
end
if ~isempty(tor2)
	if onoffa2 == 1
		t = cat(1,t,tor2(:,1));
	end
end
if ~isempty(phybtot)
	t = cat(1,t,phybtot(:,1));
end
if ~isempty(phybc1)
	t = cat(1,t,phybc1(:,1));
end
if ~isempty(phybc2)
	t = cat(1,t,phybc2(:,1));
end
if ~isempty(ahybtot)
	t = cat(1,t,phybtot(:,1));
end
if ~isempty(ahybc1)
	t = cat(1,t,phybc1(:,1));
end
if ~isempty(ahybc2)
	t = cat(1,t,phybc2(:,1));
end

t = cat(1,t,1000);
t = sort(t);
while (any(diff(t)<=0))
	t(find(diff(t)==0)+1)=t(find(diff(t)==0)+1) + 1e-6;
end



% calcul de la puissance LH
if ~isempty(phybtot)
	phybtoti=interp1c([phybtot(:,1);1000],[phybtot(:,2);phybtot(size(phybtot,1),2)],t);
	phybc1i=interp1c([phybc1(:,1);1000],[phybc1(:,2);phybc1(size(phybc1,1),2)],t);
	phybc2i=interp1c([phybc2(:,1);1000],[phybc2(:,2);phybc2(size(phybc2,1),2)],t);
	phyb=(phyb_dis(1)+phyb_dis(2))*phybtoti+phyb_dis(3)*phybc1i+phyb_dis(4)*phybc2i;
	ahybtoti=interp1c([ahybtot(:,1);1000],[ahybtot(:,2);ahybtot(size(ahybtot,1),2)],t);
	ahybc1i=interp1c([ahybc1(:,1);1000],[ahybc1(:,2);ahybc1(size(ahybc1,1),2)],t);
	ahybc2i=interp1c([ahybc2(:,1);1000],[ahybc2(:,2);ahybc2(size(ahybc2,1),2)],t);
	nparc1i = 1.83 + ahybc1i ./ 200;
	nparc2i = 2.03 + (ahybc2i + 90) ./ 310;
	npari   = (phyb_dis(1)*phybtoti	+ phyb_dis(3)*phybc1i) .* nparc1i + ...
		  (phyb_dis(2)*phybtoti	+ phyb_dis(4)*phybc2i) .* nparc2i;
	npari   = trapz(t,npari) ./ max(1e-6,trapz(t,phyb));
else
	phyb = zeros(size(t));
	npari =[];
end

% calcul de la puissance FCI
pfci = zeros(size(t));
if ~isempty(pfci1)
	pfci = pfci + interp1c([pfci1(:,1);1000],[pfci1(:,2);pfci1(size(pfci1,1),2)],t);
end
if ~isempty(pfci4)
	pfci = pfci + interp1c([pfci4(:,1);1000],[pfci4(:,2);pfci4(size(pfci4,1),2)],t);
end
if ~isempty(pfci5)
	pfci = pfci + interp1c([pfci5(:,1);1000],[pfci5(:,2);pfci5(size(pfci5,1),2)],t);
end


% calcul de la puissance FCE
pfce = zeros(size(t));
if ~isempty(pfce1)
	if onoffa1 == 1
		pfce1 = 0.3 .* interp1c([pfce1(:,1);1000],[pfce1(:,2);pfce1(size(pfce1,1),2)],t);
		tor1  = interp1c([tor1(:,1);1000],[tor1(:,2);tor1(size(tor1,1),2)],t);
		pol1  = interp1c([pol1(:,1);1000],[pol1(:,2);pol1(size(pol1,1),2)],t);
		duty1  = interp1c([duty1(:,1);1000],[duty1(:,2);duty1(size(duty1,1),2)],t);
		freq1  = interp1c([freq1(:,1);1000],[freq1(:,2);freq1(size(freq1,1),2)],t);
		pfce = pfce + pfce1;
	end
end
if ~isempty(pfce2)
	if onoffa2 == 1
		pfce2 = 0.3 .* interp1c([pfce2(:,1);1000],[pfce2(:,2);pfce2(size(pfce2,1),2)],t);
		tor2  = interp1c([tor2(:,1);1000],[tor2(:,2);tor2(size(tor2,1),2)],t);
		pol2  = interp1c([pol2(:,1);1000],[pol2(:,2);pol2(size(pol2,1),2)],t);
		duty2  = interp1c([duty2(:,1);1000],[duty2(:,2);duty2(size(duty2,1),2)],t);
		freq2  = interp1c([freq2(:,1);1000],[freq2(:,2);freq2(size(freq2,1),2)],t);
		pfce = pfce + pfce2;
	end
end

Rlim=Lim(1); 
Zlim=max(Lim(find(Lim <0))); 
C1=hyb(1); 
C2=hyb(2); 
Q1=fci(1); 
Q4=fci(2); 
Q5=fci(3);

Ipi=interp1c([Ip(:,1);1000],[Ip(:,2);Ip(size(Ip,1),2)],t);
Ri=interp1c([R(:,1);1000],[R(:,2);R(size(R,1),2)],t);
Zi=interp1c([Z(:,1);1000],[Z(:,2);Z(size(Z,1),2)],t);
Eli=interp1c([El(:,1);1000],[El(:,2);El(size(El,1),2)],t);
Tri=interp1c([Tr(:,1);1000],[Tr(:,2);Tr(size(Tr,1),2)],t);
Nli=interp1c([Nl(:,1);1000],[Nl(:,2);Nl(size(Nl,1),2)],t);

dIp = diff(Ipi)./diff(t)*1000; dIp=[0;dIp];

Appui = ['   PPI';'   LPV';'   LPM';'   C1 ';'   C2 ';'   Q1 ';'   Q4 ';'   Q5 '];
Rppi=1.556;
dobjet = [Ri-Rppi Zlim-Zi Rlim-Ri C1-Ri C2-Ri Q1-Ri Q4-Ri Q5-Ri];
a=abs(min(dobjet'))'; 
for i=1:max(size(t))
    appui(i)=find(dobjet(i,:)==min(dobjet(i,:))); 
end

% calcul de K et d
K = 1 + Eli ./ a;
d = Tri ./ a;

% gaz
if D2 
	gaz = 2;
elseif He 
	gaz = 4;
else
	gaz = 1;
end


% creation de la base temps
tip = cat(2,0:0.01:1,1.1:0.1:1000)';
ip  = interp1(cat(1,-1,t),cat(1,0,Ipi),tip);

% lecture des donnees
tt = tip(tip >=0 & ip > 3e-2);
if isempty(tt)
	disp('No plasma current');
        z0dinput = [];
	return
end


if isempty(temps)
	dt    = max(mean(diff(tt)),0.1);	
	temps = (min(tt):dt:max(tt))';
	if length(temps) > 1001
		temps = linspace(min(tt),max(tt),1001)';
	end	
	
	if isfinite(dtu)
		% modulation fce detectee
		tp = (min(tt(pfce>0)):dtu:max(tt(pfce>0)))';	
		temps = cat(1,temps,tp);
		temps = sort(temps);
		ind   = find(diff(temps) < 1e-3) + 1;
		while(~isempty(ind))
			temps(ind) = [];
			ind   = find(diff(temps) < 1e-3) + 1;
		end
		dtu	
	end
end

z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = max(1,interp10d(tip,ip,temps,'nearest') .* 1e6);
z0dinput.cons.flux     = zeros(size(temps));
z0dinput.geo.a         = interp10d(t,a,temps,'nearest');
z0dinput.geo.R         = interp10d(t,Ri,temps,'nearest');
z0dinput.geo.z0        = interp10d(t,Zi,temps,'nearest');
z0dinput.geo.K         = interp10d(t,K,temps,'nearest');
z0dinput.geo.d         = interp10d(t,d,temps,'nearest');

rb0                    = (4*pi*1e-7) .* 18 .* 2028 .* Itor ./ 2 ./ pi
z0dinput.geo.b0        = rb0 ./ z0dinput.geo.R;
z0dinput.cons.nbar     = pchip(t,max(0,Nli.*1e19),temps)./ 2 ./ z0dinput.geo.a;

if chocfci > 0
	z0dinput.cons.picrh    =  max(0,interp10d(t,pfci,temps,'nearest')) .* 1e6;
else
	z0dinput.cons.picrh    =  zeros(size(temps)); 
end
if chochyb > 0
	z0dinput.cons.plh      = max(0,interp10d(t,phyb,temps,'nearest')) .* 1e6;
else
	z0dinput.cons.plh      = zeros(size(temps)); 
end 
if chocidn > 0
	% letcture impossible ....
	% debut sequence IDN dans TOPIDN de APILOTE;Choc_Plasma;ACHRONO
	% chaque tir avec duree dans IDNON (plusieurs possibles) de de APILOTE;__IDN;ACHRONO
	z0dinput.cons.pnbi     = zeros(size(temps)); 
else
	z0dinput.cons.pnbi     = zeros(size(temps)); 
end
if chocfce > 0
	z0dinput.cons.pecrh    = max(0,interp10d(t,pfce,temps,'nearest')) .* 1e6;
else
	z0dinput.cons.pecrh    = zeros(size(temps)); 
end
z0dinput.cons.hmore    = ones(size(temps));
z0dinput.cons.zeff    = (2.3 + (gaz == 4)) .* ones(size(temps));
z0dinput.option.li     = 0.5;
z0dinput.option.mino   = 'H';

% securite densite
nsat = 0.06e20 .* (z0dinput.cons.ip ./ 1e6) .* z0dinput.geo.R .* sqrt(gaz) ./ z0dinput.geo.a .^ (5/2);
indbad = find(z0dinput.cons.nbar <= 3e17); 
indgood = find(cumsum(Nli) == max(cumsum(Nli)),1);
if ~isempty(indgood)
    indlast  = find(temps >= t(indgood),1);
    if indlast > 1
	z0dinput.cons.nbar(indlast:end) = nsat(indlast:end) ./nsat(indlast-1) .* z0dinput.cons.nbar(indlast-1);
	%z0dinput.cons.nbar =  z0dinput.cons.nbar ./ (0.3  + 0.7 .*  max(1,z0dinput.cons.nbar ./ nsat));
    end
end
	
% donnee calculee dans le zerod
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];
	
	
% les parametres
ffci = tsmat(shot,'DFCI;PILOTAGE;FREQUENC');
if ~isempty(ffci)
	z0dinput.option.freq = max(1,mean(ffci));
else
	z0dinput.option.freq = 57;
end
z0dinput.option.mino     = 'H';
geom.a                   = z0dinput.geo.a;
geom.r0                  = z0dinput.geo.R;
geom.b0                  = z0dinput.geo.b0;
geom.d0                  = 0.1 .* ones(size(z0dinput.geo.a));
pos                      = geom.r0 + geom.a + 0.02;
[scenar,scenstr,Rres]    = scenarzerodfci(geom,z0dinput.option.freq,pos,z0dinput.option.mino);
if isfinite(Rres)
  z0dinput.option.fwcd = 0;
else
  z0dinput.option.fwcd = 2;
end
z0dinput.option.fwcd = 0;
z0dinput.option.angle_nbi = 0;
z0dinput.option.modeh     = 0;
z0dinput.option.rip       = 1;
z0dinput.option.nphi      = 30;
z0dinput.option.etalh     = 0.7;
z0dinput.option.lhmode    = 3;
z0dinput.option.xlh       = 0.2;
z0dinput.option.dlh       = 0.3;
z0dinput.option.zmax      = 8;
z0dinput.option.zimp      = 6;
z0dinput.option.rimp      = 0.3;
z0dinput.option.zeff      = 7;
z0dinput.option.ane       = 0;
z0dinput.option.neasser   = 1;
z0dinput.option.Recycling = 0.7;
z0dinput.option.natural   = 1;

z0dinput.option.cmin      = 0.05;

z0dinput.option.gaz       = gaz;
z0dinput.option.matthews  = 1;
z0dinput.option.qdds      = 0.95;
z0dinput.option.kidds     = 3;
z0dinput.option.xiioxie   = 0.5;
z0dinput.option.scaling   = 0;


if (shot <= 28353) & (shot > 0)
	z0dinput.option.signe = 1;
else
	z0dinput.option.signe = 1;	
end

z0dinput.option.npar0 = 1.8;
if ~isempty(npari)
	z0dinput.option.npar0 = npari;	
end

z0dinput.machine     = 'TS';
z0dinput.shot        = shot;

% reglage de la postion de fce		
xece = zeros(size(temps));
apol = zeros(size(temps));
ator = zeros(size(temps));
if ~isempty(pfce1)
	if onoffa1 == 1
		rA1  = 3.5300;
		zA1  = 0.2000;
		pfce1 = interp10d(t,pfce1,temps,'nearest');
		tor1 = interp10d(t,tor1,temps,'nearest')./180*pi;
		pol1 = interp10d(t,pol1,temps,'nearest')./180*pi;
		freq1 = interp10d(t,freq1,temps,'nearest');
		duty1 = interp10d(t,duty1,temps,'nearest');
		[xece1,nh1,apol1] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,0.1,z0dinput.geo.a, ...
		rA1,zA1,pol1,tor1,z0dinput.geo.b0,118);
		if angvara1 == 0
			xece1(:) = xece1(1); 
			apol1(:) = apol1(1); 
			tor1(:)  = tor1(1); 
		end
		xece  = xece + xece1 .* pfce1;
		apol  = apol + abs(apol1) .* pfce1;
		ator  = ator + tor1 .* pfce1;
	end
end
if ~isempty(pfce2)
	if onoffa2 == 1
		rA2  = 3.5300;
		zA2  = 0;
		pfce2 = interp10d(t,pfce2,temps,'nearest');
		tor2 = interp10d(t,tor2,temps,'nearest')./180*pi;
		pol2 = interp10d(t,pol2,temps,'nearest')./180*pi;
		freq2 = interp10d(t,freq2,temps,'nearest');
		duty2 = interp10d(t,duty2,temps,'nearest');
		[xece2,nh2,apol2] = r2xece(z0dinput.geo.R,z0dinput.geo.z0,0.1,z0dinput.geo.a, ...
		rA2,zA2,pol2,tor2,z0dinput.geo.b0,118);
		if angvara2 == 0
			xece2(:) = xece2(1); 
			apol2(:) = apol2(1); 
			tor2(:)  = tor2(1); 
		end
		xece  = xece + xece2 .* pfce2;
		apol  = apol + abs(apol2) .* pfce2;
		ator  = ator + tor2 .* pfce2;
	end
end
pfce = interp10d(t,pfce,temps,'nearest');
z0dinput.cons.xece   = max(0,min(1,xece ./ max(eps,pfce)));
z0dinput.option.sens = sign(trapz(temps,ator)) ;
z0dinput.option.angle_ece = trapz(temps,apol) ./ max(1,trapz(temps,pfce));

% prise en compte modulation fce
pfcem = zeros(size(temps));
if ~isempty(pfce1)
	if onoffa1 == 1
		pfcem1 = 0 * pfce1;
		nexton = [];
		nextoff = [];
		fprintf('ECRH1:');
		for k =1:length(temps);
			if pfce1(k) > 0
				if duty1(k) == 1
					pfcem1(k) = 1;	
					fprintf('#');
				elseif (pfce1(k) > 0) & isempty(nexton)
					fprintf('^');
					pfcem1(k) = 1;
					nexton = temps(k) + 1./ freq1(k);
					nextoff = min(temps(k) + duty1(k)./ freq1(k),nexton -eps);
				elseif (temps(k) >= nextoff) & (pfcem1(k-1) > 0) & (k > 1)
					pfcem1(k) = 0;
					fprintf('v');
				elseif (temps(k) >= nexton) & (pfcem1(k-1) == 0) & (k > 1)
					pfcem1(k) = 1;
					fprintf('^');
					nexton = temps(k) + 1./ freq1(k);
					nextoff = min(temps(k) + duty1(k)./ freq1(k),nexton -eps);
				elseif 	(pfce1(k) > 0) & (k >1)	
					pfcem1(k) = pfcem1(k-1);
					fprintf('=');
				elseif 	(pfce1(k) >0)	 & (k == 1)
					pfcem1(k) = 1;
					fprintf('^');
				else 
					pfcem1(k) = 0;
					fprintf('v');
				end
			else
				pfcem1(k) = 0;
				fprintf('_');
				nexton = [];
				nextoff = [];
			end
		end
		fprintf('\n');	
		pfcem = pfcem1 .* pfce1 .* 1e6 + pfcem;
	end
end
if ~isempty(pfce2)
	if onoffa2 == 1
		pfcem2 = 0 * pfce2;
		nexton = [];
		nextoff = [];
		fprintf('ECRH2:');
		for k =1:length(temps);
			if pfce2(k) > 0
				if duty2(k) == 1
					pfcem2(k) = 1;	
					fprintf('#');
				elseif (pfce2(k) > 0) & isempty(nexton)
					fprintf('^');
					pfcem2(k) = 1;
					nexton = temps(k) + 1./ freq2(k);
					nextoff = min(temps(k) + duty2(k)./ freq2(k),nexton -eps);
				elseif (temps(k) >= nextoff) & (pfcem2(k-1) > 0) & (k > 1)
					pfcem2(k) = 0;
					fprintf('v');
				elseif (temps(k) >= nexton) & (pfcem2(k-1) == 0) & (k > 1)
					pfcem2(k) = 1;
					fprintf('^');
					nexton = temps(k) + 1./ freq2(k);
					nextoff = min(temps(k) + duty2(k)./ freq2(k),nexton -eps);
				elseif 	(pfce2(k) > 0) & (k >1)	
					pfcem2(k) = pfcem2(k-1);
					fprintf('=');
				elseif 	(pfce2(k) >0)	 & (k == 1)
					pfcem2(k) = 1;
					fprintf('^');
				else 
					pfcem2(k) = 0;
					fprintf('v');
				end
				
			else
				pfcem2(k) = 0;
				fprintf('_');
				nexton = [];
				nextoff = [];
			end				
			%disp([temps(k),nextoff,nexton,pfcem2(k)])

		end
		fprintf('\n');	
		pfcem = pfcem2 .* pfce2 .* 1e6 + pfcem;
	end
end
if ~all(pfcem == 0) & (chocfce > 0)
	z0dinput.cons.pecrh  = pfcem;
end
	
z0dinput.option.available_flux = 10.3 - (-10.1);  % usual 9.8  - (-7.82)
	
% choc de test 34733
if refshot > 0 
	void = zerod_init(1,refshot,((z0dinput.option.gaz == 4) + 1),z0dinput.cons.temps);
	if isfield(void.exp0d,'Rsepa')
		void.exp0d = rmfield(void.exp0d,'Rsepa');
		void.exp0d = rmfield(void.exp0d,'Zsepa');
	end
	z0dinput.exp0d = void.exp0d;
	z0dinput.cons.zeff = interp10d(void.cons.temps,void.cons.zeff,z0dinput.cons.temps,'nearest');
	ipbad = find(void.cons.ip < 3e5);
	if ~isempty(ipbad)
		[n,z]  = hist(void.cons.zeff(find(void.cons.ip > 3e5)));
		zeff0  = z(find(n==max(n),1));
		z0dinput.cons.zeff(ipbad) = zeff0;
	end
	z0dinput.option.cmin = void.option.cmin ;
	if isfield(z0dinput.exp0d,'XDURt') & isfield(z0dinput.exp0d,'XDURx') & isfield(z0dinput.exp0d,'XDURv')
		indok            = find(all(z0dinput.exp0d.XDURv>=0,2));
		if length(indok) >3
			z0dinput.option.dlh       = 0;
		end
	end
	
	% etalonage consigne sur reference cons = etalon(exp,cons_exp,cons)
	if shot ~= refshot
		etal = zerod_init(11,refshot,[],z0dinput.cons.temps);
	else 
		etal = z0dinput;
	end
	fprintf('correction nbar = ')
	% ajout de la densite naturelle
	z0dinput.cons.nbar  =  etalon(void.cons.nbar,etal.cons.nbar,z0dinput.cons.nbar);
	z0dinput.cons.nbar  = max(min(void.cons.nbar),z0dinput.cons.nbar);		       
	fprintf('correction plh = ')
	z0dinput.cons.plh   =  etalon(void.cons.plh,etal.cons.plh,z0dinput.cons.plh);
	fprintf('correction pecrh = ')
	z0dinput.cons.pecrh =  etalon(void.cons.pecrh,etal.cons.pecrh,z0dinput.cons.pecrh);
	fprintf('correction picrh = ')
	z0dinput.cons.picrh =  etalon(void.cons.picrh,etal.cons.picrh,z0dinput.cons.picrh); 
	fprintf('correction ip = ')
	z0dinput.cons.ip    =  etalon(void.cons.ip,etal.cons.ip,z0dinput.cons.ip);
	
	% on recopie la consigne de idn
	% provisoirement on recopie le choc de reference
	if any(void.cons.pnbi > 10e3) 
	    %z0dinput.cons.pnbi = interp10d(void.cons.temps,void.cons.pnbi,z0dinput.cons.temps,'nearest');
	    z0dinput.option.rtang = 1.47;
	    z0dinput.option.angle_nbi = 90;
	    z0dinput.option.einj = void.option.einj; 
	    z0dinput.cons.ftnbi = ones(size(z0dinput.cons.temps)) .* mean(void.cons.ftnbi);    
	end
end
	

%============================================================
% creation d'un jeux de donnees vide pour machine quelconque
%============================================================
function z0dinput = zerod_init_empty(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
% pas de donnees
if ispc
	user = getenv('USERNAME');
else
	user = getenv('USER');
end
prompt={'begin time (s):','end time (s):','time step (s):','Device name','Shot number'};
def={'0','30','0.1',upper(sprintf('%s_Tokamak',user)),fix(rand(1) .* 1e6 -1)};
dlgTitle='METIS time slices for new device';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);
if isempty(answer)
    mem = z0dinput;
    try
	    z0dinput = evalin('base','z0dinput');
    catch
	    try
		    z0dinput = evalin('base','post.z0dinput');
	    catch
		    z0dinput = [];
	    end
    end
    return
end
tdeb  = str2num(answer{1});
tfin  = str2num(answer{2});
pas   = str2num(answer{3});
temps = (tdeb:pas:tfin)';
vt    = ones(size(temps));
par = zerod;
z0dinput.info        = par.info;
z0dinput.option      = par.valeur;
z0dinput.zsinfo      = zero1t;
z0dinput.langue      =  lower(getappdata(0,'langue_cronos'));
z0dinput.cons.temps   = temps;
z0dinput.cons.ip      = 1e6    .* vt;
z0dinput.cons.flux    = 0      .* vt;
z0dinput.cons.nbar    = 2.5e19 .* vt;
z0dinput.cons.picrh   = 5e6    .* vt;
z0dinput.cons.plh     = 3e6    .* vt;
z0dinput.cons.pnbi    = 1e6    .* vt;
z0dinput.cons.pecrh   = 0.7e6  .* vt;
z0dinput.cons.zeff    = 3      .* vt;
z0dinput.cons.xece    = 0      .* vt;
z0dinput.option.li    = 1;
librute               = 1.7 .* vt;
z0dinput.cons.hmore   = 1      .* vt;

z0dinput.geo.a        = 0.72   .* vt;
z0dinput.geo.R        = 2.4    .* vt; 
z0dinput.geo.K        = 1      .* vt;  
z0dinput.geo.d        = 0      .* vt;
z0dinput.geo.z0        = 0      .* vt;
z0dinput.geo.b0       = 3.8    .* vt;
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];
z0dinput.mode_exp     = -1;
z0dinput.exp          = [];
z0dinput.machine     = answer{4};
z0dinput.shot        = fix(str2num(answer{5}));


% defaut iter pour ce cas
z0dinput.option.etalh     = 0.8;
z0dinput.option.lhmode    = 3;
z0dinput.option.xlh        = 0.7;
z0dinput.option.dlh        = 0.4;

%============================================================
% creation d'un jeux de donnees pour metis evolution
%============================================================
function z0dinput = zerod_init_evolution(mode_exp,shot,gaz,temps,z0dinput)
		
langue                 =  'anglais';
% pour evolution
vt    = ones(size(temps));
par = zerod;
z0dinput.info        = par.info;
z0dinput.option      = par.valeur;
z0dinput.zsinfo      = zero1t;
z0dinput.langue      =  lower(getappdata(0,'langue_cronos'));
z0dinput.cons.temps   = temps;
z0dinput.cons.ip      = 1e6    .* vt;
z0dinput.cons.flux    = 0      .* vt;
z0dinput.cons.nbar    = 2.5e19 .* vt;
z0dinput.cons.picrh   = 5e6    .* vt;
z0dinput.cons.plh     = 3e6    .* vt;
z0dinput.cons.pnbi    = 1e6    .* vt;
z0dinput.cons.pecrh   = 0.7e6  .* vt;
z0dinput.cons.zeff    = 3      .* vt;
z0dinput.cons.xece    = 0      .* vt;
z0dinput.option.li    = 1;
librute               = 1.7 .* vt;
z0dinput.cons.hmore   = 1      .* vt;

z0dinput.geo.a        = 0.72   .* vt;
z0dinput.geo.R        = 2.4    .* vt; 
z0dinput.geo.K        = 1      .* vt;  
z0dinput.geo.d        = 0      .* vt;
z0dinput.geo.z0        = 0      .* vt;
z0dinput.geo.b0       = 3.8    .* vt;
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];
z0dinput.mode_exp     = -1;
z0dinput.exp          = [];
if evalin('base','exist(''z0dinput'')')
    refinput = evalin('base','z0dinput');
    if ~isempty(refinput)
        z0dinput.machine     = refinput.machine;
        z0dinput.shot        = refinput.shot;
    else
        z0dinput.machine     = 'none';
        z0dinput.shot        = NaN;
    end
else
    z0dinput.machine     = 'none';
    z0dinput.shot        = NaN;
end

% defaut iter pour ce cas
z0dinput.option.etalh     = 0.8;
z0dinput.option.lhmode    = 3;
z0dinput.option.xlh        = 0.7;
z0dinput.option.dlh        = 0.4;


%============================================================
% lecture des donnees dans l'UAL
%============================================================
function z0dinput = zerod_init_ual(mode_exp,shot,gaz,temps,z0dinput)

langue                 =  'anglais';
%prompt={'shot number :','Run :' ,'Occurrence:',};
prompt={'shot number :','Run :' ,'Occurrence:','Tokamak:','User:','Version:'};
def={'1','1','','','',''};
dlgTitle='Readig data from UAL';
lineNo=1;
answer=zinputdlg(prompt,dlgTitle,lineNo,def);
if isempty(answer)
    z0dinput = [];
    return
end
shot  = str2num(answer{1});
run   = str2num(answer{2});
occ   = answer{3};
ual_tok   = answer{4};
ual_user   = answer{5};
ual_ver    = answer{6};
if ~isempty(ual_tok) && ~isempty(ual_tok(ual_tok > ' '))
	ual_occ = occ;
        clear occ
        occ.tokamak = ual_tok(ual_tok>' ');
        occ.user    = ual_user(ual_user> ' ');
        occ.dataversion = ual_ver(ual_ver> ' ');
        occ.occurrence = ual_occ(ual_occ> ' ');
end

[error_flag,z0dinput] = metis4itm(shot,run,occ,'metis_from_ual');
z0dinput.shot =  z0dinput.shot(1);
z0dinput.mode_exp = 117;
z0dinput.run = run;
if isstruct(occ)
   	z0dinput.tokamak = occ.tokamak;
    	z0dinput.user = occ.user;
    	z0dinput.dataversion = occ.dataversion;
    	z0dinput.occ = occ.occurrence;
else
	z0dinput.occ = occ;
end
z0dinput.shot = shot;
    

function out = chg2cell(in)

out = {};
if isempty(in)
	return
end
in(in <= ' ') = [];
if isempty(in)
	return
end
while ~isempty(in)
	[f,in] = strtok(in,',');
	out{end+1} = f;	
end


function z0dinput = zerod_init_asdex_u(mode_exp,shot,gaz,temps,z0dinput)

disp('Not yet implanted !');
z0dinput = [];
		
