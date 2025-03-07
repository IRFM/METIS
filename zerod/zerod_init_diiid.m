%============================================================
% acces au donnees de DIIID via MSD+
%============================================================
function z0dinput = zerod_init_diiid(mode_exp,shot,gaz,temps,z0dinput)

listenbi ={'30l','30r','15l','15r','21l','21r','33l','33r'};

langue                 =  'anglais';
% cas donnees DIIID
z0dinput.exp          = [];
% 1 - numero du choc
if isempty(shot) || isempty(gaz)
    prompt={'shot number :','charge of main gas:','number of time slices multiplier','NBI1 list','NBI2 list','NBI in ICRH list'};
    % def={'147634','1','1','15l,15r','30l,30r,33l,33r','21l,21r'};
    def={'147634','1','1','30l,30r,33l,33r','15l,15r','21l,21r'};
    % def={'147634','1','1','30l,30r,15l,15r,33l,33r','','21l,21r'};
    dlgTitle='Access to DIIID data';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
	disp('Operation cancelled');
	z0dinput = [];
	return
    end
    shot  = str2num(answer{1});
    gaz  = str2num(answer{2});
    nbtfact  = str2num(answer{3});
    list_nbi1  = chg2cell(answer{4});
    list_nbi2  = chg2cell(answer{5});
    list_icrh  = chg2cell(answer{6});

end
%
% lecture des donnees
%
home       = getenv('HOME');
racine     = fullfile(home,'zineb/data/diiid');
directory  = fullfile(racine,int2str(shot));
filetemp   = [directory,'/diiidtemp.mat'];
fileprof   = [directory,'/diiidprof.mat'];
filediag   = [directory,'/diiiddiag.mat'];
filetransp = [directory,'/diiidtransp.mat'];
fileonetwo = [directory,'/diiidonetwo.mat'];
filebord   = [directory,'/diiidbord.mat'];
%
load(filetemp)
load(fileprof)
if exist(filebord,'file');
  load(filebord)
  if isfield(diiidbord,'xbord')
    disp([' the LSMF comes from diiidbord, with x = ',num2str(diiidbord.xbord)])
  else
    disp([' the LSMF comes from diiidbord '])
  end
else
  diiidbord = [];
end
if exist(filediag,'file');
  load(filediag)
else
  diiiddiag = [];
end
if exist(fileonetwo,'file');
  load(fileonetwo)
else
  diiidonetwo = [];
end

if ~isempty(diiidtemp.ip)
	tip   = diiidtemp.tip;
	ip    = diiidtemp.ip;
	ip    = ip(tip >=0);
	temps = tip(tip >=0);
	if nbtfact > 1
		temps = linspace(min(temps),max(temps),length(temps) .* nbtfact)';
		ip    = interp1(diiidtemp.tip,diiidtemp.ip,temps,'linear');
	end
else
	disp('No current plasma measurement available')
	return
end

% Plasma current : positive if counter-clockwise viewed from above the torus
signe_ip               = sign(mean(ip));
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = abs(ip);

% DM 04/03/2013 remove 2*pi and next 3 lines
z0dinput.cons.flux     = signe_ip .*  interp1(diiidtemp.t,-diiidtemp.psia,temps,'linear'); 
% if z0dinput.cons.flux(end) > z0dinput.cons.flux(1)
% 	z0dinput.cons.flux = - z0dinput.cons.flux;
% end
z0dinput.geo.a         = interp1(diiidtemp.t,diiidtemp.a,temps);
z0dinput.geo.z0        = 0*temps;
z0dinput.geo.R         = interp1(diiidtemp.t,diiidtemp.R0,temps);
z0dinput.geo.K         = interp1(diiidtemp.t,diiidtemp.e,temps);
z0dinput.geo.d         = 0*temps;
z0dinput.geo.vp        = 0*temps;
z0dinput.geo.sp        = 0*temps;
z0dinput.geo.sext      = 0*temps;

% Toroidal field : positive if counter-clockwise viewed from above the torus
signe_B0               = sign(mean(diiidtemp.B0));
z0dinput.geo.b0        = interp1(diiidtemp.t,abs(diiidtemp.B0),temps);

if isempty(diiidbord)
  indok = find(all(isfinite(diiidtemp.rext) & isfinite(diiidtemp.zext),2));
  if length(indok) > 3
    % Correction DM du 18/05/16
  	% z0dinput.exp0d.Rsepa  = interp1(diiidtemp.trext(indok),diiidtemp.rext(indok,:),temps);
 	% z0dinput.exp0d.Zsepa  = interp1(diiidtemp.trext(indok),diiidtemp.zext(indok,:),temps);
    disp('Reading separatrix from diiidtemp');
  	z0dinput.exp0d.Rsepa  = interp1(diiidtemp.trext(indok),diiidtemp.rext(indok,:),temps,'nearest','extrap');
 	z0dinput.exp0d.Zsepa  = interp1(diiidtemp.trext(indok),diiidtemp.zext(indok,:),temps,'nearest','extrap');
  else
	fprintf('too many NaN in Rext & Z ext, separatrix data not usable\n'); 
  end
else
  % Correction DM du 18/05/16
  % z0dinput.exp0d.Rsepa   = interp1(diiidbord.t,diiidbord.Rext,temps);
  % z0dinput.exp0d.Zsepa   = interp1(diiidbord.t,diiidbord.Zext,temps);
  disp('Reading separatrix from diiidbord');
  z0dinput.exp0d.Rsepa   = interp1(diiidbord.t,diiidbord.Rext,temps,'nearest','extrap');
  z0dinput.exp0d.Zsepa   = interp1(diiidbord.t,diiidbord.Zext,temps,'nearest','extrap');
end
if ~isempty(z0dinput.exp0d.Rsepa)
	Rext = z0dinput.exp0d.Rsepa;
	Zext = z0dinput.exp0d.Zsepa;
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
	rb0                  = z0dinput.geo.R .* z0dinput.geo.b0;
	z0dinput.geo.R       = ra;
	z0dinput.geo.b0      = rb0 ./ z0dinput.geo.R;
	z0dinput.geo.z0      = za;
	z0dinput.geo.K       = k;
	z0dinput.geo.d       = d;
end

% z0dinput.cons.nbar     = max(1e13,trapz(diiidprof.x',interp1(diiidprof.t,diiidprof.ne,temps),2) ./ trapz(diiidprof.x',ones(size(diiidprof.x')),2));
% DM 04/03/2013 : remplacer par nbar from EFIT
% DM 15/04/2014 : subtract beam fueling from nbar (use diiidtemp.dnbar_nb ~= 0 and diiidtemp.reycling ~0 when neasser = 1)
% Parameters: Density
if ~isfield(diiidtemp,'dnbar_nb')
    z0dinput.cons.nbar = interp1(diiidtemp.tnbar,diiidtemp.nbar,temps,'linear');
    z0dinput.option.neasser   = 0;
    z0dinput.option.Recycling = 0;
else
    % Source term in the ODE (see zero1t.m)
    z0dinput.cons.nbar = interp1(diiidtemp.tnbar,diiidtemp.nbar + diiidtemp.dnbar_nb,temps,'linear');
    % Initial value in the ODE (see zero1t.m)
    z0dinput.cons.nbar(1) = interp1(diiidtemp.tnbar,diiidtemp.nbar,temps(1),'linear');
    z0dinput.option.neasser   = 1;
    z0dinput.option.Recycling = diiidtemp.recycling;
end
% DM 25/09/14 : New optional parameter governing n0a
z0dinput.option.fn0a_div   = 1;

% Puissance ICRH (peut etre utilisee plus bas pour simuler NBI sans generation de courant/rotation (balanced injection)
z0dinput.cons.picrh   = 0 .* temps;

% Possibilite d'utiliser LH pour simuler un 2eme ECCD avec largeur de depot variable si gyrinlh non vide
% DM 04/03/2013 : filtrer les modulations rapides avant d'interpoler
% kc  = 2*round(mean(diff(temps))/mean(diff(diiidtemp.tpec)))+1;
kc  = max(3,fix((length(diiidtemp.tpec) ./ length(temps) + 1) / 2) .* 2 + 1);
% gyrinlh=1:size(diiidtemp.pec,2);
gyrinlh=[];
gyrinec=setdiff(1:size(diiidtemp.pec,2),gyrinlh);
if ~isempty(gyrinlh)
    % defaut pour utiliser LH pour simuler un 2eme ECCD 
    z0dinput.option.lhmode = 5;
    z0dinput.option.etalh  = 1;
    z0dinput.option.xlh    = 0.5;
    z0dinput.option.dlh    = 0.15;
    z0dinput.cons.plh   = sum(interp1(diiidtemp.tpec,sgolayfilt(diiidtemp.pec(:,gyrinlh),1,kc),temps,'linear'),2);
    z0dinput.cons.pecrh = sum(interp1(diiidtemp.tpec,sgolayfilt(diiidtemp.pec(:,gyrinec),1,kc),temps,'linear'),2);
else
    z0dinput.option.lhmode = 2;
    z0dinput.option.etalh  = 1e19;
    z0dinput.option.xlh    = 0.5;
    z0dinput.option.dlh    = 0.15;
    z0dinput.cons.plh      = 0 .* temps;
    z0dinput.cons.pecrh    = sum(interp1(diiidtemp.tpec,sgolayfilt(diiidtemp.pec,1,kc),temps,'linear'),2);
end
%figure(21);clf
%plot(temps,z0dinput.cons.plh,'r',diiidtemp.tpec,sgolayfilt(diiidtemp.pec,1,kc),'b',diiidtemp.tpec,diiidtemp.pec,'g')
%plot(temps,z0dinput.cons.pecrh,'r',diiidtemp.tpec,sgolayfilt(diiidtemp.pec,1,kc),'b',diiidtemp.tpec,diiidtemp.pec,'g')
%drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifs DM 2014
z0dinput.cons.xece          = 0.5 .* ones(size(temps));
z0dinput.option.sens        = 1;
z0dinput.option.eccdmul     = 1;
z0dinput.option.synergie    = 1; % means no synergie
%
% Composition (zeff profile given by external data)
z0dinput.option.zeff = 0;

% Default zeff reference
if ~isempty(diiidtemp.zeffm)
    z0dinput.cons.zeff = interp1(diiidtemp.tzeffm,diiidtemp.zeffm,temps,'linear');
else
    z0dinput.cons.zeff = 2 .* ones(size(temps));
end
% Compute zeff reference from zeff profiles in onetwo simulations (only when the onetwo profile is not flat)
k12zf = find(diiidonetwo.kmanual & max(diiidonetwo.zeff,[],2)-min(diiidonetwo.zeff,[],2) > 0);
if ~isempty(k12zf)
    % space interpolation on the diiidtemp rho grid
    zeffp = interp1(diiidonetwo.x,diiidonetwo.zeff(k12zf,:)',diiidtemp.x)';
    % time interpolation on the diiidonetwo time grid
    xlao=interp1(diiidtemp.t,diiidtemp.xlao,diiidonetwo.t(k12zf),'linear');
    % space interpolation on a fixed metis grid and line-average using metis coordinates (Rmax - Rmin)/ 2a
    xmetis=linspace(0,1,21);
    zeffl=zeros(length(k12zf),21);
    for kk=1:length(k12zf)
        zeffl(kk,:)=interp1(xlao(kk,:),zeffp(kk,:)',xmetis)';
    end
    zeffl=trapz(xmetis,zeffl,2);
    % time interpolation
    z0dinput.cons.zeff(temps >= diiidonetwo.t(k12zf(1)) & temps <= diiidonetwo.t(k12zf(end))) = ...
        interp1([-1;diiidonetwo.t(k12zf);max(diiidonetwo.t(k12zf(end)),temps(end))+1],[zeffl(1);zeffl;zeffl(end)], ...
        temps(temps >= diiidonetwo.t(k12zf(1)) & temps <= diiidonetwo.t(k12zf(end))),'linear');
end
% figure(5);clf;
% plot(z0dinput.cons.temps,z0dinput.cons.zeff,'bo')
% if isfield(diiidtemp,'zf_factor'), z0dinput.cons.zeff = diiidtemp.zf_factor.*z0dinput.cons.zeff; end;
% hold on;plot(z0dinput.cons.temps,z0dinput.cons.zeff,'ro')
% z0dinput.cons.zeff(z0dinput.cons.zeff < 1) = 1;
% plot(z0dinput.cons.temps,z0dinput.cons.zeff,'r')

% Tension par tour et li initiaux (Breakdown & burnthrough) avec runaways (Current diffusion and Equilibrium).
% option.breakdown = electric field at starting time of the simulation in unit of Dreicer electric field if > 0
% or in Volt per turn if < 0 (used abs value, convert internally into electric field)
vsurf=interp1(diiidtemp.tvsurf,diiidtemp.vsurf,temps,'linear');
z0dinput.option.breakdown = -vsurf(1);

li = interp1(diiidtemp.tli,diiidtemp.li,temps,'linear');
z0dinput.option.li  = li(1);

z0dinput.option.runaway = 1;

% Neoclassical: facteur multiplicatif pour le courant de bootstrap
z0dinput.option.modeboot = 1; % Sauter
if ~isfield(diiidtemp,'bs_factor')
    z0dinput.option.bootmul = 1;
else
    z0dinput.option.bootmul = diiidtemp.bs_factor;
end

% Current diffusion and equilibrium: Current control or flux control and new regularisation du profil de q au centre
z0dinput.option.vloop = 0;
% z0dinput.option.vloop = 4;
z0dinput.option.tswitch = Inf;
z0dinput.option.cronos_regul = 0;

% Confinement and Transport
if ~isfield(diiidtemp,'H_factor')
    z0dinput.cons.hmore = ones(size(temps));
else
    z0dinput.cons.hmore = interp1(diiidtemp.t,diiidtemp.H_factor,temps,'linear');
end
z0dinput.option.scaling  = 5;
z0dinput.option.fpped    = 1.5;
z0dinput.option.xiioxie  = 0.454;
z0dinput.option.kishape  = 3.00;
z0dinput.option.ki_expo  = 4.01;
z0dinput.option.xieorkie = 1;
z0dinput.option.grad_ped = 1;
z0dinput.option.omega_shape = 0;
z0dinput.option.taurotmul   = 2.14;
z0dinput.option.fintrinsic  = 0.2;
% SOL
z0dinput.option.lambda_scale=0;
% NBI/NBICD
z0dinput.option.cur_nbi_time = 1;
z0dinput.option.e_shielding  = 'Honda-Sauter';

%   End of Modifs DM 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% activation des dents de scies	
z0dinput.option.qdds         = 1;
z0dinput.option.kidds        = 3;


z0dinput.option.zmax      = 18;
z0dinput.option.zimp      = 6;
z0dinput.option.rimp      = 0.1;
z0dinput.option.modeh     = 1;
z0dinput.cons.iso         =  zeros(size(temps));
% if gaz  == 1
% 	z0dinput.option.gaz       =  2;
% else
% 	z0dinput.option.gaz       =  4;
% end
switch gaz 
    case 1
        z0dinput.option.gaz =  2;
    case {5,11}
        error('Not yet implemented');
    otherwise
        z0dinput.option.gaz =  4;
end

% DM 04/03/2013 no need for pnbi and pidn
% if sum(diiidtemp.pidn(:)) < sum(diiidtemp.pinj(:))
%    diiidtemp.pidn = sum(diiidtemp.pinj,2);
% end
% if sum(diiidtemp.pidn(:)) > 0 
if sum(diiidtemp.pinj(:)) > 0 

     % MS 17/10/08:
     % NEW CALCULATION OF PARTICLE ENERGY FRACTIONS: AFTER NEUTRALIZER!!!
     % DM 04/03/2013 : filtrer les modulations rapides avant d'interpoler
     % kc           = 2*round(mean(diff(temps))/mean(diff(diiidtemp.tpidn)))+1;
     kc             = max(3,fix((length(diiidtemp.tpidn) ./ length(temps) + 1) / 2) .* 2 + 1);
     pidn_vect      = interp1(diiidtemp.tpidn,sgolayfilt(diiidtemp.pinj,1,kc),temps,'linear');
     Pright         = diiidtemp.pinj(:,2:2:end);
     Pleft          = diiidtemp.pinj(:,1:2:end);
     Ptotr          = sum(Pright);
     Ptotl          = sum(Pleft);
     indlr          = find((Ptotr>0+Ptotl>0)>0);
     if ~isempty(indlr)
       for j=1:length(indlr)
         disp([' correction factor to subtract 3% from power when two beams (l and r) '])
         indb = 2*indlr(j)-1;
         pidn_vect(:,indb:(indb+1))=pidn_vect(:,indb:(indb+1))*0.97;
       end
     end
else
     pidn_vect            = zeros(length(temps),8);
end
pidn_vect(~isfinite(pidn_vect)) = 0;

% MS 02.09.08 new variable "directivity" added, due to the two diiid configurations
% (ie with 2 sources becoming counter-current after #123419)
directivity = ones(1,8);
if shot > 123419
  directivity(5:6) = -1;
end

   
% TANGENCY RADIUS PER PINI (DIFFERENT FOR RIGHT AND LEFT SOURCES) (M)
% ORDER: 30L 30R 150L 150R 210L 210R 330L 330R
rtang_source([1 3 7]) = 114.6*1.e-2; % 30L 150L 330L (LEFT SOURCES)  (M)
rtang_source([2 4 8]) = 76.2*1.e-2;  % 30R 150R 330R (RIGHT SOURCES) (M)
if directivity(5)==1 % CO-INJECTION (shot <= 123419)
   rtang_source(5) = 114.6*1.e-2; % 210L
   rtang_source(6) = 76.2*1.e-2;  % 210R
else                 % COUNTER-INJECTION (shot > 123419)
   rtang_source(5) = 76.2*1.e-2;  % 210L
   rtang_source(6) = 114.6*1.e-2; % 210R
end
angleliste = [1,1,1,1,-1,-1,1,1] * 90;

% listenbi ={'30l','30r','15l','15r','21l','21r','33l','33r'};

[void,inbi]= intersect(listenbi,list_nbi1);
fprintf('injector 1 :');fprintf('%s ',listenbi{inbi});fprintf('\n');
if ~isempty(inbi)
	pnbi1   = sum(pidn_vect(:,inbi),2);
	einj1   = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * diiidtemp.Einj(inbi)*1e3),2),1) ./ max(1,trapz(temps,pnbi1,1));
	rtang1  = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * rtang_source(inbi)),2),1) ./ max(1,trapz(temps,pnbi1,1));
	angle_nbi1 = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * angleliste(inbi)),2),1) ./ max(1,trapz(temps,pnbi1,1));
	% zext1 ~=0 if Pnbi1 is used for off-axis beams (depends on tilt angle of 150R/L injectors)
	% zext1   = 0.4; % normalized radius
	zext1   = 0; % normalized radius
 	if any(pnbi1 > 0)
		fnbi1   = 1;
	else
		fnbi1   = 0;
	end
else
	pnbi1  = 0 * temps;
	rtang1 = 1.146;
	angle_nbi1 = 90;
	zext1  = 0;
	fnbi1  = 0;
end
[void,inbi]= intersect(listenbi,list_nbi2);
fprintf('injector 2 :');fprintf('%s ',listenbi{inbi});fprintf('\n');
if ~isempty(inbi)
	pnbi2   = sum(pidn_vect(:,inbi),2);
	einj2   = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * diiidtemp.Einj(inbi)*1e3),2),1) ./ max(1,trapz(temps,pnbi2,1));
	rtang2  = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * rtang_source(inbi)),2),1) ./ max(1,trapz(temps,pnbi2,1));
	angle_nbi2 = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * angleliste(inbi)),2),1) ./ max(1,trapz(temps,pnbi2,1));
    % zext2   = 0.4; % normalized radius
	% zext2 ~=0 if Pnbi2 is used for off-axis beams (depends on tilt angle of 150R/L injectors)
    zext2   = 0; % normalized radius
 	if any(pnbi2 > 0)
		fnbi2   = 1;
	else
		fnbi2   = 0;
	end
else
	pnbi2  = 0 * temps;
	rtang2 = 1.146;
	angle_nbi2 = 90;
	zext2  = 0;
	fnbi2  = 0;
end

z0dinput.option.nb_nbi     = 1;
z0dinput.option.einj2      = 1e5;
z0dinput.option.rtang2     = 2;		
z0dinput.option.angle_nbi2 = 90;		
z0dinput.option.zext2      = 0;		


[void,inbi]= intersect(listenbi,list_icrh);
fprintf('icrh use as an NBI :');fprintf('%s ',listenbi{inbi});fprintf('\n');
if ~isempty(inbi)
	pnbi3   = sum(pidn_vect(:,inbi),2);
	einj3   = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * diiidtemp.Einj(inbi)*1e3),2),1) ./  max(1,trapz(temps,pnbi3,1));
	rtang3  = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * rtang_source(inbi)),2),1) ./  max(1,trapz(temps,pnbi3,1));
	angle_nbi3 = trapz(temps,sum(pidn_vect(:,inbi) .* (ones(size(temps)) * angleliste(inbi)),2),1) ./  max(1,trapz(temps,pnbi3,1));
	zext3   = 0;
 	if any(pnbi3 > 0)
		fnbi3   = 1;
	else
		fnbi3   = 0;
	end
else
	pnbi3  = 0 * temps;
	rtang3 = 1.146;
	angle_nbi3 = 90;
	zext3  = 0;
	fnbi3   = 0;
        einj3   =1e5;
end
if (fnbi1 == 1) && (fnbi2 == 1) && (fnbi3 == 1)
	z0dinput.option.einj       = einj1;
	z0dinput.option.rtang      = rtang1;		
	z0dinput.option.angle_nbi  = angle_nbi1;		
	z0dinput.option.zext       = zext1;
	z0dinput.option.nb_nbi     = 2;
	z0dinput.option.einj2      = einj2;
	z0dinput.option.rtang2     = rtang2;		
	z0dinput.option.angle_nbi2 = angle_nbi2;		
	z0dinput.option.zext2      = zext2;		
	z0dinput.cons.pnbi         = pnbi1 + sqrt(-1) .* pnbi2;

	% icrh pour simuler le 3ieme injecteur
	if rtang3 > 0
		% position de la resonance
		rres = rtang3 - zext3 .* mean(z0dinput.geo.a);
		% champs a la resonnance
		bres = mean(z0dinput.geo.R .* z0dinput.geo.b0) ./ rres;
		% minoritaire hydrogene
		zg = 1;
		ag =1;
		z0dinput.option.mino = 'H';
		z0dinput.option.fwcd = 0;
		% frequence en MHz
		z0dinput.option.freq = bres ./(2.* pi .* 1e6) .* (95.5e6 .*zg ./ ag);
		% la largeur est fixe par nphi (quasi proportionnel)
		z0dinput.option.nphi = 25;
		% la balance entre chauffage aux ions et aux electrons est fixe par la fraction de minoritaire:
		% cmin -> 0 = chauffage electronique
		% cmin ~ 0.1 = chauffage majoritaire au ions
		z0dinput.option.cmin = 0.1;
		% puissance
		z0dinput.cons.picrh = pnbi3;
	end
elseif (fnbi1 == 1) && (fnbi2 == 1)
	z0dinput.option.einj       = einj1;
	z0dinput.option.rtang      = rtang1;		
	z0dinput.option.angle_nbi  = angle_nbi1;		
	z0dinput.option.zext       = zext1;
	z0dinput.option.nb_nbi     = 2;
	z0dinput.option.einj2      = einj2;
	z0dinput.option.rtang2     = rtang2;		
	z0dinput.option.angle_nbi2 = angle_nbi2;		
	z0dinput.option.zext2      = zext2;		
	z0dinput.cons.pnbi         = pnbi1 + sqrt(-1) .* pnbi2;

elseif (fnbi1 == 1) && (fnbi3 == 1)
	z0dinput.option.einj       = einj1;
	z0dinput.option.rtang      = rtang1;		
	z0dinput.option.angle_nbi  = angle_nbi1;		
	z0dinput.option.zext       = zext1;
	z0dinput.option.nb_nbi     = 2;
	z0dinput.option.einj2      = einj3;
	z0dinput.option.rtang2     = rtang3;		
	z0dinput.option.angle_nbi2 = angle_nbi3;		
	z0dinput.option.zext2      = zext3;		
	z0dinput.cons.pnbi         = pnbi1 + sqrt(-1) .* pnbi3;

elseif (fnbi2 == 1) && (fnbi3 == 1)

	z0dinput.option.einj       = einj2;
	z0dinput.option.rtang      = rtang2;		
	z0dinput.option.angle_nbi  = angle_nbi2;		
	z0dinput.option.zext       = zext2;
	z0dinput.option.nb_nbi     = 2;
	z0dinput.option.einj2      = einj3;
	z0dinput.option.rtang2     = rtang3;		
	z0dinput.option.angle_nbi2 = angle_nbi3;		
	z0dinput.option.zext2      = zext3;		
	z0dinput.cons.pnbi         = pnbi1 + sqrt(-1) .* pnbi3;

elseif (fnbi1 == 1)
	z0dinput.option.einj       = einj1;
	z0dinput.option.rtang      = rtang1;		
	z0dinput.option.angle_nbi  = angle_nbi1;		
	z0dinput.option.zext       = zext1;
	z0dinput.option.nb_nbi     = 1;
    z0dinput.cons.pnbi         = pnbi1;

elseif (fnbi2 == 1)
	z0dinput.option.einj       = einj2;
	z0dinput.option.rtang      = rtang2;		
	z0dinput.option.angle_nbi  = angle_nbi2;		
	z0dinput.option.zext       = zext2;
	z0dinput.option.nb_nbi     = 1;
    z0dinput.cons.pnbi         = pnbi2;

elseif (fnbi2 == 3)
	z0dinput.option.einj       = einj3;
	z0dinput.option.rtang      = rtang3;		
	z0dinput.option.angle_nbi  = angle_nbi3;		
	z0dinput.option.zext       = zext3;
	z0dinput.option.nb_nbi     = 1;
    z0dinput.cons.pnbi         = pnbi3;

else
    z0dinput.cons.pnbi         = 0 * temps;
end
%
z0dinput.option.einj      = max(1,z0dinput.option.einj);
z0dinput.option.einj2     = max(1,z0dinput.option.einj2);

z0dinput.cons.ftnbi       = 0 .* z0dinput.cons.iso;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0dinput.option.ane       = 0;

% estimation grossière du piquage moyen avec les profils de densité.
neav=sum(diiidprof.ne.*(ones(length(diiidprof.t),1)*[0 diff(diiidprof.x.*diiidprof.x)']),2);
z0dinput.option.vane      = mean(diiidprof.ne(:,1)./neav);
% il faut affiner ce parametre

z0dinput.machine     = 'DIIID';
z0dinput.shot        = shot;
z0dinput.option.configuration = 2;

% signe B-T * Ip	
z0dinput.option.signe        = signe_B0*signe_ip;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% donnees experimentale pour comparaison
% DM 06/03/2013

z0dinput.exp0d.temps  = temps;
z0dinput.exp0d.pnbi   = [real(z0dinput.cons.pnbi) imag(z0dinput.cons.pnbi)];
z0dinput.exp0d.pecrh  = z0dinput.cons.pecrh;
z0dinput.exp0d.picrh  = z0dinput.cons.picrh;
z0dinput.exp0d.plh    = z0dinput.cons.plh;
z0dinput.exp0d.flux   = z0dinput.cons.flux;
%z0dinput.exp0d.pin   = ;

z0dinput.exp0d.zeff   = z0dinput.cons.zeff;
%z0dinput.exp0d.vp     = ;
%z0dinput.exp0d.sp     = ;
%z0dinput.exp0d.sext   = ;
%z0dinput.exp0d.nem    = ;
%z0dinput.exp0d.ne0    = interp1(diiidprof.t,diiidprof.ne(:,1),temps,'linear');
%z0dinput.exp0d.nhem   =  ;
%z0dinput.exp0d.nimpm  = ;
z0dinput.exp0d.nTm    = 0.* temps;
%z0dinput.exp0d.n1m    = ;
%z0dinput.exp0d.nim    = ;
%z0dinput.exp0d.ni0    = ;

%z0dinput.exp0d.te0    = interp1(diiidprof.t,diiidprof.Te(:,1),temps,'linear');
%z0dinput.exp0d.tem    = ;
%z0dinput.exp0d.tite   = ;
%z0dinput.exp0d.ate    = ;
%z0dinput.exp0d.ape    = ; 
  
% z0dinput.exp0d.nebord = interp1(diiidtemp.tnesep,diiidtemp.nesep,temps,'linear');
% z0dinput.exp0d.tebord = interp1(diiidtemp.ttesep,diiidtemp.tesep,temps,'linear');
z0dinput.exp0d.neped  = interp1(diiidtemp.tneped,diiidtemp.neped,temps,'linear');
z0dinput.exp0d.teped  = interp1(diiidtemp.tteped,diiidtemp.teped,temps,'linear');
z0dinput.exp0d.wdia   = interp1(diiidtemp.twdia,diiidtemp.wdia,temps,'linear');
z0dinput.exp0d.wbp    = interp1(diiidtemp.twbp,diiidtemp.wbp,temps,'linear');
%z0dinput.exp0d.w      = interp1(diiidtemp.tw,diiidtemp.w,temps,'linear');

k12                   = find(diiidonetwo.kmanual);
z0dinput.exp0d.ne0    = interp1(diiidonetwo.t(k12),diiidonetwo.ne(k12,1),temps,'linear');
z0dinput.exp0d.te0    = interp1(diiidonetwo.t(k12),diiidonetwo.te(k12,1),temps,'linear');
z0dinput.exp0d.ti0    = interp1(diiidonetwo.t(k12),diiidonetwo.ti(k12,1),temps,'linear');
z0dinput.exp0d.nebord = interp1(diiidonetwo.t(k12),diiidonetwo.ne(k12,end),temps,'linear');
z0dinput.exp0d.tebord = interp1(diiidonetwo.t(k12),diiidonetwo.te(k12,end),temps,'linear');
z0dinput.exp0d.w      = interp1(diiidonetwo.t(k12),diiidonetwo.w(k12),temps,'linear');
z0dinput.exp0d.dwdt   = pdederive(temps,z0dinput.exp0d.w,2,2,1,1); 
z0dinput.exp0d.wth    = interp1(diiidonetwo.t(k12),diiidonetwo.wth(k12),temps,'linear');
z0dinput.exp0d.dwthdt = pdederive(temps,z0dinput.exp0d.wth,2,2,1,1); 
z0dinput.exp0d.esup_nbi = interp1(diiidonetwo.t(k12),diiidonetwo.wsupnb(k12),temps,'linear');
z0dinput.exp0d.taue   = interp1(diiidonetwo.t(k12),diiidonetwo.taueth(k12),temps,'linear'); 
%z0dinput.exp0d.tauhe  = ;

z0dinput.exp0d.ip     = ip;
z0dinput.exp0d.iboot  = interp1(diiidonetwo.t(k12),diiidonetwo.iboot(k12),temps,'linear');
z0dinput.exp0d.inbicd = interp1(diiidonetwo.t(k12),diiidonetwo.inbicd(k12),temps,'linear');
z0dinput.exp0d.ieccd  = interp1(diiidonetwo.t(k12),diiidonetwo.ieccd(k12),temps,'linear');
z0dinput.exp0d.icd    = interp1(diiidonetwo.t(k12),diiidonetwo.icd(k12),temps,'linear');
z0dinput.exp0d.ini    = interp1(diiidonetwo.t(k12),diiidonetwo.ini(k12),temps,'linear');
z0dinput.exp0d.iohm   = z0dinput.exp0d.ip - z0dinput.exp0d.ini;

z0dinput.exp0d.betan  = interp1(diiidonetwo.t(k12),diiidonetwo.betan(k12),temps,'linear');
z0dinput.exp0d.betaptot  = interp1(diiidonetwo.t(k12),diiidonetwo.betaptot(k12),temps,'linear');
z0dinput.exp0d.pin    = interp1(diiidonetwo.t(k12),diiidonetwo.pin(k12),temps,'linear');
z0dinput.exp0d.pohm   = interp1(diiidonetwo.t(k12),diiidonetwo.pohm(k12),temps,'linear');
% z0dinput.exp0d.betan  = interp1(diiidtemp.tbetan,diiidtemp.betan,temps,'linear');
% z0dinput.exp0d.betaptot  = interp1(diiidtemp.tbetap,diiidtemp.betap,temps,'linear');
% z0dinput.exp0d.pohm   = interp1(diiidtemp.tpoh,diiidtemp.poh,temps,'linear');
z0dinput.exp0d.vloop  = vsurf;
z0dinput.exp0d.vmes   = interp1(diiidtemp.tvmes,diiidtemp.vmes,temps,'linear');
%z0dinput.exp0d.qa     = ;
z0dinput.exp0d.q95    = interp1(diiidtemp.tq95,diiidtemp.q95,temps,'linear');
z0dinput.exp0d.q0     = interp1(diiidtemp.tq0,diiidtemp.q0,temps,'linear');
z0dinput.exp0d.qmin   = interp1(diiidtemp.tqmin,diiidtemp.qmin,temps,'linear');

% z0dinput.exp0d.nbar  = z0dinput.cons.nbar;
z0dinput.exp0d.nbar  = interp1(diiidtemp.tnbar,diiidtemp.nbar,temps,'linear');
z0dinput.exp0d.li    = li;

%z0dinput.exp0d.modeh =; 
z0dinput.exp0d.ndd        = interp1(diiidonetwo.t(k12),diiidonetwo.ndd(k12),temps,'linear');
z0dinput.exp0d.ndd_th     = interp1(diiidonetwo.t(k12),diiidonetwo.ndd_th(k12),temps,'linear');
z0dinput.exp0d.ndd_nbi_th = interp1(diiidonetwo.t(k12),diiidonetwo.ndd_nbi_th(k12),temps,'linear');

% volume-averaged rotation 
%z0dinput.exp0d.wrad =  ;

z0dinput.exp0d.prad   = interp1(diiiddiag.tprad,diiiddiag.prad,temps,'linear');
%z0dinput.exp0d.ploss  =  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiplicateur pour la puissance effective dans z0dinput.cons.pnbi (from fast ion losses in onetwo simulations to metis)
% DM 31/03/2014
% The reference nbi power is reduced but the original injected power is stored in z0dinput.exp0d.pnbi
if isfield(diiidtemp,'pnb_factor')
    z0dinput.cons.pnbi = z0dinput.cons.pnbi.*interp1(diiidtemp.t,diiidtemp.pnb_factor,temps,'linear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reactivate fast ions ionisation of neutral beam
z0dinput.option.fast_ion_sbp   = 1;


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

