%============================================================
% acces au donnees de EAST via MDS+
%============================================================
function z0dinput = zerod_init_east(mode_exp,shot,gaz,temps,z0dinput)

listenbi ={'A_l','A_r','F_l','F_r'};

langue                 =  'anglais';
% cas donnees EAST
z0dinput.exp          = [];
% 1 - numero du choc
if isempty(shot) || isempty(gaz)
    prompt={'shot number :','charge of main gas:','number of time slices multiplier','NBI1 list','NBI2 list'};
    % D3D: def={'147634','1','1','15l,15r','30l,30r,33l,33r','21l,21r'};
    % def={'54828','1','1','co','cnt'};
    def={'62946','1','1','A_l,A_r','F_l,F_r'};
    dlgTitle='Access to EAST data';
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
end

%
% lecture des donnees
%
home       = getenv('HOME');
racine     = fullfile(home,'zineb/data/east');
directory  = fullfile(racine,int2str(shot));
filetemp   = [directory,'/easttemp.mat'];
fileprof   = [directory,'/eastprof.mat'];
filefit   = [directory,'/eastfit.mat'];
filetransp = [directory,'/easttransp.mat'];
fileeq   = [directory,'/easteq.mat'];
%
load(filetemp)

if ~isfield(easttemp,'rext') || ~isfield(easttemp,'zext')
    if exist(fileeq,'file');
        load(fileeq)
        if isfield(easteq,'xbord')
            disp([' the LMFS comes from easteq, with x = ',num2str(easteq.xbord)])
        else
            disp([' the LMSF comes from easteq (rext, zext)'])
        end
    else
        easteq = [];
    end
end

load(fileprof)

if exist(filefit,'file');
  load(filefit)
else
  eastfit = [];
end

if exist(filetransp,'file');
  load(filetransp)
else
  easttransp = [];
end

if ~isempty(easttemp.ip)
	tip   = easttemp.tip;
	ip    = easttemp.ip;
	ip    = ip(tip >=0);
	temps = tip(tip >=0);
	if nbtfact > 1
		temps = linspace(min(temps),max(temps),length(temps) .* nbtfact)';
		ip    = interp1(easttemp.tip,easttemp.ip,temps,'linear');
	end
else
	disp('No plasma current available')
	return
end
Ntemps=length(temps);
dtemps=mean(diff(temps));

% Plasma current : positive if counter-clockwise viewed from above the torus
signe_ip               = sign(mean(ip));
z0dinput.cons.temps    = temps;
z0dinput.cons.ip       = abs(ip);

z0dinput.cons.flux     = signe_ip .*  interp1(easttemp.t,-easttemp.psia,temps,'linear'); 

z0dinput.geo.a         = interp1(easttemp.t,easttemp.a,temps);
z0dinput.geo.z0        = 0*temps;
z0dinput.geo.R         = interp1(easttemp.t,easttemp.R0,temps);
z0dinput.geo.K         = interp1(easttemp.t,easttemp.e,temps);
z0dinput.geo.d         = 0*temps;
z0dinput.geo.vp        = 0*temps;
z0dinput.geo.sp        = 0*temps;
z0dinput.geo.sext      = 0*temps;
z0dinput.cons.iso      = 0*temps;

% Toroidal field : positive if counter-clockwise viewed from above the torus
signe_B0               = sign(mean(easttemp.b0));
z0dinput.geo.b0        = interp1(easttemp.t,abs(easttemp.b0),temps);

if isfield(easttemp,'rext') && isfield(easttemp,'zext')
  indok = find(all(isfinite(easttemp.rext) & isfinite(easttemp.zext),2));
  if length(indok) > 3
    disp('Reading separatrix from easttemp');
  	z0dinput.exp0d.Rsepa  = interp1(easttemp.trext(indok),easttemp.rext(indok,:),temps,'nearest','extrap');
 	z0dinput.exp0d.Zsepa  = interp1(easttemp.trext(indok),easttemp.zext(indok,:),temps,'nearest','extrap');
  else
	fprintf('too many NaN in Rext & Z ext, separatrix data not usable\n'); 
  end
elseif ~isempty(easteq)
  indok = find(all(isfinite(easteq.rext) & isfinite(easteq.zext),2));
  if length(indok) > 3
    disp('Reading separatrix from easttemp');
  	z0dinput.exp0d.Rsepa  = interp1(easteq.temps(indok),easteq.rext(indok,:),temps,'nearest','extrap');
 	z0dinput.exp0d.Zsepa  = interp1(easteq.temps(indok),easteq.zext(indok,:),temps,'nearest','extrap');
  else
	fprintf('too many NaN in rext & zext, separatrix data not usable\n'); 
  end
else
    disp('Separatrix data not found: Provide easttemp.rext and easttemp.zext (or easteq.rext and easteq.zext) vs time.');
    return
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

% z0dinput.cons.nbar     = max(1e13,trapz(eastprof.x',interp1(eastprof.t,eastprof.ne,temps),2) ...
%     ./ trapz(eastprof.x',ones(size(eastprof.x')),2));
% DM 15/04/2014 : subtract beam fueling from nbar (use easttemp.dnbar_nb ~= 0 and easttemp.reycling ~0 when neasser = 1)
% Parameters: Density
% DM 17/06/16 : filtrer les modulations rapides avant d'interpoler
knbar=find(easttemp.tnbar >= temps(1)-dtemps & easttemp.tnbar <= temps(end)+dtemps);
easttemp.tnbar=easttemp.tnbar(knbar);
easttemp.nbar=easttemp.nbar(knbar);
% kc  = 2*round(mean(diff(temps))/mean(diff(easttemp.tnbar)))+1;
kc  = max(3,fix((length(easttemp.tnbar) ./ Ntemps + 1) / 2) .* 2 + 1);
z0dinput.cons.nbar    = interp1(easttemp.tnbar,sgolayfilt(easttemp.nbar,1,kc),temps,'linear');
z0dinput.cons.nbar(~isfinite(z0dinput.cons.nbar) | z0dinput.cons.nbar < 0) = 0;
% figure(21);clf
% plot(easttemp.tnbar,easttemp.nbar,'g',easttemp.tnbar,sgolayfilt(easttemp.nbar,1,kc),'b',temps,z0dinput.cons.nbar,'r')

if ~isfield(easttemp,'dnbar_nb')
    z0dinput.option.neasser   = 0;
    z0dinput.option.Recycling = 0;
else
    nbar0 = z0dinput.cons.nbar(1);
    % Source term in the ODE (see zero1t.m)
    z0dinput.cons.nbar = z0dinput.cons.nbar + interp1(easttemp.tnbar,easttemp.dnbar_nb,temps,'linear');
    % Initial value in the ODE (see zero1t.m)
    z0dinput.cons.nbar(1) = nbar0;
    z0dinput.option.neasser   = 1;
    z0dinput.option.Recycling = easttemp.recycling;
end
% DM 25/09/14 : New optional parameter governing n0a
z0dinput.option.fn0a_div   = 1;

% ICRH Power (peut etre utilisee plus bas pour simuler NBI sans generation de courant/rotation (balanced injection)
z0dinput.cons.picrh   = 0 .* temps;
% From Miaohui test case
z0dinput.option.fwcd = 0;
z0dinput.option.mino = 'H';
z0dinput.option.cmin = 0.07;
z0dinput.option.nphi = 25;
z0dinput.option.freq = 28;
z0dinput.option.icrh_width = 1;
z0dinput.option.fact_mino = 0;

% LH Power @ 4.6 GHz (W)
% z0dinput.option.lhmode = 2;
%
z0dinput.option.dlh    = 0.5;
% From Miaohui test case
z0dinput.option.lhmode = 0;
z0dinput.option.etalh  = 0.75; % *1e19 ?;
z0dinput.option.wlh = 0;
z0dinput.option.npar_neg = -6.11;
z0dinput.option.upshiftmode = 'linear';
z0dinput.option.npar0 = 2.04;
z0dinput.option.xlh    = 0.45;
z0dinput.option.angle_ecrh2 = 90;
z0dinput.option.fupshift = 1.8;
z0dinput.option.freqlh = 4.6;
% z0dinput.option.dlh    = 0.07;
% DM 17/06/16 : filtrer les modulations rapides avant d'interpoler
kplh=find(easttemp.tplh >= temps(1)-dtemps & easttemp.tplh <= temps(end)+dtemps);
easttemp.tplh=easttemp.tplh(kplh);
easttemp.plh=easttemp.plh(kplh,:);
% kc  = 2*round(mean(diff(temps))/mean(diff(easttemp.tplh)))+1;
kc  = max(3,fix((length(kplh) ./ Ntemps + 1) / 2) .* 2 + 1);
z0dinput.cons.plh = sum(interp1(easttemp.tplh,sgolayfilt(easttemp.plh,1,kc),temps,'linear'),2);
z0dinput.cons.plh(~isfinite(z0dinput.cons.plh) | z0dinput.cons.plh < 0) = 0;
%figure(21);clf
% plot(easttemp.tplh,sum(easttemp.plh,2),'g',easttemp.tplh,sgolayfilt(easttemp.plh,1,kc),'b',temps,z0dinput.cons.plh,'r')

% ECRH power used for LH Power @ 2.45 GHz (W)
% z0dinput.cons.pecrh    = 0. * temps;
% DM 04/03/2013 : filtrer les modulations rapides avant d'interpoler
kpec=find(easttemp.tpec >= temps(1)-dtemps & easttemp.tpec <= temps(end)+dtemps);
easttemp.tpec=easttemp.tpec(kpec);
easttemp.pec=easttemp.pec(kpec,:);
% kc  = 2*round(mean(diff(temps))/mean(diff(easttemp.tpec)))+1;
kc  = max(3,fix((length(easttemp.tpec) ./ Ntemps + 1) / 2) .* 2 + 1);
z0dinput.cons.pecrh    = sum(interp1(easttemp.tpec,sgolayfilt(easttemp.pec,1,kc),temps,'linear'),2);
z0dinput.cons.pecrh(~isfinite(z0dinput.cons.pecrh) | z0dinput.cons.pecrh < 0) = 0;
% figure(21);clf
% plot(easttemp.tpec,sum(easttemp.pec,2),'g',easttemp.tpec,sgolayfilt(easttemp.pec,1,kc),'b',temps,z0dinput.cons.pecrh,'r')
% Modifs DM 2014
z0dinput.cons.xece          = 0.45 .* ones(size(temps));
z0dinput.option.synergie    = 1; % means no synergie
z0dinput.option.angle_ecrh  = 90;
z0dinput.option.sens        = 1; % co-current
% From Miaohui test case:
% z0dinput.cons.xece          = 0 .* temps;
% z0dinput.option.angle_ecrh  = 90;
% z0dinput.option.synergie    = 0; % means auto
% z0dinput.option.sens        = 0; % no current drive
z0dinput.option.eccdmul     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NBI Power
if sum(easttemp.pinj(:)) > 0 
     % DM 04/03/2013 : filtrer les modulations rapides avant d'interpoler
    kpinj=find(easttemp.tpinj >= temps(1)-dtemps & easttemp.tpinj <= temps(end)+dtemps);
    easttemp.tpinj=easttemp.tpinj(kpinj);
    easttemp.pinj=easttemp.pinj(kpinj,:);
     % kc = 2*round(mean(diff(temps))/mean(diff(easttemp.tpinj)))+1;
     kc = max(3,fix((length(easttemp.tpinj) ./ Ntemps + 1) / 2) .* 2 + 1);
     pinj_vect = interp1(easttemp.tpinj,sgolayfilt(easttemp.pinj,1,kc),temps,'linear');
else
     pinj_vect = zeros(Ntemps,4);
end
pinj_vect(~isfinite(pinj_vect) | pinj_vect < 0) = 0;
% figure(21);clf
% plot(easttemp.tpinj,easttemp.pinj,'g',easttemp.tpinj,sgolayfilt(easttemp.pinj,1,kc),'b',temps,pinj_vect,'r')

% MS 02.09.2008 new variable "directivity" added with counter-current sources.
% directivity = [1 1 -1 -1];

% TANGENCY RADIUS PER PINI (DIFFERENT FOR RIGHT AND LEFT SOURCES) (m)
% ORDER: co-L co-R cnt-L cnt-R
rtang_source(1) = 1.261; % co-L (LEFT SOURCE)  (m)
rtang_source(2) = 0.606; % co-R (RIGHT SOURCE) (m)
rtang_source(3) = 0.606; % cnt-L (LEFT SOURCE)  (m)
rtang_source(4) = 1.261; % cnt-R (RIGHT SOURCE) (m)
angleliste = [90,90,-90,-90]; % degrees

% listenbi ={'co','ctr'};

[void,inbi]= intersect(listenbi,list_nbi1);
fprintf('injector 1 :');fprintf('%s ',listenbi{inbi});fprintf('\n');
if ~isempty(inbi)
	pnbi1   = sum(pinj_vect(:,inbi),2);
	einj1   = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * easttemp.Einj(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi1,1));
	rtang1  = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * rtang_source(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi1,1));
	angle_nbi1 = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * angleliste(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi1,1));
	% zext1 ~=0 if Pnbi1 is off-axis
	zext1   = 0; % normalized radius
 	if any(pnbi1 > 0)
		fnbi1   = 1;
	else
		fnbi1   = 0;
	end
else
	pnbi1  = 0 * temps;
	rtang1 = 1.261;
	angle_nbi1 = 90;
	zext1  = 0;
	fnbi1  = 0;
end
[void,inbi]= intersect(listenbi,list_nbi2);
fprintf('injector 2 :');fprintf('%s ',listenbi{inbi});fprintf('\n');
if ~isempty(inbi)
	pnbi2   = sum(pinj_vect(:,inbi),2);
	einj2   = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * easttemp.Einj(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi2,1));
	rtang2  = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * rtang_source(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi2,1));
	angle_nbi2 = trapz(temps,sum(pinj_vect(:,inbi) .* (ones(size(temps)) * angleliste(inbi)),2),1) ...
        ./ max(1,trapz(temps,pnbi2,1));
	% zext2 ~=0 if Pnbi2 is off-axis
    zext2   = 0; % normalized radius
 	if any(pnbi2 > 0)
		fnbi2   = 1;
	else
		fnbi2   = 0;
	end
else
	pnbi2  = 0 * temps;
	rtang2 = 1.261;
	angle_nbi2 = -90;
	zext2  = 0;
	fnbi2  = 0;
end

z0dinput.option.einj2      = 8e4;
z0dinput.option.rtang2     = 1.261;		
z0dinput.option.angle_nbi2 = -90;		
z0dinput.option.zext2      = 0;		

if (fnbi1 == 1) && (fnbi2 == 1)
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

else
    z0dinput.option.nb_nbi     = 1;
    z0dinput.cons.pnbi         = 0 * temps;
end
%
z0dinput.option.einj      = max(1,z0dinput.option.einj);
z0dinput.option.einj2     = max(1,z0dinput.option.einj2);

z0dinput.cons.ftnbi       = 0 .* z0dinput.cons.iso;
% Modifs DM 2014
z0dinput.option.cur_nbi_time = 1;
% z0dinput.option.e_shielding  = 'Honda-Sauter';
%
% From Miaohui test case
z0dinput.option.e_shielding  = 'Lin-Liu';
z0dinput.option.nbicdmul  = 6.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Composition (zeff profile given by external data)
z0dinput.option.zeff = 0; % Given by reference and flat profile

% Default zeff reference
if ~isempty(easttemp.zeffm)
    % DM 17/06/16 : filtrer les modulations rapides avant d'interpoler
    kzeffm=find(easttemp.tzeffm >= temps(1)-dtemps & easttemp.tzeffm <= temps(end)+dtemps);
    easttemp.tzeffm=easttemp.tzeffm(kzeffm);
    easttemp.zeffm=easttemp.zeffm(kzeffm);
    % kc  = 2*round(mean(diff(temps))/mean(diff(easttemp.tzeffm)))+1;
    kc  = max(3,fix((length(easttemp.tzeffm) ./ Ntemps + 1) / 2) .* 2 + 1);
    z0dinput.cons.zeff    = interp1(easttemp.tzeffm,sgolayfilt(easttemp.zeffm,1,kc),temps,'linear');
    z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff) | z0dinput.cons.zeff < 0) = 0;
    % figure(21);clf
    % plot(easttemp.tzeffm,easttemp.zeffm,'g',easttemp.tzeffm,sgolayfilt(easttemp.zeffm,1,kc),'b',temps,z0dinput.cons.zeff,'r')
else
    z0dinput.cons.zeff = 2 .* ones(size(temps));
end
% Compute zeff reference from zeff profiles in transp simulations (only when the transp profile is not flat)
kzf = find(max(easttransp.zeff,[],2)-min(easttransp.zeff,[],2) > 0);
if ~isempty(kzf)
    % space interpolation on the easttemp rho grid
    zeffp = interp1(easttransp.x,easttransp.zeff(kzf,:)',easttemp.x)';
    % time interpolation on the easttransp time grid
    xlao=interp1(easttemp.t,easttemp.xlao,easttransp.t(kzf),'linear');
    % space interpolation on a fixed metis grid and line-average using metis coordinates (Rmax - Rmin)/ 2a
    xmetis=linspace(0,1,21);
    zeffl=zeros(length(kzf),21);
    for kk=1:length(kzf)
        zeffl(kk,:)=interp1(xlao(kk,:),zeffp(kk,:)',xmetis)';
    end
    zeffl=trapz(xmetis,zeffl,2);
    % time interpolation
    z0dinput.cons.zeff(temps >= easttransp.t(kzf(1)) & temps <= easttransp.t(kzf(end))) = ...
        interp1([-1;easttransp.t(kzf);max(easttransp.t(kzf(end)),temps(end))+1],[zeffl(1);zeffl;zeffl(end)], ...
        temps(temps >= easttransp.t(kzf(1)) & temps <= easttransp.t(kzf(end))),'linear');
end
% figure(5);clf;
% plot(z0dinput.cons.temps,z0dinput.cons.zeff,'bo')
% if isfield(easttemp,'zf_factor'), z0dinput.cons.zeff = easttemp.zf_factor.*z0dinput.cons.zeff; end;
% hold on;plot(z0dinput.cons.temps,z0dinput.cons.zeff,'ro')
% z0dinput.cons.zeff(z0dinput.cons.zeff < 1) = 1;
% plot(z0dinput.cons.temps,z0dinput.cons.zeff,'r')

% Tension par tour et li initiaux (Breakdown & burnthrough) avec runaways (Current diffusion and Equilibrium).
% option.breakdown = electric field at starting time of the simulation in unit of Dreicer electric field if > 0
% or in Volt per turn if < 0 (used abs value, convert internally into electric field)
% vsurf=interp1(easttemp.tvsurf,easttemp.vsurf,temps,'linear');
% z0dinput.option.breakdown = -vsurf(1);
% DM 17/06/16 : filtrer les modulations rapides avant d'interpoler
kvmes=find(easttemp.tvmes >= temps(1)-dtemps & easttemp.tvmes <= temps(end)+dtemps);
easttemp.tvmes=easttemp.tvmes(kvmes);
easttemp.vmes=easttemp.vmes(kvmes);
% kc  = 2*round(mean(diff(temps))/mean(diff(easttemp.tvmes)))+1;
kc  = max(3,fix((length(easttemp.tvmes) ./ Ntemps + 1) / 2) .* 2 + 1);
vmes    = interp1(easttemp.tvmes,sgolayfilt(easttemp.vmes,1,kc),temps,'linear');
vmes(~isfinite(vmes)) = 0;
% figure(21);clf
% plot(easttemp.tvmes,easttemp.vmes,'g',easttemp.tvmes,sgolayfilt(easttemp.vmes,1,kc),'b',temps,vmes,'r')
z0dinput.option.breakdown = -vmes(1);

li = interp1(easttemp.tli,easttemp.li,temps,'linear');
z0dinput.option.li  = li(1);

z0dinput.option.runaway = 1;

% Neoclassical: facteur multiplicatif pour le courant de bootstrap
z0dinput.option.modeboot = 1; % Sauter
if ~isfield(easttemp,'bs_factor')
    z0dinput.option.bootmul = 1;
else
    z0dinput.option.bootmul = easttemp.bs_factor;
end

% Current diffusion and equilibrium: Current control or flux control and new regularisation du profil de q au centre
z0dinput.option.vloop = 0;
% z0dinput.option.vloop = 4;
z0dinput.option.tswitch = Inf;
z0dinput.option.cronos_regul = 0;

% Confinement and Transport
if ~isfield(easttemp,'H_factor')
    z0dinput.cons.hmore = ones(size(temps));
else
    z0dinput.cons.hmore = interp1(easttemp.t,easttemp.H_factor,temps,'linear');
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
%   End of Modifs DM 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% activation des dents de scies	
z0dinput.option.qdds         = 1;
z0dinput.option.kidds        = 3;


z0dinput.option.zmax      = 18;
z0dinput.option.zimp      = 6;
z0dinput.option.rimp      = 0.1;
z0dinput.option.modeh     = 1;
switch gaz 
    case 1
        z0dinput.option.gaz =  2;
    case {5,11}
        error('Not yet implemented');
    otherwise
        z0dinput.option.gaz =  4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0dinput.option.ane = 4;
% estimation grossiere du piquage moyen avec les profils de densité.
% neav=sum(eastprof.ne.*(ones(length(eastprof.t),1)*[0 diff(eastprof.x.*eastprof.x)']),2);
% z0dinput.option.vane      = mean(eastprof.ne(:,1)./neav);
z0dinput.option.vane = mean((sum(easttransp.dv,2).*easttransp.ne(:,1))./sum(easttransp.dv.*easttransp.ne,2));

% From Miaohui test case
% z0dinput.option.ane = 4;
% z0dinput.option.vane = 1.7;

% z0dinput.option.configuration = 2;
z0dinput.machine     = 'EAST';
z0dinput.shot        = shot;

% signe B-T * Ip	
z0dinput.option.signe = signe_B0*signe_ip;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% donnees experimentale pour comparaison

z0dinput.exp0d.temps  = temps;
z0dinput.exp0d.pnbi   = [real(z0dinput.cons.pnbi) imag(z0dinput.cons.pnbi)];
z0dinput.exp0d.plh    = z0dinput.cons.plh;
z0dinput.exp0d.pecrh  = z0dinput.cons.pecrh;
z0dinput.exp0d.picrh  = z0dinput.cons.picrh;
z0dinput.exp0d.flux   = z0dinput.cons.flux;

z0dinput.exp0d.ip            = ip;
z0dinput.exp0d.betaptot  = interp1(easttemp.t,easttemp.betap,temps,'linear');
z0dinput.exp0d.betan     = interp1(easttemp.t,easttemp.betan,temps,'linear');
% z0dinput.exp0d.vloop  = vsurf;
z0dinput.exp0d.vmes   = interp1(easttemp.tvmes,easttemp.vmes,temps,'linear');

z0dinput.exp0d.q95    = interp1(easttemp.tq95,easttemp.q95,temps,'linear');
z0dinput.exp0d.q0     = interp1(easttemp.tq0,easttemp.q0,temps,'linear');

z0dinput.exp0d.nbar  = interp1(easttemp.tnbar,easttemp.nbar,temps,'linear');
z0dinput.exp0d.li    = li;

z0dinput.exp0d.prad   = interp1(easttemp.tprad,easttemp.prad,temps,'linear');

ktr=zeros(1,length(easttransp.t));
for i=1:length(ktr)
    [~,ktr(i)]=min(abs(temps-easttransp.t(i)));
end
if isempty(easttemp.zeffm)
    z0dinput.exp0d.zeff      = zeros(Ntemps,1);
    z0dinput.exp0d.zeff(ktr) = sum(easttransp.dv.*easttransp.zeff,2)./sum(easttransp.dv,2);
else
    z0dinput.exp0d.zeff      = interp1(easttemp.t,easttemp.zeffm,temps,'linear');
end
z0dinput.exp0d.ne0      = NaN(Ntemps,1);
z0dinput.exp0d.ne0(ktr) = easttransp.ne(:,1);
z0dinput.exp0d.te0      = NaN(Ntemps,1);
z0dinput.exp0d.te0(ktr) = easttransp.te(:,1);
z0dinput.exp0d.ti0      = NaN(Ntemps,1);
z0dinput.exp0d.ti0(ktr) = easttransp.ti(:,1);
z0dinput.exp0d.w        = NaN(Ntemps,1);
z0dinput.exp0d.w(ktr)   = easttransp.w;
z0dinput.exp0d.wth      = NaN(Ntemps,1);
z0dinput.exp0d.wth(ktr) = easttransp.wth;
z0dinput.exp0d.esup_nbi = NaN(Ntemps,1);
z0dinput.exp0d.esup_nbi(ktr) = easttransp.wsupnb;
z0dinput.exp0d.taue     = NaN(Ntemps,1);
z0dinput.exp0d.taue(ktr) = easttransp.taueth; 

z0dinput.exp0d.iboot       = NaN(Ntemps,1);
z0dinput.exp0d.iboot(ktr)  = easttransp.iboot;
z0dinput.exp0d.inbicd      = NaN(Ntemps,1);
z0dinput.exp0d.inbicd(ktr) = easttransp.inbicd;
z0dinput.exp0d.ilh         = NaN(Ntemps,1);
z0dinput.exp0d.ilh(ktr)    = easttransp.ilh;
z0dinput.exp0d.ieccd       = NaN(Ntemps,1);
% z0dinput.exp0d.ieccd(ktr)  = easttransp.ieccd;
z0dinput.exp0d.icd         = NaN(Ntemps,1);
z0dinput.exp0d.icd(ktr)    = easttransp.icd;
z0dinput.exp0d.ini         = NaN(Ntemps,1);
z0dinput.exp0d.ini(ktr)    = easttransp.ini;
z0dinput.exp0d.iohm        = z0dinput.exp0d.ip - z0dinput.exp0d.ini;

z0dinput.exp0d.pin         = NaN(Ntemps,1);
z0dinput.exp0d.pin(ktr)    = easttransp.pin;
z0dinput.exp0d.pohm        = NaN(Ntemps,1);
z0dinput.exp0d.pohm(ktr)   = easttransp.pohm;
% volume-averaged rotation 
z0dinput.exp0d.wrad        = NaN(Ntemps,1);
z0dinput.exp0d.wrad(ktr)   = sum(easttransp.dv.*easttransp.vtor,2)./sum(easttransp.dv,2)./z0dinput.geo.R(ktr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiplicateur pour la puissance effective dans z0dinput.cons.pnbi (from fast ion losses in transp simulations to metis)
% DM 31/03/2014
% The reference nbi power is reduced but the original injected power is stored in z0dinput.exp0d.pnbi
if isfield(easttemp,'pnb_factor')
    z0dinput.cons.pnbi = z0dinput.cons.pnbi.*interp1(easttemp.t,easttemp.pnb_factor,temps,'linear');
end

% reactivate fast ions ionisation of neutral beam
z0dinput.option.fast_ion_sbp   = 1;



return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

return