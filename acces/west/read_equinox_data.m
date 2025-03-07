function [Rsepa,Zsepa,ip,li,betap,qax,q95,wmhd,R0,Z0,Rax,Zax,aminor,K,d,dlow,dup,vloop,flux_edge]=read_equinox_data(shot,time,varargin)

% try with NICE + polarimeter
time_protected_ = time;
varargin_protected = varargin;

try
  equi_  = imas_west_get(shot,'equilibrium',0,1);
  noms_ = fieldnames(equi_);
  for k_=1:length(noms_)
    eval(sprintf('%s = equi_.%s;',noms_{k_},noms_{k_}));
  end
catch
  equi_ = [];
end
time = time_protected_;
varargin = varargin_protected;
% old fassion 
if isempty(equi_)
    tpn = tempname;
    pwd_mem = pwd;
    mkdir(tpn)
    cd(tpn)

    cmd = fullfile(fileparts(which('zerod_init_west')),'py','equilibrium_python_to_matlab.py');
    s = unix(sprintf('python %s %d',cmd,shot));
    if s ~= 0
      cd(pwd_mem);
      unix(sprintf('rmdir %s',tpn));  
      error('No EQUINOX data available')
    end
    filename = sprintf('data_vactheqx_Shot%d_Run0_Occ0_imas_public_west.mat',shot);
    if exist(filename,'file')
	load(filename);
    else
	filename = fullfile(fileparts(which('zerod_init_west')),'tmp','matlab_equilibria',filename);
	load(filename);
    end
    delete(filename);
    cd(pwd_mem);
    unix(sprintf('rmdir %s',tpn));  
end

temps = timeEquiIDS';
% rescale on ignitron
temps = temps - tsbase(shot,'RIGNITRON');


ip    = interp1(temps,ip',time,'pchip','extrap');
li    = interp1(temps,li_3',time,'pchip','extrap');
betap = interp1(temps,beta_pol',time,'pchip','extrap');
qax   = interp1(temps,q_axis',time,'pchip','extrap');
q95   = interp1(temps,q_95',time,'pchip','extrap');
wmhd  = interp1(temps,w_mhd',time,'pchip','extrap');
Rax   = interp1(temps,mag_ax_R',time,'pchip','extrap');
Zax   = interp1(temps,mag_ax_Z',time,'pchip','extrap');
K     = interp1(temps,elong(:,end),time,'pchip','extrap');
dlow  = interp1(temps,triang_low(:,end),time,'pchip','extrap');
dup   = interp1(temps,triang_up(:,end),time,'pchip','extrap');
d     = (dlow + dup) / 2;


nbp = 201;
Rsepa    = NaN * ones(size(time,1),nbp);
Zsepa    = NaN * ones(size(time,1),nbp);
aminor   = NaN * ones(size(time));
R0       = NaN * ones(size(time));
Z0       = NaN * ones(size(time));
K        = NaN * ones(size(time));
d        = NaN * ones(size(time));

Rsepa_eq    = NaN * ones(size(temps,1),nbp);
Zsepa_eq   = NaN * ones(size(temps,1),nbp);
aminor_eq   = NaN * ones(size(temps));
R0_eq       = NaN * ones(size(temps));
Z0_eq       = NaN * ones(size(temps));
K_eq       = NaN * ones(size(temps));
d_eq        = NaN * ones(size(temps));

for k=1:length(temps)
  Rext = squeeze(boundPlasma(k,1,:));
  Zext = squeeze(boundPlasma(k,2,:));
  [Rsepa_eq(k,:),Zsepa_eq(k,:),aminor_eq(k),R0_eq(k),Z0_eq(k),K_eq(k),d_eq(k)] = reshape_LCFS(Rext,Zext,nbp);
end

for l=1:nbp
    Rsepa(:,l)      = interp1(temps,Rsepa_eq(:,l),time,'pchip','extrap');
    Zsepa(:,l)      = interp1(temps,Zsepa_eq(:,l),time,'pchip','extrap');
end

aminor  = interp1(temps,aminor_eq,time,'pchip','extrap');
R0      = interp1(temps,R0_eq,time,'pchip','extrap');
Z0      = interp1(temps,Z0_eq,time,'pchip','extrap');
K       = interp1(temps,K_eq,time,'pchip','extrap');
d       = interp1(temps,d_eq,time,'pchip','extrap');

% compute Vloop
vloop = z0dxdt(psi_prof(:,end),temps);
vloop = interp1(temps,vloop,time,'pchip','extrap');
flux_edge = interp1(temps,psi_prof(:,end),time,'pchip','extrap');

function [Rsepa,Zsepa,a,ra,za,k,d] = reshape_LCFS(Rext,Zext,nbp)


Rsepa = [];
Zsepa = [];
a = [];
ra = [];
za = [];
k = [];
d = [];

indbad = 0;
indok  = find(isfinite(Rext) & isfinite(Zext));
if length(indok) < 5
    indbad = 1;
else
    r = Rext(indok);
    r = r(:);
    z = Zext(indok);
    z = z(:);		    
    KH = sort(unique(convhull(r,z)));
    if (length(KH) ~= length(r))
	    r = r(KH);
	    z = z(KH);
	    r = r(:);
	    z = z(:);
    end
    if length(r) < 5
	indbad = 1;
    else
	r0   = (min(r) + max(r)) ./ 2;
	z0   = (min(z) + max(z)) ./ 2;
	cc   = (r - r0) + sqrt(-1) .* (z - z0);
	thc  = unwrap(angle(cc));
	thc(thc <0) = thc(thc<0) + 2 .* pi;
	rhoc = abs(cc);
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(1,rhoc,rhoc,rhoc);
	thc = cat(1,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	while (length(thc) > 3) && ~isempty(indnok)
	  thc(indnok) =[];
	  rhoc(indnok)  = [];
	  indnok = find(diff(thc)<=0);
	end  
	teta  = linspace(0,2*pi,nbp)';
	rho = spline(thc,rhoc,teta);
	R  = r0 + rho .* cos(teta);
	Z  = z0 + rho .* sin(teta);
	R(end) = (R(1)+ R(end)) ./ 2;
	R(1)   = R(end);
	Z(end) = (Z(1)+ Z(end)) ./ 2;
	Z(1)   = Z(end);
	
	if any(~isfinite(R)) || any(~isfinite(Z))		    
	  indbad = 1;
	end
	figure(21);clf;plot(Rext,Zext,'b',R,Z,'.r');drawnow
	%pause(0.1);
    end
end

if indbad == 1
  return
end

Rext = R';
Zext = Z';		
if any(~isfinite(Rext(:))) || any(~isfinite(Zext(:)))
  keyboard
end 


% recalcul des parametres pour verification
ve    = ones(1,size(Rext,2));
rmin  = min(Rext,[],2);
rmax  = max(Rext,[],2);
ra    = max(0.1,0.5 .* (rmin + rmax));
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
Rsepa = Rext;
%Zsepa = Zext - za * ones(1,size(Zext,2));
Zsepa = Zext;


