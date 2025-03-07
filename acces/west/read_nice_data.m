function [Rsepa,Zsepa,ip,li,betap,qax,q95,wmhd,R0,Z0,Rax,Zax,aminor,K,d,dlow,dup,vloop,flux_edge]=read_nice_data(shot,time,varargin)

% init data
Rsepa			= [];
Zsepa			= [];
ip			= [];
li			= [];
betap			= [];
qax			= [];
q95			= [];
wmhd			= [];
R0			= [];
Z0			= [];
Rax			= [];
Zax			= [];
aminor			= [];
K			= [];
d			= [];
dlow			= [];
dup			= [];
vloop			= [];
flux_edge		= [];

% read NICE data
equi = imas_west_get(shot,'equilibrium',0,1);
if isempty(equi)
  disp('NICE + Polarimetry not available: try NICE with magnetic alone');
  equi = imas_west_get(shot,'equilibrium');
end
if isempty(equi)
  fprintf('NICE data are not available for shot %d\n',shot);
  return
end 

temps = equi.time';
% rescale on ignitron
temps = temps - tsbase(shot,'RIGNITRON');

% good time only
good_time = find(equi.code.output_flag >= 0);
if ~isempty(good_time)
  temps = temps(good_time);
else
  disp('No valid time for this shot')
  return
end

ip    =  abs(interp1(temps,equi.ip(good_time)',time,'linear','extrap'));
li    =  abs(interp1(temps,equi.li_3(good_time)',time,'linear','extrap'));
betap =  abs(interp1(temps,equi.beta_pol(good_time)',time,'linear','extrap'));
qax   =  abs(interp1(temps,equi.q_axis(good_time)',time,'linear','extrap'));
q95   =  abs(interp1(temps,equi.q_95(good_time)',time,'linear','extrap'));
wmhd  =  abs(interp1(temps,equi.w_mhd(good_time)',time,'linear','extrap'));
Rax   =  abs(interp1(temps,equi.mag_ax_R(good_time)',time,'linear','extrap'));
Zax   =  abs(interp1(temps,equi.mag_ax_Z(good_time)',time,'linear','extrap'));
K     =  abs(interp1(temps,equi.elong(good_time,end),time,'linear','extrap'));
dlow  =  abs(interp1(temps,equi.triang_low(good_time,end),time,'linear','extrap'));
dup   =  abs(interp1(temps,equi.triang_up(good_time,end),time,'linear','extrap'));
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
  Rext = squeeze(equi.boundPlasma(good_time(k),1,:));
  Zext = squeeze(equi.boundPlasma(good_time(k),2,:));
  [Rsepa_eq(k,:),Zsepa_eq(k,:),aminor_eq(k),R0_eq(k),Z0_eq(k),K_eq(k),d_eq(k)] = reshape_LCFS(Rext,Zext,nbp);
end

for l=1:nbp
    Rsepa(:,l)      = interp1(temps,Rsepa_eq(:,l),time,'linear',NaN);
    Zsepa(:,l)      = interp1(temps,Zsepa_eq(:,l),time,'linear',NaN);
end

aminor  = interp1(temps,aminor_eq,time,'linear','extrap');
R0      = interp1(temps,R0_eq,time,'linear','extrap');
Z0      = interp1(temps,Z0_eq,time,'linear','extrap');
K       = interp1(temps,K_eq,time,'linear','extrap');
d       = interp1(temps,d_eq,time,'linear','extrap');

% compute Vloop
vloop = z0dxdt(equi.psi_prof(good_time,end),temps);
vloop = interp1(temps,vloop,time,'linear','extrap');
flux_edge = interp1(temps,equi.psi_prof(good_time,end),time,'linear','extrap');

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
%      KH = sort(unique(convhull(r,z)));
%      if (length(KH) ~= length(r))
%  	    r = r(KH);
%  	    z = z(KH);
%  	    r = r(:);
%  	    z = z(:);
%      end
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
	%figure(21);clf;plot(Rext,Zext,'b',R,Z,'.r');drawnow
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


