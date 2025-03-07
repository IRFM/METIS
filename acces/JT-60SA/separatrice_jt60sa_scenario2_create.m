% separatrice pour Ramp_up JT-60SA d'apres donnees CREATE
% exemple :




function [R,Z] = separatrice_jt60sa_scenario2_create(ip,temps)



info = load('LCFS_JT-60SA_Scenario2_from_CREATE.mat');
% tableau des ips
urz      = info.LCFS;
ip_ref   = NaN .* ones(1,length(urz));
time_ref = NaN .* ones(1,length(urz));
for k= 1:length(urz)
	ip_ref(k)   = urz(k).ip;
	time_ref(k) = urz(k).time;
end
uni =  linspace(0, 2.* pi,201)';
r0  = NaN * ones(1,length(ip_ref));
z0  = NaN * ones(1,length(ip_ref));
rho = NaN * ones(length(uni),length(ip_ref));
for k= 1:length(ip_ref)

	% extraction de la sepa
	r = urz(k).r_lim;
	z = urz(k).z_lim;
	r = r(:);
	z = z(:);
	KH = sort(unique(convhull(r,z)));
	if (length(KH) ~= length(r))
	    index_full = 1:length(r);
	    r = r(KH);
	    z = z(KH);
	    r = interp1(KH,r,index_full,'linear');
	    z = interp1(KH,z,index_full,'linear');
	    indbad_lcfs = find(~isfinite(r) | ~isfinite(z));
	    if ~isempty(indbad_lcfs)
		r(indbad_lcfs) = [];
		z(indbad_lcfs) = [];
	    end
	    r = r(:);
	    z = z(:);
	end

	r0(k)   = (min(r) + max(r)) ./ 2;
	mask    = (r == max(r));
	z0(k)   = sum(z .* mask) ./ max(1,sum(mask));
	cc   = (r - r0(k)) + sqrt(-1) .* (z - z0(k));
	thc  = unwrap(angle(cc));
	thc(thc <0) = thc(thc<0) + 2 .* pi;
	rhoc = abs(cc);
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(1,rhoc,rhoc,rhoc);
	thc = cat(1,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	thc(indnok) =[];
	rhoc(indnok)  = [];
	rho(:,k) = spline(thc,rhoc,uni);
end

if max(ip) > ip_ref(3)
    ip_ref(4:5) = max(ip);
else
    ip = ip ./ max(ip) .* max(ip_ref);
end

% ajoute une separatrice pour les bas courants
a_ini     = 1.043;
R_ini     = 2.8;
ip_ini    = 1e-3;
time_ini  = 0;
r0        = cat(2,R_ini,r0);
z0        = cat(2,0,z0);
ip_ref    = cat(2,ip_ini,ip_ref);
time_ref  = cat(2,time_ini,time_ref);
rho       = cat(2,a_ini*ones(size(rho,1),1),rho);

% interpolation sur ip ramp-up et flat-top
rho_out = pchip(ip_ref(1:end-2),rho(:,1:end-2),ip(:));
r0_out = pchip(ip_ref(1:end-2),r0(1:end-2),ip(:))';
z0_out = pchip(ip_ref(1:end-2),z0(1:end-2),ip(:))';

%
indlow  = find(ip <= min(ip_ref));
if ~isempty(indlow)
	rho_out(:,indlow) = rho(:,1) * ones(1,length(indlow));
	r0_out(indlow)    = r0(1);
	z0_out(indlow)    = z0(1);	
end
indhigh  = find(ip >= max(ip_ref));
if ~isempty(indhigh)
	rho_out(:,indhigh) = rho(:,end-2) * ones(1,length(indhigh));
	r0_out(indhigh)    = r0(end-2);
	z0_out(indhigh)    = z0(end-2);	
end

% interpolation between start and end of flattop
indipmax = find(abs(ip - max(ip)) < (1e-2 .*max(ip)));
tt       = cat(1,temps(indipmax(1)),temps(indipmax(end)));
rho_out(:,indipmax)  = pchip(tt,rho(:,end-2:end-1),temps(indipmax));
r0_out(indipmax)  = pchip(tt,r0(end-2:end-1),temps(indipmax))';
z0_out(indipmax)  = pchip(tt,z0(end-2:end-1),temps(indipmax))';


% ajout du point final
a_end     = 0.6;
R_end     = 2.51;
Z_end     = -1.2869;
ip_end    = 0.24;
time_end  = 140;
r0        = cat(2,r0,R_end);
z0        = cat(2,z0,Z_end);
ip_ref    = cat(2,ip_ref,ip_end);
time_ref  = cat(2,time_ref,time_end);
rho       = cat(2,rho,a_end*ones(size(rho,1),1));

% detection de dip/dt < 0
dipdt = z0dxdt(ip,temps);
ind_rampdown = NaN;
for k=1:length(dipdt)
  if all(sign(dipdt(k:end)) < 0)
      ind_rampdown = k;
      break;
  end
end
% interpolation sur ip ramp-down 
if isfinite(ind_rampdown)
  rho_out(:,ind_rampdown:end) = pchip(ip_ref(end-2:end),rho(:,end-2:end),ip(ind_rampdown:end));
  r0_out(ind_rampdown:end) = pchip(ip_ref(end-2:end),r0(end-2:end),ip(ind_rampdown:end))';
  z0_out(ind_rampdown:end) = pchip(ip_ref(end-2:end),z0(end-2:end),ip(ind_rampdown:end))';
end

R  = (ones(size(uni)) * r0_out + rho_out .* cos(uni * ones(1,size(rho_out,2))))';
Z  = (ones(size(uni)) * z0_out + rho_out .* sin(uni * ones(1,size(rho_out,2))))';

%  figure
%  plot(R',Z');
% keyboard

