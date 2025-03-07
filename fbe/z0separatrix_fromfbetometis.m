%  Z0SEPARATRIX_FROMFBETOMETIS  function used  to importe LCFS from FBE codes
%------------------------------------------------------------------------------- 
% file:  z0separatrix_fromfbetometis.m  ->  z0separatrix_fromfbetometis
% 
% 
% function Matlab2012b: 
% 
% z0separatrix_fromfbetometis is a function used  to importe LCFS from FBE codes
%  
% syntax:  
%   z0dinput = z0separatrix_fromfbetometis(z0dinput,fbe_lcfs,{mode_resampling,nbtheta,plotonoff}) 
%  
% input:  
%  z0dinput        = input structure from METIS
%  fbe_lcfs        = structure containing LCFS information from FBE codes
%  mode_resampling = method to remsample LCFS,
%  nbtheta         = number of points in LCFS outputed by the function
%  plotonoff       = flag to turn "on" (1) or "off" (0) diagnostic graphs. 
%
% fields of fbe_lcfs:
%  fbe_lcfs = array of cells
%  fbe_lcfs{k}.time = time in s
%  fbe_lcfs{k}.ip   = plasma current in A
%  fbe_lcfs{k}.R    = vector of R coordinate of points of the LCFS
%  fbe_lcfs{k}.Z    = vector of Z coordinate of points of the LCFS
%
% possible values of mode_resampling:
%   'input'= kept time slices from fbe_lcfs
%   'time' = resampling of time
%   'ip'   = resampling of current value with handling of ramp-up and ramp-down
%   'none' = no resampling (Ip and time vector must be the same in z0dinput and fbe_lcfs)
%  
% sorties :  
%   z0dinput  = modified input structure from METIS containing new LCFS and time slices vetor or Ip waveform.
%  
% function writed by J-F Artaud
% version  1.0  (20180413)
%  
%-------------------------------------------------------------------------------  
%  
function z0dinput = z0separatrix_fromfbetometis(z0dinput,fbe_lcfs,mode_resampling,nbtheta,plotonoff)

% test input
if nargin < 3
  mode_resampling = 'time';
elseif isempty(mode_resampling)
  mode_resampling = 'time';
end
if nargin < 4
  nbtheta = 201;
elseif isempty(nbtheta)
  nbtheta = 201;
end
if nargin < 5
  plotonoff = 0;
elseif isempty(plotonoff)
  plotonoff = 0;
end

% mode_resampling    = 'ip' or 'time'

% iso angle LCFS
ip_ref   = NaN .* ones(1,length(fbe_lcfs));
time_ref = NaN .* ones(1,length(fbe_lcfs));
for k= 1:length(fbe_lcfs)
	ip_ref(k)   = fbe_lcfs{k}.ip;
	time_ref(k) = fbe_lcfs{k}.time;
end
uni =  linspace(0, 2.* pi,nbtheta)';
r0  = NaN * ones(1,length(ip_ref));
z0  = NaN * ones(1,length(ip_ref));
aminor  = NaN * ones(1,length(ip_ref));
elong  = NaN * ones(1,length(ip_ref));
delta  = NaN * ones(1,length(ip_ref));
rho = NaN * ones(length(uni),length(ip_ref));
for k= 1:length(ip_ref)

	% extraction de la sepa
	r = fbe_lcfs{k}.R;
	z = fbe_lcfs{k}.Z;
	r = r(:);
	z = z(:);
    indbad = find(r <= 0);
    if ~isempty(indbad)
        r(indbad) = [];
        z(indbad) = [];
    end
% 	KH = sort(unique(convhull(r,z)));
% 	if (length(KH) ~= length(r))
% 	    index_full = 1:length(r);
% 	    r = r(KH);
% 	    z = z(KH);
% 	    r = interp1(KH,r,index_full,'linear');
% 	    z = interp1(KH,z,index_full,'linear');
% 	    indbad_lcfs = find(~isfinite(r) | ~isfinite(z));
% 	    if ~isempty(indbad_lcfs)
% 		r(indbad_lcfs) = [];
% 		z(indbad_lcfs) = [];
% 	    end
% 	    r = r(:);
% 	    z = z(:);
% 	end

	r0(k)   = (min(r) + max(r)) ./ 2;
	%mask    = (r == max(r));
	%z0(k)   = sum(z .* mask) ./ max(1,sum(mask));
	z0(k)   = (min(z) + max(z)) ./ 2;
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
    
    % calcul des moments
    xu    = linspace(0,1,length(r))';
    dRdx  = pdederive(xu,r,2,2,1,1);
    % recalcul des parametres sur le vecteur final
    rmin  = min(r);
    rmax  = max(r);
    aminor(k) = 0.5 .* (rmax - rmin);
    zmin  = min(z);
    zmax  = max(z);
    elong(k)    = (abs(trapz(xu,z .*  dRdx) ./ pi ./ aminor(k)) + (zmax - zmin)) ./ 3 ./ aminor(k);
    rzmax = r(min(find(z == zmax)));
    rzmin = r(min(find(z== zmin)));
    uu   =  angle(rzmax - r0(k) + sqrt(-1) .* zmax);
    ul   =  angle(rzmin - r0(k) + sqrt(-1) .* zmin);
    tu   =  abs((acos((rzmax - r0(k)) ./ aminor(k)) - acos(cos(uu))) ./ sin(uu));
    tl   =  abs((acos((rzmin - r0(k)) ./ aminor(k)) - acos(cos(ul))) ./ sin(ul));
    tm   =  (tl + tu) ./ 2;
    d    =   abs(rzmax + rzmin -  2 .* r0(k)) ./ 2 ./ aminor(k);
    delta(k) =  0.6 .* d + 0.4  .* sin(tm);
    
end

if ~all(diff(time_ref)> 0)
    disp('Warning : Time slices of imported LCFS set are not in an causal chronological order');
end

% figure;plot(time_ref,r0,time_ref,z0);drawnow
% figure;plot(time_ref,rho);drawnow
%keyboard
% figure
% for k= 1:length(ip_ref)
%     r = fbe_lcfs{k}.R;
%     z = fbe_lcfs{k}.Z;
%     plot(r0(k) + rho(:,k) .* cos(uni), z0(k) + rho(:,k) .* sin(uni),'b',r,z,'r.')
%     pause
% end

% anti spike before resampling
% if length(r0) >=3
%     rho = sgolayfilt(rho',1,3)';
%     r0  = sgolayfilt(r0',1,3)';
%     z0  = sgolayfilt(z0',1,3)';
%     aminor  = sgolayfilt(aminor',1,3)';
%     elong  = sgolayfilt(elong',1,3)';
%     delta  = sgolayfilt(delta',1,3)';
% end


% rescale ip if input waveform has greater IP than ip_ref
switch mode_resampling
case 'input'
  % nothing: controlled by input
otherwise
  ip_ref = ip_ref ./   max(ip_ref) .* max(z0dinput.cons.ip);
end

% add backup value if start and end of METIS simulation is not included in given LCFS
switch mode_resampling
case 'input'
  % search for time offset
  indmax = find(z0dinput.cons.ip == max(z0dinput.cons.ip),1);
  indoffset = find(z0dinput.cons.ip >= ip_ref(1),1);
  if indoffset <= indmax
    time_offset = time_ref(1) - z0dinput.cons.temps(indoffset);
    fprintf('Time offset for synchronising METIS simulation and FBE is = %g\n',time_offset);
  else
    time_offset = 0;
    disp('unable to synchronise METIS simulation on FBE input');
  end
  % kept input time
  noms  = fieldnames(z0dinput.cons);
  temps = z0dinput.cons.temps;
  for k=1:length(noms)
    z0dinput.cons.(noms{k}) = interp1(temps + time_offset,z0dinput.cons.(noms{k}),time_ref','linear','extrap');
  end
  noms  = fieldnames(z0dinput.geo);
  for k=1:length(noms)
    z0dinput.geo.(noms{k}) = interp1(temps + time_offset,z0dinput.geo.(noms{k}),time_ref','linear','extrap');
  end
  noms  = fieldnames(z0dinput.exp0d);
  for k=1:length(noms)
    z0dinput.exp0d.(noms{k}) = interp1(temps + time_offset,z0dinput.exp0d.(noms{k}),time_ref','linear','extrap');
  end
  % resampling in time
  rho_out = rho;
  r0_out  = r0;
  z0_out  = z0;
  a_out   = aminor;
  K_out   = elong;
  d_out   = delta;
  ip_out  = ip_ref';
  time_out = time_ref';
 
case 'time'
  if time_ref(1) > z0dinput.cons.temps(1)
      % adding LCFS for first time slice
      a_ini     = min(z0dinput.geo.a(1),aminor(1));
      R_ini     = z0dinput.geo.R(1);
      d_ini     = min(z0dinput.geo.d(1),delta(1));
      K_ini     = min(z0dinput.geo.K(1),elong(1));
      z0_ini    = z0dinput.geo.z0(1);
      ip_ini    = z0dinput.cons.ip(1);
      time_ini  = z0dinput.cons.temps(1);
      r0        = cat(2,R_ini,r0);
      z0        = cat(2,z0_ini,z0);
      aminor    = cat(2,a_ini,aminor);
      elong     = cat(2,K_ini,elong);
      delta     = cat(2,d_ini,delta);
      ip_ref    = cat(2,ip_ini,ip_ref);
      time_ref  = cat(2,time_ini,time_ref);
      rho       = cat(2,a_ini*ones(size(rho,1),1),rho);  
  end
  if time_ref(end) < z0dinput.cons.temps(end)
      % adding LCFS for last time slice
      a_ini     = min(z0dinput.geo.a(end),aminor(end));
      R_ini     = z0dinput.geo.R(end);
      d_ini     = min(z0dinput.geo.d(end),delta(end));
      K_ini     = min(z0dinput.geo.K(end),elong(end));
      z0_ini    = z0dinput.geo.z0(end);
      ip_ini    = z0dinput.cons.ip(end);
      time_ini  = z0dinput.cons.temps(end);
      r0        = cat(2,r0,R_ini);
      z0        = cat(2,z0,z0_ini);
      aminor    = cat(2,aminor,a_ini);
      elong     = cat(2,elong,K_ini);
      delta     = cat(2,delta,d_ini);
      ip_ref    = cat(2,ip_ref,ip_ini);
      time_ref  = cat(2,time_ref,time_ini);
      rho       = cat(2,rho,a_ini*ones(size(rho,1),1));  
  end
  % resampling in time
  rho_out  = interp1(time_ref',rho',z0dinput.cons.temps,'linear','extrap')';
  r0_out   = interp1(time_ref',r0',z0dinput.cons.temps,'linear','extrap')';
  z0_out   = interp1(time_ref',z0',z0dinput.cons.temps,'linear','extrap')';
  a_out    = interp1(time_ref',aminor',z0dinput.cons.temps,'linear','extrap')';
  K_out    = interp1(time_ref',elong',z0dinput.cons.temps,'linear','extrap')';
  d_out    = interp1(time_ref',delta',z0dinput.cons.temps,'linear','extrap')';
  ip_out   = interp1(time_ref',ip_ref',z0dinput.cons.temps,'linear','extrap');
  time_out = z0dinput.cons.temps;
 
case 'ip'

  rmin = min(min(z0dinput.geo.R - z0dinput.geo.a),min(r0 + min(rho .* cos(uni * ones(1,size(rho,2))),[],1)));
  rmax = max(max(z0dinput.geo.R + z0dinput.geo.a),max(r0 + max(rho .* cos(uni * ones(1,size(rho,2))),[],1)));

  
  if ip_ref(1) > z0dinput.cons.ip(1)
      % adding LCFS for initil current
      a_ini     = min(z0dinput.geo.a(1),aminor(1));
      R_ini     = z0dinput.geo.R(1);
      z0_ini    = z0dinput.geo.z0(1);
      K_ini     = min(z0dinput.geo.K(1),elong(1));
      d_ini     = min(z0dinput.geo.d(1),delta(1));
      ip_ini    = z0dinput.cons.ip(1);
      time_ini  = z0dinput.cons.temps(1);
      r0        = cat(2,R_ini,r0);
      z0        = cat(2,z0_ini,z0);
      aminor    = cat(2,a_ini,aminor);
      elong     = cat(2,K_ini,elong);
      delta     = cat(2,d_ini,delta);
      ip_ref    = cat(2,ip_ini,ip_ref);
      time_ref  = cat(2,time_ini,time_ref);
      rho       = cat(2,a_ini*ones(size(rho,1),1),rho);    
  end
  if ip_ref(end) > z0dinput.cons.ip(end)
      % adding LCFS for initil current
      a_ini     = min(z0dinput.geo.a(end),aminor(end));
      R_ini     = z0dinput.geo.R(end);
      z0_ini    = z0dinput.geo.z0(end);
      K_ini     = min(z0dinput.geo.R(end),elong(end));
      d_ini     = min(z0dinput.geo.z0(end),delta(end));
      ip_ini    = z0dinput.cons.ip(end);
      time_ini  = z0dinput.cons.temps(end);
      r0        = cat(2,r0,R_ini);
      z0        = cat(2,z0,z0_ini);
      aminor    = cat(2,aminor,a_ini);
      elong     = cat(2,elong,K_ini);
      delta     = cat(2,delta,d_ini);
      ip_ref    = cat(2,ip_ref,ip_ini);
      time_ref  = cat(2,time_ref,time_ini);
      rho       = cat(2,rho,a_ini*ones(size(rho,1),1));    
  end
  % resampling on ip
  % make monotone variable for resampling
  dipin  = abs(diff(z0dinput.cons.ip));
  dipout = abs(diff(ip_ref));
  dipmin = min(min(dipin(dipin > eps)),min(dipin(dipout > eps)));
  ipinm  = z0dinput.cons.ip(1) + cat(1,0,cumsum(max(dipmin,abs(diff(z0dinput.cons.ip)))));
  ipoutm = ip_ref(1) + cat(1,0,cumsum(max(dipmin,abs(diff(ip_ref(:))))));
  %  
  rho_out   = interp1(ipoutm,rho',ipinm,'linear','extrap')';
  r0_out    = interp1(ipoutm,r0',ipinm,'linear','extrap')';
  z0_out    = interp1(ipoutm,z0',ipinm,'linear','extrap')';
  a_out     = interp1(ipoutm,aminor',ipinm,'linear','extrap')';
  K_out     = interp1(ipoutm,elong',ipinm,'linear','extrap')';
  d_out     = interp1(ipoutm,delta',ipinm,'linear','extrap')';
  time_out  = interp1(ipoutm,time_ref',ipinm,'linear','extrap');
  ip_out    = z0dinput.cons.ip;
  % prevent small time step
  dt_min = min(min(abs(diff(z0dinput.cons.temps))),abs(min(diff(time_ref))));
  for k = 2:length(time_out)
    if diff(time_out(k-1:k)) < dt_min
        time_out(k) =  time_out(k -1) + dt_min;       
    end 
  end
  % security to keep plasma in vacuum chamber
  uni2d = uni * ones(1,size(rho_out,2));
  rho_test = sqrt(max(cos(uni2d) .* rho_out,rmin - ones(size(uni)) * r0_out) .^ 2 + (sin(uni2d) .* rho_out) .^ 2);
  rho_test = sqrt(min(cos(uni2d) .* rho_test,rmax - ones(size(uni)) * r0_out) .^ 2 + (sin(uni2d) .* rho_test) .^ 2);
  frho     = ones(size(rho_out,1),1) * min(rho_test ./ rho_out,[],1); 
  rho_out = frho .* rho_out;  
  a_out   = frho(1,:) .* a_out;
  K_out   = frho(1,:) .* K_out;
  d_out   = frho(1,:) .* d_out;
  
otherwise
  % no resampling
  if ~same(time_ref,z0dinput.cons.temps) || ~same(ip_ref,z0dinput.cons.ip)
	error('data is not compatible to be used without resampling')
  else
      rho_out   = rho;
      r0_out    = r0;
      z0_out    = z0;
      a_out     = aminor;
      K_out     = elong;
      d_out     = delta;
      time_out  = time_ref';
      ip_out    = z0dinput.cons.ip';  
  end 
end

% new LCFS for METIS
sepa.R  = (ones(size(uni)) * r0_out + rho_out .* cos(uni * ones(1,size(rho_out,2))))';
sepa.Z  = (rho_out .* sin(uni * ones(1,size(rho_out,2))))';

% correction of change of height of most external point (must be 0 at this
% stage as it is the definition of z0).
% Rsepa_max = max(sepa.R,[],2);
% maskrmax  = (sepa.R == Rsepa_max(:,ones(1,size(sepa.R,2))));
% z0_error  = sum(sepa.Z .* maskrmax,2) ./ sum(maskrmax,2);
% sepa.Z      = sepa.Z - z0_error(:,ones(1,size(sepa.Z,2)));


% mise a jour des moments
geo.a    = a_out(:);
geo.R    = r0_out(:);
geo.z0   = z0_out(:);
geo.K    = K_out(:);
geo.d    = d_out(:);
% conservation de la raideur magnetique
sepa.b0 = z0dinput.geo.b0 .* z0dinput.geo.R ./ geo.R; 



z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
z0dinput.geo.z0      = geo.z0;     % centre geometrique du plasma en Z (m)
z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
z0dinput.geo.K       = geo.K;     % elongation (b/a)
z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
z0dinput.geo.b0      = sepa.b0;  % champ toroidal a R0

z0dinput.cons.temps = time_out(:);
z0dinput.cons.ip = ip_out(:);

z0dinput.exp0d.temps = time_out(:);
z0dinput.exp0d.ip    = ip_out(:);
z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
z0dinput.exp0d.Zsepa = sepa.Z;       % vecteur Z des points de la separtrice (m)


if plotonoff == 0
  return
end

t  = asin(max(0,min(1,z0dinput.geo.d)));
u  = linspace(0,2.*pi,nbtheta);
vu = ones(size(u));
Rtest  = z0dinput.geo.R *vu + (z0dinput.geo.a * vu) .* cos(ones(size(geo.R,1),1) * u + t * sin(u));
Ztest = (z0dinput.geo.a .* z0dinput.geo.K) * sin(u);

h = findobj(0,'type','figure','tag','z0geosepa');
if isempty(h)
h=figure('tag','z0geosepa');
else
figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
zplotprof(gca,z0dinput.cons.temps,Rtest,Ztest+z0dinput.geo.z0 * vu,'color','b','marker','.','linestyle','none');
zplotprof(gca,z0dinput.cons.temps,z0dinput.exp0d.Rsepa,z0dinput.exp0d.Zsepa+z0dinput.geo.z0 * vu,'color','r','marker','none','linestyle','-');
axis('square')
axis('equal')


% figure
% time_ref = NaN .* ones(length(fbe_lcfs),1);
% for k= 1:length(fbe_lcfs)
% 	time_ref(k) = fbe_lcfs{k}.time;
% end
% 
% for k= 1:length(z0dinput.cons.temps)
%     l = find(time_ref >= z0dinput.cons.temps(k),1);
%     if isempty(l)
%         l = length(time_ref);
%     end
%     r = fbe_lcfs{l}.R;
%     z = fbe_lcfs{l}.Z;
%     plot(z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:)+z0dinput.geo.z0(k),'b',r,z,'r.')
%     pause
% end
% 
% 
% keyboard


