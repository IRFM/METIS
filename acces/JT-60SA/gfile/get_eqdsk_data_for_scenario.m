function [time,ip,R_LCFS,Z_LCFS,RB0,Rsepa,Zsepa,R0,Z0,aminor,K,d,psi,psin,q,ptot,fdia,pprim,ffprim] = get_eqdsk_data_for_scenario(dirname)

cd_mem = pwd;
cd(dirname);
rep = dir('.');
name_list = {};
time = [];
for k=1:length(rep)
    switch rep(k).name
        case {'.','..'}
           % rien     
        otherwise
          name_list{end+1} = rep(k).name;
          time(end+1) = str2num(rep(k).name) ./ 1e3;
    end
end
[time,index] = sort(time);
name_list    = name_list(index);
time         = time(:);
% read daat from gfiles
ip = NaN * ones(size(time));
RB0 = NaN * ones(size(time));
R0 = NaN * ones(size(time));
Z0 = NaN * ones(size(time));
aminor= NaN * ones(size(time));
K = NaN * ones(size(time));
d = NaN * ones(size(time));
nbx = 129;
psi    = NaN * ones(length(time),nbx);
q      = NaN * ones(length(time),nbx);
pprim  = NaN * ones(length(time),nbx);
ffprim = NaN * ones(length(time),nbx);
fdia   = NaN * ones(length(time),nbx);
ptot   = NaN * ones(length(time),nbx);
nbp = 201;
Rsepa = NaN * ones(length(time),nbp);
Zsepa = NaN * ones(length(time),nbp);
R_LCFS = {};
Z_LCFS = {};
for k=1:length(time)
    [gdata,ireadok] = read_gfile_jt60sa(name_list{k});
    R_LCFS{k}       = gdata.rbbbs;
    Z_LCFS{k}       = gdata.zbbbs;
    ip(k)           = gdata.cpasma;
    RB0(k)          = gdata.bzero .* gdata.rzero;
    psi(k,:)        = linspace(gdata.psimag,gdata.psibnd,nbx);
    q(k,:)          = gdata.qpsi;
    pprim(k,:)      = gdata.pprime;
    ffprim(k,:)     = gdata.ffprim;
    ptot(k,:)       = gdata.pres;
    fdia(k,:)       = gdata.fpol;
    % put all LCFS on same number of points
    [Rsepa(k,:),Zsepa(k,:),aminor(k),R0(k),Z0(k),K(k),d(k)] = reshape_LCFS(R_LCFS{k},Z_LCFS{k},nbp);
end
psin = ones(size(time)) * linspace(0,1,nbx);
cd(cd_mem);




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

