%-------------------------------------------------------------------------------  
%  
function [data,gdata]=zgfile_jt60sa(file,modesepa)


if (nargin < 1) ||  (isempty(file))
  
  error('[data,gdata]=zgfile_jt60sa(file)');

end
if (nargin < 2) ||  isempty(modesepa)
    modesepa = 0;
end
if isletter(file(end))
    file(end)=[];
end
ind = find(isletter(file));
if sum(ind) > 1
  file=file(ind(1):ind(2)-1);
end
ind = find(isstrprop(file,'cntrl'));
if ~isempty(ind)
    file=file(1:ind(1)-1);
end
data.temps = str2num(file(end-5:end))/1e3;
[gdata,ireadok] = read_gfile_jt60sa(file);
[rr,zz] = meshgrid(gdata.rgefit,gdata.zgefit);
cc=contour(rr,zz,gdata.psirz,gdata.psibnd);

ntot = length(cc);
nco  = 1;
Rcl = {};
Zcl ={};

while nco < ntot

  nc   = cc(2,nco);
  ndeb = nco+1;
  nfin = ndeb+nc-1;
  Rcl{end+1} = cc(1,ndeb:nfin);
  Zcl{end+1} = cc(2,ndeb:nfin);

  nco = nco+nc+1;
  
end

% select the right contour
flag = zeros(length(Rcl),1);
for k=1:length(Rcl)
    flag(k) = inpolygon(gdata.rmaxis,gdata.zmaxis,Rcl{k},Zcl{k});
end
Rc = [];
Zc = [];
for  k=1:length(flag)
   if flag(k) && (length(Rcl{k}) > length(Rc))
       Rc = Rcl{k};
       Zc = Zcl{k};
   end
end

qpsi = gdata.qpsi;
xpsi = linspace(0,1,length(qpsi));
rho  = cumtrapz(xpsi,qpsi);
xphi = sqrt(rho/max(rho));


[R,Z]        = fitsepa(Rc,Zc,201);
figure;plot(Rc,Zc,'r',R,Z,'k');
if modesepa == 1
    Rc = R;
    Zc = Z;
elseif (Rc(1) ~= Rc(end)) || (Zc(1) ~= Zc(end))
    Rc(end+1) = Rc(1);
    Zc(end+1) = Zc(1);
end

data.geo.R   = Rc;
data.geo.Z   = Zc;
data.geo.b0  = gdata.bzero;
data.geo.a   = (max(data.geo.R(:))-min(data.geo.R(:)))/2;
data.geo.r0  = (max(data.geo.R(:))+min(data.geo.R(:)))/2;
data.ip = gdata.cpasma;
data.prof.jmoy(1,:) = interp1(xphi,jmoy,param.gene.x);
data.prof.ptot(1,:) = interp1(xphi,gdata.pres,param.gene.x);

data.geo.b0 = gdata.bzero;
data.psibnd = gdata.psibnd;
data.psimag = gdata.psimag;

keyboard
