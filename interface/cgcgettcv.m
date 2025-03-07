function [data,prof,equi]=cgcgettcvbug(numchoc,liste)

if nargin < 2
	error('two argument as input  ...');
end
if isempty(numchoc)
	error('shot number is needed')
end
if isempty(liste)
	error('no data asked')
elseif ~iscell(liste)
        liste = cellstr(liste);
end

try
 	cr = mdsconnect('frank.epfl.ch');
	if isempty(cr)
	  cr =0;
                             mdsopen('tcv_shot',numchoc);
	end

catch
	cr = -1;
end
if cr ~= 0
	disp('Probleme de connexion a TCV');
	return
end

   Ip = mdsvalue(eval(sprintf(['''\\',liste{1},''''])));
   tIp = mdsvalue(eval(sprintf(['''dim_of(\\',liste{1},')'''])));
   R0 = mdsvalue(eval(sprintf(['''\\',liste{2},''''])));
   Z0 = mdsvalue(eval(sprintf(['''\\',liste{3},''''])));
   tR0 = tIp;
   tZ0 = tIp;
   RB0 = mdsvalue(eval(sprintf(['''\\',liste{4},''''])));
   tRB0 = mdsvalue(eval(sprintf(['''dim_of(\\',liste{4},')'''])));
   li = mdsvalue(eval(sprintf(['''\\',liste{5},''''])));
   tli = mdsvalue(eval(sprintf(['''dim_of(\\',liste{5},')'''])));
   vl = mdsvalue(eval(sprintf(['''\\',liste{6},''''])));
   tvl = mdsvalue(eval(sprintf(['''dim_of(\\',liste{6},')'''])));
   bp = mdsvalue(eval(sprintf(['''\\',liste{7},''''])));
   tbp = mdsvalue(eval(sprintf(['''dim_of(\\',liste{7},')'''])));
   jtot = mdsvalue(eval(sprintf(['''\\',liste{8},''''])));
   rjtot = mdsvalue(eval(sprintf(['''dim_of(\\',liste{8},'0)'''])));
   zjtot = mdsvalue(eval(sprintf(['''dim_of(\\',liste{8},'1)'''])));
   qpsi = mdsvalue(eval(sprintf(['''\\',liste{9},''''])));
   tqpsi  =mdsvalue(eval(sprintf(['''dim_of(\\',liste{9},')'''])));
   rc = mdsvalue(eval(sprintf(['''\\',liste{10},''''])));
   trc  =tIp;
   zc = mdsvalue(eval(sprintf(['''\\',liste{11},''''])));
   tzc  =tIp;
%
% profil precis
%  
   timeth = mdsvalue(eval(sprintf(['''\\',liste{12},''''])));
   if ~isempty(timeth)
     Teth =  mdsvalue(eval(sprintf(['''\\',liste{13},''''])));
     neth =  mdsvalue(eval(sprintf(['''\\',liste{14},''''])));
     psiteth =  mdsvalue(eval(sprintf(['''\\',liste{15},''''])));

     Te = mdsvalue(eval(sprintf(['''\\',liste{16},''''])));
     rho = mdsvalue(eval(sprintf(['''dim_of(\\',liste{16},')'''])));
     tTe = mdsvalue(eval(sprintf(['''\\',liste{17},''''])));
     if ischar(tTe)
       tTe = timeth';
     else
       tTe=tTe(1,:);
     end
     ne = mdsvalue(eval(sprintf(['''\\',liste{18},''''])));
     tne  =tTe;
   else
     disp('no avialable data from Thomson for this shot')
     return
   end
%
% Puissance ECRH injectee
%
   freqecrh = mdsvalue(eval(sprintf(['''\\',liste{19},''''])));
   phiecrh = mdsvalue(eval(sprintf(['''\\',liste{20},''''])));
   thetaecrh = mdsvalue(eval(sprintf(['''\\',liste{21},''''])));
   Pecrh = mdsvalue(eval(sprintf(['''\\',liste{22},''''])));
   tPecrh  =mdsvalue(eval(sprintf(['''dim_of(\\',liste{22},')'''])));
%
% flux poloidal
%
  psi = mdsvalue(eval(sprintf(['''\\',liste{23},''''])));
  rpsi  =mdsvalue(eval(sprintf(['''dim_of(\\',liste{23},',0)'''])));
  zpsi  =mdsvalue(eval(sprintf(['''dim_of(\\',liste{23},',1)'''])));
  tpsi = tIp;
%
% donnees geometriques
%
  kappa95 = mdsvalue(eval(sprintf(['''\\',liste{24},''''])));
  kappaedge = mdsvalue(eval(sprintf(['''\\',liste{25},''''])));
  delta95 = mdsvalue(eval(sprintf(['''\\',liste{26},''''])));
  deltaedge = mdsvalue(eval(sprintf(['''\\',liste{27},''''])));
  tkappa95 = tIp;
%
% impuretes
%
  zeffm = mdsvalue(eval(sprintf(['''\\',liste{28},''''])));
  tzeffm   = mdsvalue(eval(sprintf(['''dim_of(\\',liste{28},')'''])));
  zeffmerr = mdsvalue(eval(sprintf(['''\\',liste{29},''''])));
  zeff    =  mdsvalue(eval(sprintf(['''\\',liste{30},''''])));
  if ~ischar(zeff)
    tzeff   = mdsvalue(eval(sprintf(['''dim_of(\\',liste{30},')'''])));
    zzeff  = mdsvalue(eval(sprintf(['''\\',liste{31},''''])));
  else
  
    zeff = zeros(length(tIp),10);
    tzeff = tIp;
    zzeff = zeff;

  end
%
% psitoolbox
%
  rs = mdsvalue(eval(sprintf(['''\\',liste{32},''''])));
  zs = mdsvalue(eval(sprintf(['''\\',liste{33},''''])));
  rhotbx = linspace(0,1,41);
  ttbx = mdsvalue(eval(sprintf(['''dim_of(\\',liste{33},',2)'''])));
  tetatbx = mdsvalue(eval(sprintf(['''dim_of(\\',liste{33},',1)'''])));
  vol = mdsvalue(eval(sprintf(['''\\',liste{34},''''])));
  wtot = mdsvalue(eval(sprintf(['''\\',liste{35},''''])));
  twtot = tIp;
  psia = mdsvalue(eval(sprintf(['''\\',liste{36},''''])));
  tpsia = tIp;

%
% mise en forme donnees zeff
%


%
% mise en forme donnees TCV ECCD
%
if ~ischar(tPecrh)
  data.tPec = tPecrh;
   for k=1:9
       ind =find(~isnan(freqecrh(k,:)));
       if ~isempty(ind)
         data.ecfreq(k) = mean(freqecrh(k,ind));
         data.ecphi(k)   = mean(phiecrh(k,ind));
         data.ectheta(k)   = mean(thetaecrh(k,ind));
         data.Pec(:,k)     = Pecrh(:,k)'*1e3;
         indnan = find(isnan(data.Pec(:,k)));
         data.Pec(indnan,k)=0;
        else
          data.ecfreq(k) = nan;
          data.ecphi(k) = nan;
          data.ectheta(k) = nan;
         data.Pec(:,k)     = zeros(size(tPecrh));
        end
   end
else
  disp(' no data avialable for ECRH')
  data.tPec = tIp;
  data.Pec = zeros(length(tIp),9);
  data.ecfreq = zeros(1,9);
  data.ecphi = zeros(1,9);
  data.ectheta = zeros(1,9);
end
%
% mis en forme separatrice
%
  rc(isnan(rc))=0;
  zc(isnan(zc))=0;
   for k=1:size(rc,2)
      indzero = min(find(rc(:,k)==0));
      if ~isempty(indzero)
        indz(k) = indzero-1;
     else
       indz(k)= size(rc,1);
     end
   end
%
% meme nombre de points pour la separatrice
%
   nteta = linspace(0,2*pi,250)';
  for k=1:size(rc,2)
    rca =rc(1:indz(k),k);
    zca =zc(1:indz(k),k);
    rc0=mean(rca);
    zc0=mean(zca);
    teta = unwrap(angle((rca-rc0)+sqrt(-1)*(zca-zc0)));
    if min(teta) <0
      teta = teta-min(teta);
    end
    indp = find(diff(teta) == 0);
    if ~isempty(indp)
       teta(indp) = [];
       rca(indp) = [];
       zca(indp) = [];
    end
    nrc(k,:) = interp1(teta,rca,nteta,'linear')';
    nzc(k,:) = interp1(teta,zca,nteta,'linear')';
%
% calcul du petit rayon
%
    amin(k) = (max(rca)-min(rca))/2;
  end
%
% profil sur 101 points radiaux + elimination des temps avec NaN
%
xtcv = linspace(0,1,101);
nc=1;
clear Tetcv netcv tTetcv
for k=1:size(tTe,2)
    indnan=find(~isnan(Te(:,k)));
    if ~isempty(indnan) & ~isnan(tTe(k))
      Tetcv(nc,:) = interp1(rho(indnan)',...
                                           Te(indnan,k)',xtcv,'spline');
      netcv(nc,:) = interp1(rho(indnan)',...
                                           ne(indnan,k)',xtcv,'spline');
       tTetcv(nc)=tTe(k);
      nc = nc+1;
     end
end


data.ip = Ip;
data.tip = tIp;
data.rb0= RB0;
data.trb0 = tRB0;
data.li = li;
data.tli = tli;
data.vl = vl;
data.tvl = tvl;
data.vl = bp;
data.tvl = tbp;
data.j = jtot;
data.tj = jtot;
data.qpsi = qpsi;
data.tqpsi = tqpsi;
data.rext = nrc;
data.trext = trc;
data.zext = nzc;
data.tzext = tzc;
data.amin = amin;
data.zeffm = zeffm;
data.tzeffm = tzeffm
data.zeffmerr = zeffmerr;
prof.Te = Tetcv;
prof.ne = netcv;
prof.tTe = tTetcv;
prof.rhotcv = xtcv;

equi.a = interp1(trc,data.amin,ttbx);
equi.t = ttbx;
equi.R0 = (max(nrc,[],2)+min(nrc,[],2))/2;
equi.B0 = interp1(data.trb0,data.rb0,trc)./equi.R0;
equi.R0=interp1(trc,equi.R0,ttbx);
equi.e = (max(nzc,[],2)-min(nzc,[],2))./(max(nrc,[],2)-min(nrc,[],2));
equi.e=interp1(trc,equi.e,ttbx);

equi.B0=interp1(trc,equi.B0,ttbx);
equi.x  = rhotbx;
equi.vol = vol';
equi.tpbx = ttbx;
npsi=interp1(tpsi,shiftdim(psi,2),ttbx);
[volloc,psiloc]=zpsitcv(rpsi,zpsi,shiftdim(npsi,1));
equi.psi = npsi;
equi.rpsi = rpsi;
equi.zpsi=zpsi;
equi.volloc=volloc;
equi.psiloc=psiloc;
equi.kappa95 = kappa95;
equi.tkappa95 = tkappa95;
equi.wdia = wtot;
equi.twdia = twtot;
equi.psia = psia;
equi.tpsia = tpsia;

mdsclose(numchoc);

