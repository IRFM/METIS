
function [temps_out,Rout,Zout,Rsepa,Zsepa]=z0make2dsurface(post,temps)



% pas
profli = post.profil0d;
zs     = post.zerod;
geo    = post.z0dinput.geo;
if ~isfield(geo,'z0')
	disp('undefined z0 : z0 set to 0');
	geo.z0 = 0.*geo.a;
end

% (R,Z) separatrice
if isfield(post.z0dinput.exp0d,'Rsepa') & isfield(post.z0dinput.exp0d,'Zsepa')
	Rs = post.z0dinput.exp0d.Rsepa;
	Zs = post.z0dinput.exp0d.Zsepa +  geo.z0 * ones(1,size(Rs,2));
        maskrmax  = (Rs == (max(Rs,[],2) * ones(1,size(Rs,2))));
        z0   = sum(Zs .* maskrmax,2) ./ sum(maskrmax,2);
else
	td  = asin(max(0,min(1,geo.d)));
	u   = linspace(0,2.*pi,201);
	vu  = ones(size(u));
	vt  = ones(size(td));
	Rs  = geo.R *vu + (geo.a * vu) .* cos(vt * u + td * sin(u));
	Zs  = (geo.a .* geo.K) * sin(u) + geo.z0 * ones(1,size(Rs,2));
	z0  = geo.z0;
end




% echelle
Rmin = min(Rs(:)) .* 0.95;
Rmax = max(Rs(:)) .* 1.05;
Zmin = min(Zs(:)) .* 1.05;
Zmax = max(Zs(:)) .* 1.05;


% si 2 entrees
if nargin < 2
	% nombre de frames
	nbf = length(profli.temps);
	nbini = 1;
	temps = profli.temps;
else
	nbini = find(profli.temps >= min(temps),1);
	nbf   = max(find(profli.temps <= max(temps)));
        nbini = min(nbf,nbini);
end
llist = [2,4,6,8,10,12,14,16,18,20,21];
u    = linspace(0,2.*pi,35);


% reservasion memoire
nbt_eff   = length(nbini:nbf);
Rout      = NaN * ones(nbt_eff,length(llist),length(u)+1);
Zout      = NaN * ones(nbt_eff,length(llist),length(u)+1);
temps_out = NaN * ones(nbt_eff,1);
Rsepa     = NaN * ones(nbt_eff,size(Rs,2));
Zsepa     = NaN * ones(nbt_eff,size(Zs,2));

% boucle sur les temps
first = 1;
count = 1;
fprintf('2D:');
for kf =nbini:nbf
  fprintf('.');
  % temps courant
  tc = profli.temps(kf);
  k  = min(find(tc <= post.zerod.temps));
  if isempty(k)
  	k= length(post.zerod.temps);
  end

   %clear Rfond Zfond
  for ll = 1:length(llist)
  	l =llist(ll);

	td  = asin(max(0,min(1,geo.d(k) .* (profli.xli(l) .^ 3))));
	R   = profli.Raxe(kf,l) + geo.a(k) .* profli.xli(l).* cos(u + td * sin(u));
	Z   = geo.a(k) .* profli.xli(l) .* profli.kx(kf,l) * sin(u)+ z0(k) ;

	% deformation
	cx   = (R - profli.Raxe(kf,l)) + sqrt(-1) .* (Z - z0(k));
	thx  = unwrap(angle(cx),[],2);
	rhox = abs(cx);
	thx(thx<0) = thx(thx<0) + 2 .* pi;
	[thx,indx] = sort(thx);
	rhox       = rhox(indx);

	% separtrice a ce temps version analytique
	tdl  = asin(max(0,min(1,geo.d(k))));
	Rl   = profli.Raxe(kf,end) + geo.a(k) .* cos(u + tdl * sin(u));
	Zl   = geo.a(k).* profli.kx(kf,end) * sin(u)+ z0(k);
	cl   = (Rl - profli.Raxe(kf,l)) + sqrt(-1) .* (Zl - z0(k));
	thl  = unwrap(angle(cl),[],2);
	rhol = abs(cl);
	thl(thl<0) = thx(thl<0) + 2 .* pi;
	[thl,indl] = sort(thl);
	rhol       = rhol(indl);
	rhol = cat(2,rhol,rhol,rhol);
	thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
	indnok = find(diff(thl)<=0);
	thl(indnok) =[];
	rhol(indnok)  = [];

	% separtrice a ce temps complete
	cc   = (Rs(k,:) - profli.Raxe(kf,l)) + sqrt(-1) .* (Zs(k,:) - z0(k));
	thc  = unwrap(angle(cc),[],2);
	rhoc = abs(cc);
	thc(thc<0) = thc(thc<0) + 2 .* pi;
	[thc,indc] = sort(thc);
	rhoc       = rhoc(indc);
	rhoc = cat(2,rhoc,rhoc,rhoc);
	thc = cat(2,thc -2.*pi,thc,thc+2.*pi);
	indnok = find(diff(thc)<=0);
	thc(indnok) =[];
	rhoc(indnok)  = [];


	% regle de trois, mise a l'echelle
	if post.z0dinput.option.morphing > 0
		morf    = (1 - profli.xli(l).^ post.z0dinput.option.morphing) +  ...
	          		profli.xli(l) .^ post.z0dinput.option.morphing .* pchip(thc',rhoc',thx')'  ./ pchip(thl',rhol',thx')';
	else
		morf = ones(size(thx));
	end
	rhox = rhox .* morf;
	Rx   = rhox .* cos(thx) + profli.Raxe(kf,l);
	Zx   = rhox .* sin(thx)+ z0(k);
	Rx(end+1) = Rx(1);
	Zx(end+1) = Zx(1);

        Rfond(ll,:) = Rx;
        Zfond(ll,:) = Zx;

  end
  
  temps_out(count) = tc;
  Rout(count,:,:) = shiftdim(Rfond,-1);
  Zout(count,:,:) = shiftdim(Zfond,-1);
  Rsepa(count,:)  = Rs(k,:);
  Zsepa(count,:)  = Zs(k,:);
  
  first = 0;
  count = count + 1;
end
fprintf('\n');
