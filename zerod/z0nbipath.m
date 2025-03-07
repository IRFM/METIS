% calcul de l'bsorbtion de NBI pour le zerod
% lm1 = sigma * ni
% zout = sqrt(PHI/pi/B0) normalized
% pvout = power density (arbitrary unit)
% pitch_out = cos(beam_direction,toroidal_direction)
% shinethrough = fraction of power losses
function [zout,pvout,shinethrough,pitch_out] = z0nbipath(z,lm1,R,a,d0,rtang,angle_nbi,vol,Raxe,dext)

% alongement (couche de plus a bord)
nbp  = size(z,2)+1;
dvol = vol(:,end) - vol(:,end-1);
vol  = cat(2,vol,vol(:,end) + dvol);
dz   = z(2) - z(1);
z    = cat(2,z,z(end) + dz);
if ~isempty(Raxe)
	dRaxe = Raxe(:,end) - Raxe(:,end-1);
	Raxe =cat(2,Raxe,Raxe(:,end) + dRaxe);
	Raxe    = cat(2,Raxe(:,end:-1:2),Raxe);
end

lm1  = cat(2,lm1,1e-38 .* ones(size(lm1,1),1));
lm1  = cat(2,lm1(:,end:-1:2),lm1);
lm1(:,end) = 1e38;

z      = cat(2,z(end:-1:2),-z);
vol    = cat(2,vol(:,end:-1:2),vol);

ve = ones(1,size(z,2));
vt = ones(size(R,1),1);

%Rext = (R + a.* max(z(:))) * ve;
Rext = (R + a.* max(z(:)));
Rext = Rext(:,ones(size(ve)));
if (rtang > 0) && (angle_nbi ~= 0)
	v = min(pi/2 - eps,abs(asin(rtang ./ Rext)));
else
	v = abs(angle_nbi ./ 180 .* pi) * ones(size(Rext));
end
ux  = 1 - z .^ 2;

% injetcion off axis
%az  = max(0.1 .* (a * ve),sqrt((a * ve) .^ 2 - (((vt .* Zext) * ve) .^ 2))); 
%az   = a * ve;
az   = a(:,ones(size(ve)));
% rayon des surface de flux
if isempty(Raxe)
	%Rc  = R * ve + az .* (vt * z)  + d0 * ux;
	Rc  = R * ve + az .* z(ones(size(vt)),:)  + d0 * ux;
else
	%Rc  = Raxe + az .* (vt * z);
	Rc  = Raxe + az .* z(ones(size(vt)),:);
end

% intersection cercle et droite
c2 = ones(size(Rc));
c1 = 2 .* (Rext .* cos(v) - dext .* sin(v));
c0 = Rext .^ 2 + dext .^ 2 - Rc .^ 2;
[tp,tm]= zpolyroot(c0,c1,c2);
xm = Rext  + tm .* cos(v);
ym = dext  - tm .* sin(v);
xp = Rext  + tp .* cos(v);
yp = dext  - tp .* sin(v);

% concanetation
zz   = cat(2,z,z);
vvol = cat(2,vol,vol);
lm1   = cat(2,lm1,lm1);
vee  = ones(1,size(zz,2));
xx   = cat(2,xp,xm);
yy   = cat(2,yp,ym);
rr   = sqrt(xx .^2 +yy .^2);
Rext   = cat(2,Rext,Rext);
%lm1(rr < cat(2,(R *ve - az),(R *ve - az))) = 1e38;
lm1(rr < cat(2,(R(:,ones(size(ve))) - az),(R(:,ones(size(ve))) - az))) = 1e38;
% indice de couche
lpath     = sqrt((xx- Rext) .^2 + (yy - dext) .^2);
% pitch angle
pitch_angle   = abs((xx - Rext) .* yy  - (yy - dext) .* xx) ./ max(eps,lpath) ./ max(eps,rr); 
%zzl  = vt * zz;
zzl  = zz(ones(size(vt)),:);

% selcetion des point dans le plasma
xx(~isfinite(xx)) = sqrt(-1);
yy(~isfinite(yy)) = sqrt(-1);
mask = ~((imag(xx)) | (imag(yy)));
lpath(find(~mask)) = 1e308;

% tri des donnees
[lpath,iks] = sort(lpath,2,'ascend');
for l = 1:size(iks,1)
	iqz            = iks(l,:);
        pitch_angle(l,:) = pitch_angle(l,iqz);
	zzl(l,:)       = zzl(l,iqz);
	vvol(l,:)      = vvol(l,iqz);
	xx(l,:)        = xx(l,iqz);
	yy(l,:)        = yy(l,iqz);
	lm1(l,:)       = lm1(l,iqz);
	masks(l,:)     = mask(l,iqz);
end

% longueur de corde
maskd  = min(masks(:,1:(end-1)),masks(:,2:end));
dl   = diff(lpath,1,2);
dl(find(~maskd)) = 0;
% facteur d'attenuation de la couche
dinl  = exp(- (lm1(:,1:(end-1)) + lm1(:,2:end)) .* dl ./ 2);
dinl(find(~maskd)) = 0;
inl   = cat(2,ones(size(dinl,1),1),cumprod(dinl,2));
ll    = cat(2,zeros(size(dl,1),1),cumsum(dl,2));
pvi   = - diff(inl,1,2) ./ max(eps,abs(diff(vvol,1,2)));
pvi(find(~maskd)) = 0;
pvi(diff(vvol,1,2) == 0) = 0;
pvi   = cat(2,zeros(size(inl,1),1),pvi);

% pour le cacul des pertes
inol = inl;
inol(inol <=  eps) = Inf;
shinethrough = min(inol,[],2);


% tri final
[zzl,iko] = sort(zzl,2,'descend');
for l = 1:size(iko,1)
	iqz            = iko(l,:);
	pitch_angle(l,:) = pitch_angle(l,iqz);
	inol(l,:)       = inol(l,iqz);
	pvi(l,:)        = pvi(l,iqz);
end

% mise en forme
zzl     = (zzl(:,1:2:end-1) +zzl(:,2:2:end))./2;
pitch_angle = (pvi(:,1:2:end-1) .* pitch_angle(:,1:2:end-1) + pvi(:,2:2:end) .* pitch_angle(:,2:2:end));
pvi     = (pvi(:,1:2:end-1) +pvi(:,2:2:end));
pitch_angle = pitch_angle ./ max(eps,pvi);
pitch_angle(find(pvi <= eps)) = 0;

zzf     = zzl(:,1:nbp);
pvout  = pvi(:,1:nbp);
pvout(:,1:(nbp-1)) = pvout(:,1:(nbp-1)) + pvi(:,end:-1:(nbp+1));
pitch_out  = pitch_angle(:,1:nbp) .* pvi(:,1:nbp);
pitch_out(:,1:(nbp-1)) = pitch_out(:,1:(nbp-1)) + pitch_angle(:,end:-1:(nbp+1)) .* pvi(:,end:-1:(nbp+1));
pitch_out  = pitch_out ./ max(eps,pvout);
pitch_out(find(pvout == 0)) = 0;

% creation des donnees de sorties
pvout = pvout(:,end:-1:2);
zout  = zzf(:,end:-1:2);
pitch_out  = pitch_out(:,end:-1:2);

pitch_out(~isfinite(pvout) | ~isfinite(pvout)) = 0;
pvout(~isfinite(pvout)) = 0;


if 1>2

	figure(26);clf
	s = linspace(0,2*pi);
	for k=1:size(Rc,2)
		plot(mean(Rc(:,k)) .* cos(s),mean(Rc(:,k)) .* sin(s),'b');
		hold on
	end
	lx = linspace(-Rext(1),Rext(1));
	plot(Rext(1)+lx .* cos(mean(v(:))),dext(1) - lx.*sin(mean(v(:))),'r');
	plot(xx(find(masks)),yy(find(masks)),'+g');
	plot(Rext(1),dext,'ok');
	hold off
	drawnow
	
end

