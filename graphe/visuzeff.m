% visualisation du Zeff
%
% syntaxe : visuzeff({choc,Z,zimp,publi})
%  
% entrees :
%
%    choc = numero du choc (defaut = dernier choc)
%    Z    = charge du gaz principal (defaut = 1)
%    zimp = charge de l'impurete principale (defaut = 6)
%    publi = eppaisseur des lignes (defaut = 0.3)
%
function visuzeff(choc,Z,zimp,publi)

if nargin < 1
	choc = tsdernier_choc;
end

if nargin < 2
   Z = 1;
end
if nargin < 3
	zimp = 6;
elseif isempty(zimp)
	zimp = 6;
end
if nargin < 4 
   publi = 0.5;
elseif isempty(publi)
   publi = 3;
end


[zeff,tzeff]               = tsbase(choc,'szfbrm');
[zeff_tg,tzeff_tg]         = tsbase(choc,'szfbrmtan');
[ip,tip]        = tsbase(choc,'sipmes');
[r0,tr0]        = tsbase(choc,'srmaj');
[a,ta]            = tsbase (choc,'samin');
[gnl,tnl]       = tsbase(choc,'gnl');
nl                 = max(gnl,[],2);
[prad,tprad] = tsbase(choc,'sprad');
if isempty(prad)
	prad = NaN .* tip;
	tprad = tip;
end

[itor,titor ]   = tsbase(choc,'sitor');
[ptot,tptot ]   = tsbase(choc,'gbilan%10');
[pohm,tpohm ]   = tsbase(choc,'gbilan%1');
[pfci,tfci ]   = tsbase(choc,'spuiss');
if isempty(pfci)
	pfci = NaN .* tip;
	tfci = tip;
end
[plh,tlh ]   = tsbase(choc,'GHYB%3');
if isempty(plh)
	plh = NaN .* tip;
	tlh = tip;
end
[qa,tqa,c] =  tsbase(choc,'sqpsi');


% largeur de la couche rayonnante
[centre,largeur,cprad] = zbolo(choc,0);
if ~isfinite(largeur)
	largeur = 0.1;
end


if ~isempty(tzeff)
   temps  = (0:0.1:max(max(tip(ip >0.1)),max(tzeff)))';
else
   temps  = (0:0.1:max(tip))';
end
r0        = interp1r(tr0,r0,temps,'nearest');
[ta,ind]       = sort(tr0);
a          = interp1r(ta,a,temps,'nearest');
if ~isempty(nl)
   nl        = interp1r(tnl,nl,temps,'nearest');
end
if ~isempty(prad)
   prad    = interp1r(tprad,prad,temps,'nearest');
end
itor     = interp1r(titor,itor,temps,'nearest');
rb0          = (4*pi*1e-7) .* 18 .* 2028 .* itor ./ 2 ./ pi;
ip        = interp1r(tip,ip,temps,'nearest');
ploss   = interp1r(tptot,ptot,temps,'nearest');
pohm   = interp1r(tptot,pohm,temps,'nearest');
qa = interp1r(tqa,qa,temps,'nearest');


if ~isempty(nl)
   nbar      = nl  ./ 2 ./ a ./ 1e20;
   Sp         = (4*pi^2) .* r0 .* a;
   if ~isempty(prad)
        %zeffrad   = Z + 7 .* prad ./ nbar .^ 2 ./ Sp; 
        zeffrad   = Z + 4.5 .* zimp .^ 0.12 .* prad ./ nbar .^ 1.89 ./ Sp .^ 0.94; 
        zeffrad_cor   = Z + 4.5 .* zimp .^ 0.12 .* prad ./ nbar .^ 1.89 ./ (Sp .* (1.1 - largeur)) .^ 0.94; 
  else
         zeffrad   = [];
	 zeffrad_cor = [];
  end
   if Z ==1 | Z == 2   
       zeffscl    = zeffscaling(nbar*10,ploss,ip,a,r0,itor,Z);
   else
       zeffscl1    = zeffscaling(nbar*10,ploss,ip,a,r0,itor,1);
       zeffscl2    = zeffscaling(nbar*10,ploss,ip,a,r0,itor,2);
       dz =Z - 1;
       zeffscl     = zeffscl1 .* (1 -dz) + zeffscl2 .* dz;
    end
    
 else
   zeffrad   = [];
   zeffscl    = [];
end
zeffmean = [];
zeffmin   = [];
zeffmax   = [];
nb = 0;


% zeff deduit du vloop
try
	sortie = vloop2zeff(choc,1);
	zeff_vloop = sortie.zeff;
	temps_zeff_vloop = sortie.temps;
catch
	zeff_vloop = [];
	temps_zeff_vloop = [];
end

% recherche dans les 30 dernier chocs
for k = 1:30
    [zeff_,tzeff_]  = tsbase(choc - k,'szfbrm');
    mem(k).zeff   = zeff_;
    mem(k).tzeff = tzeff_;
    if ~isempty(zeff_)
          if ~all(zeff_ == mean(zeff_))
                zz = interp1(tzeff_,medfilt1(zeff_,11),temps,'nearest');
                indnok = find(~isfinite(zz)|(zz < 1) | (zz > 6));
                indok   = find(isfinite(zz)& (zz >= 1) & (zz <= 6));
               if ~isempty(indok)
                  zz(indnok) = mean(zz(indok));
               end
               if isempty(indok)
               elseif isempty(zeffmean)
                     zeffmean = zz;
                     zeffmin   = zz;
                     zeffmax   = zz;
                     nb = 1;
               else
                     zeffmean = zeffmean + zz;
                     zeffmin   = min(zeffmin,zz);
                     zeffmax   = max(zz,zeffmax);
                     nb = nb + 1;
                end
        end
    end
    
end
if nb > 0
               zeffmean = zeffmean ./ nb;
 
else
      zeffmean  = NaN .* temps;
      zeffmin   = NaN .* temps; 
      zeffmax   = NaN .* temps; 
end
nbmem = nb;

% recherche dans les 30 dernier chocs (visee tangente)
zeffmean_tg  = [];
zeffmin_tg   = [];
zeffmax_tg   = [];
nb = 0;

for k = 1:30
    [zeff_,tzeff_]  = tsbase(choc - k,'szfbrmtan');
    mem(k).zeff_tg   = zeff_;
    mem(k).tzeff_tg = tzeff_;
    if ~isempty(zeff_)
          if ~all(zeff_ == mean(zeff_))
                zz = interp1(tzeff_,medfilt1(zeff_,11),temps,'nearest');
                indnok = find(~isfinite(zz)|(zz < 1) | (zz > 6));
                indok   = find(isfinite(zz)& (zz >= 1) & (zz <= 6));
               if ~isempty(indok)
                  zz(indnok) = mean(zz(indok));
               end
               if isempty(indok)
               elseif isempty(zeffmean_tg)
                     zeffmean_tg = zz;
                     zeffmin_tg   = zz;
                     zeffmax_tg   = zz;
                     nb = 1;
               else
                     zeffmean_tg = zeffmean_tg + zz;
                     zeffmin_tg   = min(zeffmin_tg,zz);
                     zeffmax_tg   = max(zz,zeffmax_tg);
                     nb = nb + 1;
                end
        end
    end
    
end
if nb > 0
               zeffmean_tg = zeffmean_tg ./ nb;
 
else
      zeffmean_tg = NaN .* temps;
      zeffmin_tg   = NaN .* temps; 
      zeffmax_tg   = NaN .* temps; 
end
nbmem_tg = nb;

hf1 = findobj(0,'type','figure','tag','visuzeff');
if isempty(hf1)
     hf1 = figure('tag','visuzeff','Name','Verification du Zeff');
else
    figure(hf1);
end
clf
set(hf1,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',publi,'color',[1 1 1])

if isempty(zeff_tg)
	kg = 2;
else
	kg = 4;
end
	

subplot(kg,1,1,'align')	
plot(temps,ip,temps,nbar*10,tfci,pfci,tlh,plh,temps,prad,temps,ploss);

hold off
legend('I_p (MA)','n_b_a_r (1e19)','P_F_C_I (MW)','P_L_H (MW)','P_r_a_d (MW)','P_i_n (MW)');
title(sprintf('Tore Supra # %d',choc));
%xlabel('temps (s)')
set(gca,'ylim',[0 Inf]);
grid on
	
subplot(kg,1,2,'align')	
leg = {};
lk        = 1;
if ~isempty(zeff)
   plot(tzeff,zeff,'r')
   leg{lk} = 'Zebre';
   lk   = lk +1;
end 
hold on
if ~isempty(zeffmean)
   plot(temps,zeffmean,'k',temps,zeffmax,'b',temps,zeffmin,'b')
   leg{lk} = sprintf('moy [%d chocs]',nbmem);
   lk   = lk +1;
   leg{lk} = sprintf('max [%d chocs]',nbmem);
   lk   = lk +1;
   leg{lk} = sprintf('min [%d chocs]',nbmem);
   lk   = lk +1;
end
if ~isempty(zeffscl)
   plot(temps,zeffscl,'g')
   leg{lk} = 'scaling TS';
   lk   = lk +1;
end
if~isempty(zeffrad)
   plot(temps,real(zeffrad),'c');
   leg{lk} = 'Rad (Matthews)';
   lk   = lk +1;
   plot(temps,zeffrad_cor,'m');
   leg{lk} = 'Rad corrige (Matthews)';
   lk   = lk +1;
end
if~isempty(zeff_vloop)
   plot(temps_zeff_vloop,real(zeff_vloop),'c.');
   leg{lk} = 'Zeff(Vloop) (Ohm)';
   lk   = lk +1;
end
 
hold off
legend(leg);
%title(sprintf('Tore Supra # %d',choc));
%xlabel('temps (s)')
ylabel('Zeff')   
set(gca,'ylim',[1 6])
grid on
%zoom xon

if kg > 2
	subplot(4,1,3,'align')	
	leg = {};
	lk        = 1;
	if ~isempty(zeff_tg)
	plot(tzeff_tg,zeff_tg,'r')
	leg{lk} = 'Zebre tangent';
	lk   = lk +1;
	end 
	hold on
	if ~isempty(zeffmean)
	plot(temps,zeffmean_tg,'k',temps,zeffmax_tg,'b',temps,zeffmin_tg,'b')
	leg{lk} = sprintf('moy [%d chocs]',nbmem_tg);
	lk   = lk +1;
	leg{lk} = sprintf('max [%d chocs]',nbmem_tg);
	lk   = lk +1;
	leg{lk} = sprintf('min [%d chocs]',nbmem_tg);
	lk   = lk +1;
	end
	if ~isempty(zeffscl)
	plot(temps,real(zeffscl),'g')
	leg{lk} = 'scaling TS';
	lk   = lk +1;
	end
	if~isempty(zeffrad)
	plot(temps,zeffrad,'c');
	leg{lk} = 'Rad (Matthews)';
	lk   = lk +1;
	plot(temps,zeffrad_cor,'m');
	leg{lk} = 'Rad corrige (Matthews)';
	lk   = lk +1;
	end
	if~isempty(zeff_vloop)
   		plot(temps_zeff_vloop,real(zeff_vloop),'c.');
   		leg{lk} = 'Zeff(Vloop) (Ohm)';
   		lk   = lk +1;
	end
	
	hold off
	legend(leg);
	%title(sprintf('Tore Supra # %d',choc));
	%xlabel('temps (s)')
	ylabel('Zeff tangent')   
	set(gca,'ylim',[1 6])
	grid on
	%zoom xon
	
	
	subplot(4,1,4,'align')	
	leg = {};
	lk        = 1;
	if ~isempty(zeff)
	plot(tzeff,zeff,'r')
	leg{lk} = 'Zebre';
	lk   = lk +1;
	end 
	hold on
	if ~isempty(zeff_tg)
	plot(tzeff_tg,zeff_tg,'b')
	leg{lk} = 'Zebre tangent';
	lk   = lk +1;
	end 
	if~isempty(zeffrad)
	plot(temps,zeffrad,'c');
	leg{lk} = 'Rad (Matthews)';
	lk   = lk +1;
	plot(temps,zeffrad_cor,'m');
	leg{lk} = 'Rad corrige (Matthews)';
	lk   = lk +1;
	end
	if~isempty(zeff_vloop)
   		plot(temps,real(zeff_vloop),'c.');
   		leg{lk} = 'Hugill (Ohm)';
   		lk   = lk +1;
	end
	
	hold off
	legend(leg);
	%title(sprintf('Tore Supra # %d',choc));
	xlabel('temps (s)')
	ylabel('Zeff & Zeff tangent')   
	set(gca,'ylim',[1 6])
	grid on
end
xlabel('temps (s)')
zoom on


function zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz)
% zeffsc = zeffscaling(nbar,ptot,ip,a,R0,itor,gaz);
% nbar : 1019 m-2
% ptot : MW ( equivalent ploss)
% ip   : MA
% a    : petit rayon (m)
% R0   : grand rayon (m)
% itor : courant toroidal (kA)
% gaz  : 2 -> helium, 1 -> Deuterium
%
if gaz == 2
  alpha = [-0.6907    0.1467    0.3994   -2.0644    0.9234   -0.1183];

  c0 = 3.0160;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);
  a4 = alpha(4);
  a5 = alpha(5);
  a6 = alpha(6);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* ip.^a3 .* a.^a4 .* R0.^a5 .* itor.^a6;
end
if gaz == 1

  alpha = [-0.3931    0.1116   3.6174];
  c0 = 185.2499;

  a0 = c0;
  a1 = alpha(1);
  a2 = alpha(2);
  a3 = alpha(3);

  zeffsc = a0 * nbar.^a1 .* ptot.^a2 .* (a./R0).^a3;


end


function out = interp1r(x,y,xx,method)

ind = find(diff(x) <=0) + 1;
if ~isempty(ind)
   x(ind) =[];
   y(ind) =[];
end
out = interp1(x,y,xx,method);



% fit gaussien du profil de puissance rayonnee
% syntaxe : [centre,largeur,cprad] = bolo_width(choc,plotonoff);
%           x = linspace(0,1);   % coordonnees r/a 
%           fit = exp(-(x-centre) .^ 2 ./ largeur .^2);
% cprad = certification : 0 = ok
% choc  = numero du choc
% si plotonoff = 1 -> trace le graphe du resultat
function [centre,largeur,cprad] = zbolo(choc,plotonoff)

% test des entrees
if nargin < 2
	plotonoff = 0;
end

% valeur par defaut
centre  = NaN;
largeur = NaN;

% lecture des puissances 
[prad,tprad,iprad,cprad] = tsbase(choc,'gpbolo1');
if isempty(prad)
	return
end
[cert,vers,date,heure,unix,uniy,uniz,nomdon,type]=tsbase_cert(cprad);
if any(cert(1:3) == -2)
	return
end

% svd de filtrage
[u,s,v] = svd(prad,0);
% on garde le profil principal
forme_h = v(1:8,1)';
forme_l = v(end:-1:9,1)';

% lecture des angles
theta   = tsmat(choc,'DBOLO;CONFIG_CAM1;THETA');
if isempty(theta)
	cprad(2:(1+cprad(1))) = -2;
	return
end
pimp    = tsmat(choc,'DBOLO;CONFIG_CAM1;IMPACT');
ssol    = tsmat(choc,'DBOLO;CONFIG_CAM1;SOLIDE');

lz      = pimp ./ cos( pi/2 - abs(theta));
lz      = max(lz) -lz;
zsa     = lz ./ max(lz);
dforme_h  = gradient(forme_h);
if any(dforme_h > 0)
		indm = min(find((dforme_h>0) & (zsa < max(zsa))));
		forme_h(indm:end) = forme_h(indm);
end 
forme_h   = forme_h - min(forme_h);
forme_h   = forme_h .* abs(cos(theta));
forme_h   = forme_h ./ max(forme_h);
 
dforme_l  = gradient(forme_l);
if any(dforme_l > 0)
		indm = min(find((dforme_l>0) & (zsa < max(zsa))));
		forme_l(indm:end) = forme_l(indm);
end 
forme_l   = forme_l - min(forme_l);
forme_l   = forme_l .* abs(cos(theta));
forme_l   = forme_l ./ max(forme_l);
 

x       = linspace(0,1);
ff      = pchip(zsa,forme_h,x); 

% centre 
indc    = min(find(ff == max(ff)));
centre  = x(indc);
xc      = x(1:indc);
ffc     = ff(1:indc);
d       = abs(ffc - max(ffc)./exp(1));
lc      = xc(min(find(d == min(d))));
largeur = centre -lc;

if largeur > 0.5
   	disp('largeur non valide')
 	ff      = pchip(zsa,forme_l,x); 

	% centre 
	indc    = min(find(ff == max(ff)));
	centre  = x(indc);
	xc      = x(1:indc);
	ffc     = ff(1:indc);
	d       = abs(ffc - max(ffc)./exp(1));
	lc      = xc(min(find(d == min(d))));
	largeur = centre -lc;
	forme   = forme_l;
	if largeur > 0.5
  		disp('largeur non valide')
  		largeur = NaN;
		cprad(2:(1+cprad(1))) = -2;
	end	
	cprad(2:(1+cprad(1)))  = -1;

else
	forme   = forme_h;
end


if plotonoff
	h=figure;
	set(h,'defaultaxesfontsize',12,'defaultaxesfontweight', ...
	'bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
	fit = exp(-(x-centre) .^ 2 ./ largeur .^2);
	fit = max(ff) .* fit ./ max(fit);  
	plot(zsa,forme,'or',x,ff,'.b',x,fit,'k');
	title(int2str(choc))
	xlabel('r/a');
	ylabel('Prad (su)');
	legend('measurement','interpolation','fit','Location','NorthWest');
end

