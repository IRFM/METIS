% ZPELLETCOOR cette fonction calcule sur la droite d'ablation les grandeurs plasmas
%---------------------------------------------------------------------------------------
% fichier zpelletcoor.m ->  zpelletcoor
%
%
% fonction Matlab 7 :
%
% Cette fonction calcule pour l'injection de glacon, sur la droite d'ablation les grandeurs plasmas 
%  
% syntaxe  :
%  
%      [inout,R,Z,ne,te]=zpelletcoor(origine,pos,x,equi,prof,source,impur,neo,{plotonoff});
%
% entrees :
%
%     origine        = structure decrivant la droite d'ablation
%		       origine contient soit r0,z0,r1,z1 
%					soit r0,z0,theta
%			theta tel que : r1 = r0 + u * cos(theta) et z1 = z0 + u * sin(theta)
%
%     pos            = position le long de la droite d'ablation (vecteur)
%                      pos est compte 0 a partir de la dsmf et augment dans le plasma
%                      pos est en m
%                      si pos est complexe, pos est compris comme des coordonnees (R,Z)
%
%     x              =  param.gene.x 
%     equi           =  datak.equi
%     prof           =  datak.prof
%     source         =  datak.source
%     impur          =  datak.impur
%     neo            =  datak.neo
%     plotonoff      =  si = 1 trace le graphe 2D d'intersction (optionnel)
%
% sorties :
% 
%     inout          = 1 si le point correspondant est dans plasma (taille de pos)
%     R,Z            = coordonnees (R,Z) des point correspondant a pos  
%     ne,te,...      = valeur des champs en (R,Z)
%     
% fonction ecrite par J-F Artaud
% version 4.0, du 20/06/2007.
% 
% 
% liste des modifications :  cf cvs.
% 
%---------------------------------------------------------------------------------
%
function [inout,rho,R,Z,volume,ne,te,nions,tions,gradbob,b,indice]=zpelletcoor(origine,pos,x,equi,prof,source,impur,neo,plotonoff,methode)

% selon les entrees
if nargin < 9
	plotonoff = 0;
end
if nargin < 10
	methode = 'cubic';
end

% si pas de point dans le plasma
vpos  = ones(size(pos));
inout = 0   .* vpos;
R     = NaN .* vpos;
Z     = NaN .* vpos;
ne    = NaN .* vpos;
te    = NaN .* vpos;

% selon le contenu de pos
if ~any(imag(pos))

	% description de origine 
	% (r0,z0) et (r1,z1) ou theta  (en radian)
	% theta est tel que r1 = r0 + L * cos(theta) et z1 = z0 + L * sin(theta)
	%  datak est une structure a un temps de cronos
	
	% pos est la position en m depuis le point (rdsmf,zdsmf) sur la droite d'injection
	if isfield(origine,'theta')
		cs = cos(origine.theta);
		ss = sin(origine.theta);	
	else
		L = sqrt((origine.r1 - origine.r0) .^ 2  + (origine.z1 - origine.z0).^ 2);
		cs = (origine.r1 - origine.r0) ./ L;
		ss = (origine.z1 - origine.z0) ./ L;	
	end
	% calcul des increments sur la DSMF
	Rm = double(squeeze(equi.R(1,end,:)));
	Zm = double(squeeze(equi.Z(1,end,:)));
	
	if plotonoff
		h = findobj(0,'type','figure','tag','pelletcoor');
		if isempty(h)
			h=figure('tag','pelletcoor');
		else
			figure(h);
		end   
		clf
		set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
		'defaultlinelinewidth',1,'color',[1 1 1])
	
		lref = sqrt((min(Rm) - max(Rm)).^2 + (min(Zm) - max(Zm)).^2);
		plot(Rm,Zm,'b',[origine.r0,origine.r0 + lref .* cs],[origine.z0,origine.z0 + lref .* ss],'r', ...
		origine.r0,origine.z0,'*r');
	end
	Rp = cat(1,Rm(2:end),Rm(1));
	Zp = cat(1,Zm(2:end),Zm(1));
	% intersection avec la dsmf
	a  = (Rm - Rp) .* ss - (Zm -Zp) .* cs;
	b  = (Rm - origine.r0) .* ss  - (Zm - origine.z0) .* cs;
	a(a == 0 ) = NaN;
	xl  = b ./ a ;
	indok = find((xl >= 0) & (xl <= 1));
	if isempty(indok)
		fprintf('No point in plasma \n');
		return
	end
	if plotonoff
		hold on
		plot(Rm(indok),Zm(indok),'oc',Rp(indok),Zp(indok),'oc');
	end
	% si non tangent
	if length(indok) > 1 
		dd = sqrt((2 .* origine.r0 - Rm(indok) - Rp(indok)) .^ 2 + (2 .* origine.z0 - Zm(indok) - Zp(indok)) .^ 2); 
		indfirst = find(dd == min(dd),1);
		indok    = indok(indfirst);
	end
	if plotonoff
		plot(Rm(indok),Zm(indok),'*m',Rp(indok),Zp(indok),'*m');
	end
	Rdsmf = Rm(indok) + (Rp(indok) - Rm(indok)) .* xl(indok);
	Zdsmf = Zm(indok) + (Zp(indok) - Zm(indok)) .* xl(indok);
	if plotonoff
		plot(Rdsmf,Zdsmf,'or');
	end
	% calcul des points R,Z sur la coorde :
	R = Rdsmf + pos .* cs;
	Z = Zdsmf + pos .* ss;
	if plotonoff
		plot(R,Z,'*b');
	end
else
	R = real(pos);
	Z = imag(pos);
end	
% verification que le point est dans le plasma
inout = zinout(Rm,Zm,R,Z);

% mapping des profils sur la grille de l'equilibre
Req   = double(squeeze(equi.R(1,2:end,1:end-1)));
Zeq   = double(squeeze(equi.Z(1,2:end,1:end-1)));
rhoeq = double(squeeze(equi.rhoRZ(2:end)));

% donnees physique en sortie
vpr_eq = pchip(x .* equi.rhomax,equi.vpr,rhoeq);
voleq  = cumtrapz(rhoeq,vpr_eq,2);
rho    = griddata(Req,Zeq,rhoeq' * ones(1,size(Req,2)) ,R,Z,methode);
ne     = pchip(x .* equi.rhomax,prof.ne,rho);
te     = pchip(x .* equi.rhomax,prof.te,rho);
nions  = NaN .* ones(length(rho),size(impur.impur,3));
tions  = NaN .* ones(length(rho),size(impur.impur,3));
for k =1:size(nions,2)
	nions(:,k) =  pchip(x .* equi.rhomax,squeeze(impur.impur(1,:,k)),rho);
	tions(:,k) =  pchip(x .* equi.rhomax,prof.ti,rho);
end
indice = pchip(x .* equi.rhomax,1:length(x),rho);
volume = pchip(rhoeq,voleq,rho);

% calcule de gradbob
switch methode 
case 'nearest'
	gradbob = 1 ./ R;
	b       = pchip(x .* equi.rhomax,equi.F,rho) ./ R;
otherwise
	BReq    = double(squeeze(equi.BR(1,2:end,1:end-1)));
	BZeq    = double(squeeze(equi.BZ(1,2:end,1:end-1)));
	BPHIeq  = double(squeeze(equi.BPHI(1,2:end,1:end-1)));
	B       = sqrt(BPHIeq .^ 2 + BZeq .^ 2 + BReq .^2);
	da      = equi.rhomax ./ length(prof.ne) ./ 2;
	nbp     = length(R);
	Rx      = R(:)';
	Zx      = Z(:)';	
	Rv      = cat(1,Rx,Rx - da,Rx + da,Rx,Rx);
	Zv      = cat(1,Zx,Zx,Zx,Zx - da,Zx + da);
	bdb     = griddata(Req,Zeq,B ,Rv,Zv ,methode);
	
	bdrp    = bdb(3,:);
	bdrm    = bdb(2,:);
	bdzp    = bdb(5,:);
	bdzm    = bdb(4,:);
	b       = bdb(1,:);
	dBdR    = (bdrp - bdrm) ./ 2 ./ da; 
	dBdZ    = (bdzp - bdzm) ./ 2 ./ da; 
	gradbob = sqrt(dBdR .^ 2 + dBdZ .^ 2) ./ b; 
	if size(R,1) > 1
		gradbob = gradbob';
	end
end




