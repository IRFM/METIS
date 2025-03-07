function cr = zupdatecompo(inter)
cr = 0;

fun = evalin('base','param.fonction.impur');
if ~strcmp(fun,'zinebcompo')
   disp('Info : Le module d''impuretes selectionner n''est pas compatible avece ce module');
   cr = -1;
   return
end
 
 
if nargin > 0 
 
      ButtonName = questdlg(' Do you want to update the plasma composition and the ionic pressure and density?', ...
                         'Confirmation :', ...
                         'Yes','No','No');

      if strcmp(ButtonName,'No')
         return
      end
end

impur = evalin('base','data.impur'); 
prof  = evalin('base','data.prof'); 
       
cons  = evalin('base','param.cons.impur');
zeffm = evalin('base','data.cons.zeffm');
compo = evalin('base','param.compo');     
phys  = evalin('base','param.phys');     
x     =  evalin('base','param.gene.x');  
ve    = ones(size(x)); 
% calcul de zmax
zmax   = compo.z(end - 1) .* (1 - cons.rimp) + compo.z(end) .* cons.rimp;
zmax   = min(max(compo.z) .* 0.9,zmax);  
zmax   = max(min(compo.z) .* 2.1 , zmax);  

          

% calcul du zeff	                         
if cons.zeff ==0
	% le profil de zeff est une donnee
	impur.zeff = prof.zeff;
else	
	% si cons.zeff == 2 => zeffm est calcule avec la loi d'echelle TS
	% utiliser pour l'injection de glacon
	nbar = trapz(x,prof.ne,2)./1e19;
	if compo.z(1) == 1
		  zeffmr = 4.35 .* nbar .^ -0.4; 	
	elseif compo.z(1) ==2
		  zeffmr = 10 .* nbar .^ -0.7;      
	end
	if (cons.zeff == 2) 
      zeffm = zeffmr;
	end
   ind  = find(~isfinite(zeffm));
   if ~isempty(ind)
      zeffm(ind) = zeffmr(ind);
   end
	zeffm = min(max(zeffm,compo.z(1) .* 1.1),zmax - 0.2);    

	if cons.exposant == 0
%		impur.zeff = zeffm .* ones(size(prof.ne));
		impur.zeff = zeffm * ones(size(x));
	else	
		zeff = prof.ne .^ cons.exposant;
		zeff(end)  = 2 .* zeff(end-1) - zeff(end-2);
		g = trapz(x,zeff,2);
		alpha = (zeffm - compo.z(1)) ./ (g - 1e21 .^ cons.exposant);
		gamma = (g .* compo.z(1)  - zeffm .* 1e21 .^ cons.exposant) ./(g - 1e21 .^ cons.exposant); 
		impur.zeff = (alpha * ve).* zeff + (gamma * ve);
	end
end

% securite zeff
ind = find(impur.zeff < min(compo.z));
if ~isempty(ind)
	impur.zeff(ind) = min(compo.z);
end
ind = find(impur.zeff > zmax);
if ~isempty(ind)
	impur.zeff(ind) = zmax;
end


% calcul de la composition du plasma
[ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(prof.ne,impur.zeff,cons.cmin1,cons.cmin2,cons.rimp, ...
                                         compo.z(1),compo.z(2),compo.z(3),compo.z(4),compo.z(5));

% securite (le plasma ne peut pas contenir que des impuretes)
if any(nion <= 0)
     impur.zeff(impur.zeff > (zmax/2)) = zmax / 2;  
     [ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(prof.ne,impur.zeff,cons.cmin1,cons.cmin2,cons.rimp, ...
                                         compo.z(1),compo.z(2),compo.z(3),compo.z(4),compo.z(5));
end
% mise des resultats dans la structure impur
impur.ae    = ae .* (ae > 0.01) + 0.01 .* (ae <= 0.01);
nion        = nion  .* (nion > 0);
nmin1       = nmin1 .* (nmin1 > 0);
nmin2       = nmin2 .* (nmin2 > 0);
nimp2       = nimp2 .* (nimp2 > 0);
nimp1       = nimp1 .* (nimp1 > 0);
impur.impur = cat(3,nion,nmin1,nmin2,nimp1,nimp2);                                        


% assignation
zassignin('base','data.impur',impur)


% mise a jour de la pression ionique en conservant ti et ni et ae
pion  = prof.ti .* prof.ne .* impur.ae .* phys.e;
ni    = prof.ne .* impur.ae;
zassignin('base','data.prof.ae',impur.ae)
zassignin('base','data.prof.ni',ni)
zassignin('base','data.prof.pion',pion)
 
% calcul de la compositions du plasma 
% en entree : 
%    ne      = densite electronique
%    zeff    = profil de zeff
%    c1      = rapport de la densite du 1er minoriatire sur la densite de l'espece principale
%    c2      = rapport de la densite du 2ieme minoriatire sur la densite de l'espece principale
%    rimp    = rapport de la densite de la 2ieme impurete sur la densite de la premiere impurete
%    zion    = numero atomic (ou charge moyenne) de l'espece principale
%    zmin1   = numero atomic (ou charge moyenne) du 1er minoritaire
%    zmin2   = numero atomic (ou charge moyenne) du 2ieme minoritaire
%    zimp1   = numero atomic (ou charge moyenne) de la 1ere impurete
%    zimp2   = numero atomic (ou charge moyenne) de la 2ieme impurete
%    
% en sortie : ae,nion,nmin1,nmin2,nimp1,nimp2 
% 
% equations :
% ae*ne =nion+nmin1+nmin2+nimp1+nimp2
% ni =ae*ne
% nmin1 = c1 *nion
% nmin2 = c2*nion
% nimp2 = rimp*nimp1
% ne = nion*zion +nmin1 *zmin1 +nmin2 *zmin2 +nimp1*zimp1 +nimp2*zimp2
% ne *zeff = nion*zion^2 +nmin1 *zmin1^2 +nmin2 *zmin2^2 +nimp1*zimp1^2 +nimp2*zimp2^2
% 
function [ae,nion,nmin1,nmin2,nimp1,nimp2]=zcompo(ne,zeff,c1,c2,rimp,zion,zmin1,zmin2,zimp1,zimp2)

% variables de calcul
de =  - zimp1 .^ 2 .* zion - zimp1 .^ 2 .* c1 .* zmin1 - zimp1 .^ 2 .* c2 .* zmin2 - ...
        rimp .* zimp2 .^ 2 .* zion - rimp .* zimp2 .^ 2 .* c1 .* zmin1 -  ...
        rimp .* zimp2 .^ 2 .* c2 .* zmin2 + zion .^ 2 .* zimp1 +  ...
        zion .^ 2 .* rimp .* zimp2 + c1 .* zmin1 .^ 2 .* zimp1 + ...
        c1 .* zmin1 .^ 2 .* rimp .* zimp2 + c2 .* zmin2 .^ 2 .* zimp1 + ...
        c2 .* zmin2 .^ 2 .* rimp .* zimp2;
       
% especes principales            
nion = -ne .* (zimp1 .^ 2 + rimp .* zimp2 .^ 2 - zeff .* zimp1 - zeff .* rimp .* zimp2) ./ de;
nmin1 = c1.*nion;
nmin2 = c2.*nion;

% impuretees:
nimp1 = (zion .^ 2 + c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - ...
         c1 .* zmin1 .* zeff - c2 .* zmin2 .* zeff ) .* ne ./ de;
nimp2 = rimp.*nimp1;

% rapport somme(ni)/ne :
ae = (- zimp1 .^ 2 - rimp .* zimp2 .^ 2 + zeff .* zimp1 + zeff .* rimp .* zimp2 - ...
        zimp1 .^ 2 .* c1 - c1 .* rimp .* zimp2 .^ 2 + c1 .* zeff .* zimp1 + ...
        c1 .* zeff .* rimp .* zimp2 - c2 .* zimp1 .^ 2 - c2 .* rimp .* zimp2 .^ 2 + ...
        c2 .* zeff .* zimp1 + c2 .* zeff .* rimp .* zimp2 + zion .^ 2 + ...
        c1 .* zmin1 .^ 2 + c2 .* zmin2 .^ 2 - zion .* zeff - c1 .* zmin1 .* zeff - ...
        c2 .* zmin2 .* zeff + zion .^ 2 .* rimp + c1 .* zmin1 .^ 2 .* rimp + ...
        c2 .* zmin2 .^ 2 .* rimp - zion .* zeff .* rimp - c1 .* zmin1 .* zeff .* rimp - ...
        c2 .* zmin2 .* zeff .* rimp) ./ de;

