% diagramme LH
function [x,fpout,xlh,dlh,efficiency,rapnegpos,lc,hc,acc,landau] = ...
	 z0lhacc2lobes(flh,npar0,width,agaz,zgaz,temps,x,nep,tep,qp,Raxe,rmx, ...
	 spr,vpr,Bt,plh,xlhin,dlhin,transitoire,directivity,friplh,kx,plotonoff,upshiftmode,npar_neg)
	 
if nargin < 23
	plotonoff = 0;
end
if nargin < 24
	upshiftmode = 'linear';
end
switch upshiftmode
case {'newmodel','newmodel + tail'}
	if nargin < 25
		% d'apres les simulations completes ALOHA/C3PO/LUKE, cet effet n'existe pas 
		npar_neg = - npar0;
	elseif isempty(npar_neg) || (npar_neg> -1)
		% d'apres les simulations completes ALOHA/C3PO/LUKE, cet effet n'existe pas 
		npar_neg = - npar0;
	end
otherwise
	if nargin < 25
		% d'apres les simulations completes ALOHA/C3PO/LUKE, cet effet n'existe pas 
		npar_neg = - npar0 .* pi;
	elseif isempty(npar_neg) || (npar_neg> -1)
		% d'apres les simulations completes ALOHA/C3PO/LUKE, cet effet n'existe pas 
		npar_neg = - npar0 .* pi;
	end
end

directivity = min(1,max(0,abs(directivity)));
ve     = ones(size(x));	 	 
	 
[x,fpoutp,xlh,dlh,lc,hc,acc,landau,efficiency] = ...
	 z0lhacc(flh,npar0,width,agaz,zgaz,temps,x,nep,tep,qp,Raxe,rmx, ...
	 spr,vpr,Bt,plh .* directivity,xlhin,dlhin,transitoire,kx,plotonoff,upshiftmode);	 
% normalisation
fjlhoutp = imag(fpoutp) .* (((plh .* directivity) ./ max(1,trapz(x,vpr .* real(fpoutp),2))) * ve);
fpoutp   = real(fpoutp) .* (((plh .* directivity) ./ max(1,trapz(x,vpr .* real(fpoutp),2))) * ve);
 
[x_,fpoutn,xlh_,dlh_,lc_,hc_,acc_,landau_,rapnegpos] = ...
	 z0lhacc(flh,npar_neg,width,agaz,zgaz,temps,x,nep,tep,qp,Raxe,rmx, ...
	 spr,vpr,Bt,plh .* max(0.01,1 - directivity - friplh),xlhin,dlhin,transitoire,kx,plotonoff,upshiftmode);
% normalisation
fjlhoutn = imag(fpoutn) .* (((plh .* max(0,1 - directivity - friplh)) ./ max(1,trapz(x,vpr .* real(fpoutn),2))) * ve);
fpoutn   = real(fpoutn) .* (((plh .* max(0,1 - directivity - friplh)) ./ max(1,trapz(x,vpr .* real(fpoutn),2))) * ve);

% rapport d'efficacite pour cree le profil de courant
rapnegpos = min(1,rapnegpos ./ max(1,efficiency));

% profil de puissance (< 0 pour contre courant)
% changement si les deux depots sont superposes.
fpout  = fpoutp  + fpoutn + sqrt(-1) .* (fjlhoutp - fjlhoutn);
% efficacite globale tel que Ilh  = efficiency * Plh /R0 / <ne> = int(S'*Jlh,x=0..1)
% l'effet du Zeff est ajoute apres
efficiency = efficiency .* (directivity - rapnegpos .*  max(0,1 - directivity - friplh));
