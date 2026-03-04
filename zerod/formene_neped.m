% calcul des profils de densite en mode H nped donne
% ref pour neped : 
% [nout,fbest] = formene_neped(post.profil0d.xli,post.profil0d.vpr,post.zerod.nem,post.zerod.ane,post.zerod.nebord,post.zerod.negr,neped);
function [nout,fbest] = formene_neped(x,vpr,nem,ane,nebord,negr,neped_in)

persistent ap gap pp_gap

% vecteurs utiles
ve = ones(size(x));
vt = ones(size(vpr,1),1);
vp = trapz(x,vpr,2);
vp_noped = trapz(x(1:(end-1)),vpr(:,1:(end-1)),2);

factor = vt;
factor_mem = vt;

% volume element
dvol = diff(cumtrapz(x,vpr,2),1,2);

for k = 1:100

	% valeur en haut du piedestal
	neped = factor .* min(negr,neped_in);
	% securite
    neped = max(nebord .* 1.001,min((vp .* nem  - (vpr(:,end) .* nebord + vpr(:,end-1) .* neped) ./ 2 .* 0.05) ./ vp_noped,neped));
	% electrons residuels pour le centre
	nres = vp .* nem  - ((vpr(:,end) .* nebord + vpr(:,end-1) .* neped) ./ 2 .* 0.05 + neped .* vp_noped);

	% calcul de la valeur cible galpha
	mxp  = 1 - 0.95 .^ 2;
	warning off
	ga  = nres ./ ((1 + ane) .* nem - neped);
	gar = ga ./ vp; 
	warning on

	% valeur initiale de alpha
	if isempty(ap) 
	    ap  = logspace(-3,2,1001) - 2;
	    gap= (1 - mxp .^ (ap + 1)) ./ (ap + 1);
	    indbad = find(~isfinite(gap));
	    ap(indbad) = [];
	    gap(indbad) = [];
        pp_gap = pchip(gap,ap);
    end
    
    alpha = ppval(pp_gap,gar);
    %figure(21);clf;plot(gar,alpha,'ob',gar,alpha_alt,'.r');drawnow
    
	alpha((gar > max(gap)) | (gar < min(gap))) = 0;
	
	% profile
	nep = nebord * ve;
	nep(:,end - 1) = neped;
	nep(:,1:(end - 2)) = neped * ve(1:(end-2)) + (((1 + ane) .* nem - neped) * ve(1:(end-2))) .* (1 - (vt * x(1:(end-2))) .^ 2) .^ (alpha * ve(1:(end-2)));

	% norme
	factor = factor + (1/sqrt(k)) .* (vp .* nem  - sum((nep(:,1:end-1) + nep(:,2:end)) ./ 2 .* dvol,2)) ./ (vp .* nem);
    %figure(21);clf;plot(1:length(factor),factor,1:length(factor_alt),factor_alt,'.r');drawnow
	if all(abs(factor - factor_mem) < 1e-3)
	    break;
	end
 
end

nout  = nep;
fbest = factor;

return
figure(23);clf
subplot(2,2,1);
plot(ap,gap,alpha,gar,'or');
subplot(2,2,2);
plot(x,nep,'b',0,(1+ane) .* nem ,'or',1,nebord,'or',0.95,neped,'or');
subplot(2,2,3)
plot(factor);hold on;plot((vp .* nem  - trapz(x,nep .* vpr,2)) ./ (vp .* nem),'r');
drawnow

