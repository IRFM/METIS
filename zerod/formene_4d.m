% calcul des profils de densite en mode H avec 4 degrees de liberte
% ref pour neped : 
% [nout,fbest] = formene_4d(post.profil0d.xli,post.profil0d.vpr,post.zerod.nem,post.zerod.ane,post.zerod.nebord,post.zerod.negr);
function [nout,fbest] = formene_4d(x,vpr,nem,ane,nebord,negr,expo)

persistent ap gap pp_gap

% vecteurs utiles
ve = ones(size(x));
vt = ones(size(vpr,1),1);
vp = trapz(x,vpr,2);
vp_noped = trapz(x(1:(end-1)),vpr(:,1:(end-1)),2);

factor = vt;
factor_mem = vt;

if nargin < 7
  expo = -0.7;
elseif isempty(expo)
  expo = -0.7;
elseif expo == 0
  expo = -0.7;
end
% volume element
dvol = diff(cumtrapz(x,vpr,2),1,2);

for k = 1:100

	% valeur en haut du piedestal
	% Scaling make from data in refrence:
        % Multi-machine comparisons of H-mode separatrix densities and edge profile behaviour in the ITPA SOL and Divertor Physics Topical Group
        % A. Kallenbach a, N. Asakura b, A. Kirk c, A. Korotkov, M.A. Mahdavi , D. Mossessian , G.D. Porter
        % Journal of Nuclear Materials 337?339 (2005) 381?385
        % the best fit expo = -0.7
	neped = factor .* min(negr,nebord .* (nebord ./ negr) .^ expo);
	% securite
	%neped_old = max(nebord .* 1.001,min(vp .* nem  - ((vpr(:,end) .* nebord + vpr(:,end-1) .* neped) ./ 2 .* 0.05) ./ vp_noped,neped));	
        neped = max(nebord .* 1.001,min((vp .* nem  - (vpr(:,end) .* nebord + vpr(:,end-1) .* neped) ./ 2 .* 0.05) ./ vp_noped,neped));
        %if any(any(neped_old - neped))
	%    figure(21);clf;plot(1:length(neped),neped_old,'b', 1:length(neped),neped,'r');drawnow
        %end
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
        %alpha  = interp1(gap,ap,gar,'pchip',0);
	
	%alpha_alt = pchip(gap,ap,gar);
    alpha = ppval(pp_gap,gar);
    %figure(21);clf;plot(gar,alpha,'ob',gar,alpha_alt,'.r');drawnow
    
	alpha((gar > max(gap)) | (gar < min(gap))) = 0;
	
%  	alpha_test  = interp1(gap,ap,gar,'pchip',0);
%          if any(alpha_test(:) ~= alpha(:))
%            indnz = find(alpha_test(:) ~= alpha(:))
%            if any(abs(alpha_test(indnz) - alpha(indnz)) > sqrt(eps))
%  	    keyboard
%            end
%            alpha_test(indnz) - alpha(indnz)
%          end
        
	% profile
	nep = nebord * ve;
	nep(:,end - 1) = neped;
	nep(:,1:(end - 2)) = neped * ve(1:(end-2)) + (((1 + ane) .* nem - neped) * ve(1:(end-2))) .* (1 - (vt * x(1:(end-2))) .^ 2) .^ (alpha * ve(1:(end-2)));

	% norme
	%%%factor_alt = factor + (1/sqrt(k)) .* (vp .* nem  - trapz(x,nep .* vpr,2)) ./ (vp .* nem);
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

