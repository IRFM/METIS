function [zeff,zmszl,zeffp,n1p,nhep,nzp,nip] = z0zeff_DHe3(x,nep,tep,tip,nwp,zu1,zu2,zimp,zmax,zeff0,mode,nhem,gaz,nem,tim,frhe0,vpr,ftauhe,faccu,rimp,temps,Sn_fraction,iso,nTm,cmin)
% la coherence est sur les densites moyennes

% vecteurs utils
ve = ones(size(x));
vt = ones(size(tip,1),1);

% coherence nep et nem
nepm    = trapz(x,nep .* vpr,2) ./ trapz(x,vpr,2);
%  figure(21)
%  plot(1:length(nepm),nepm,'or',1:length(nepm),nem,'b.');
%  drawnow
factor_ave2prof = (nepm ./ nem) * ve;

% choix de l'impurete qui doit est traitee en transport neoclassique
zneo = fix(zimp);


% Pas d'accumulation de cendre dans ce cas
% % accumulation de l'helium (le zeff de reference est donne sans he du au cendre)
% zeff0   = (zeff0 .* nem + 4 .* max(0,nhem - frhe0 .* nem)) ./ nem;

% securite sur le zeff
zeffmin = max(1.01 .* vt,1 + 2 .* (nhem  + frhe0 .* nem)./ nem);
nzmax   = max(0,nem - 2 .* nhem  - 1e13) ./ zu1;
zeffmax = min(zmax - 0.1,(4 .* nhem + 1e13 + zu2 .* nzmax) ./ nem);
if any((zeffmin > zeffmax) & (nhem > 0) & (nzmax > 0))
  error('Z0ZEFF: to much helium in the plasma');
end
zeff0  = min(zeffmax,max(zeff0,zeffmin));

%  figure(21)
%  clf
%  plot(1:length(nhem),nhem .* trapz(x,vpr,2),'ob',1:length(nhem),trapz(x,nhep .* vpr,2),'r',1:length(nhem),(nhem .* trapz(x,vpr,2) - trapz(x,nhep .* vpr,2)),'g');
%  drawnow

% solve for densities
nDm =  ((nem .* (1 - 2 .* frhe0) - nTm) .* zu2 - (nem .* (zeff0 - 4 .* frhe0) - nTm) .* zu1) ./ ...
       ((1 + cmin + 2 .* iso) .* zu2 - (1 + cmin + 4 .* iso) .* zu1);
   
nhem =   iso .* nDm;
% securite sur nhem
%nhem = min(nhem,nem./2.1);
% He3 suppose homotetique a ne 
% la densite de bord est taue/tauhe fois plus faible.
% l'helium 3 est suppose entierement ionise (contrairement au impuretes)
%nhep   = ((nhem ./nem) * ve) .* nep;
nhea    = (max(0,nhem) ./ nem) .* nep(:,end) ./ ftauhe;
nhea    = min(nep(:,end) ./ 2.1, nhea);
nhep0   = max(0,nep - nep(:,end) * ve);
nhep0(find(any(nhep0(:,1:end-1) == 0,2)),:) = 1e13; 
nhem0   = trapz(x,nhep0 .* vpr,2) ./ trapz(x,vpr,2);
nhep    = nhep0 .* ((max(0,nhem - nhea) ./ max(1,nhem0)) * ve) + nhea * ve; 
nhep    = factor_ave2prof .* nhep .* ((nhem ./ max(eps,trapz(x,nhep .* vpr,2) ./ trapz(x,vpr,2))) * ve);

% He4 simply factor on nep
nhe4p  = frhe0 .* nep;
%
nHm     = cmin .* nDm;
n1m     = nDm + nHm + nTm;
nimpm = (nem .* (1 - 2 .* frhe0) - nDm - nHm - nTm - 2 .* nhem) ./ zu1; % W is encode in zu1
nmaxm = nimpm .* rimp;
% densite ionique
% nim = n1m + nhem + nimpm .* (1 + rimp) + nwm + nem .* frhe0;


% % calcul de la valeur moyenne des impuretes (seule valeur definie a partir de zeff moyen
% % nimpm =  ((zeff0 - 1) .* nem - 2 .* nhem) ./ (zu2 - zu1);
% % nmaxm = nimpm .* rimp;
% % n1m   = nem - 2 .* nhem - zu1 .* nimpm;

	
if length(zu1)> 1
	% cas avec tungsten trait? ? part
	mode = 0;
end
nerp = max(1,nep - 2 .* nhep);

% calcul coherent pour le tungstene
if length(zu1) > 1
    nwpm   = trapz(x,nwp .* vpr,2) ./ trapz(x,vpr,2);
    % electron manquant
    %nwez1m   = max(0,nem -n1m - 2 .* nhem -  (zimp + rimp * zmax) .* nimpm);
    nwz1m  = (zu1 - (zimp + rimp * zmax)) .* nimpm;
    nwz1m  = nwz1m .* (nwpm > 0);
    % les 2 calculs sont equivalents
    %figure(21);clf;plot(1:length(n1m),nwez1m,'or',1:length(n1m),nwz1m,'.b');drawnow
    nwz2m  = (zu2 - (zimp .^ 2 + rimp * zmax .^ 2)) .* nimpm ;
    nwz2m  = nwz2m.* (nwpm > 0);
    nwz1   = z0wavez(tep) .* nwp;
    if Sn_fraction > 0
        nwz1   =  (1 - Sn_fraction) .* nwz1 + Sn_fraction .* z0snavez(tep) .* nwp;
    end
    nwz1_nn = nwz1;
    nwz1    = min(nerp,factor_ave2prof .*  nwz1 .* ((nwz1m ./  max(1,trapz(x,nwz1 .* vpr,2) ./ trapz(x,vpr,2))) * ve));
    %figure(21);clf;plot(x,nerp - nip ,'r',x,nwz1,'b');drawnow;
    if Sn_fraction > 0
        factor_nn = nwz1 ./ max(1,nwz1_nn);
        nwz2   = ((1 - Sn_fraction) .* z0wavez(tep).^ 2 + Sn_fraction .* z0snavez(tep) .^ 2) .* nwp .* factor_nn;
    else
        nwz2   = z0wavez(tep) .* nwz1;
    end
    nwz2   = factor_ave2prof .* nwz2 .* ((nwz2m ./  max(1,trapz(x,nwz2 .* vpr,2) ./ trapz(x,vpr,2))) * ve);
else
    nwz1 = zeros(size(nep));
    nwz2 = zeros(size(nep));   
end
% sans le tungstene
nerp    = max(1,nerp - nwz1);




% selon le mode
switch mode
case {0,5,6,7}
	% pas d'effet de profil de zeff
	nzp = nerp;
	% normalisation de nzp pour retrouver nimpm
	nzpm    = trapz(x,nzp .* vpr,2) ./ trapz(x,vpr,2);
	nzp     = nzp .* ((nimpm ./ nzpm) * ve);
	
otherwise
	% profil d'impurete lourde 
	if length(zu1)> 1
		% pas d'effet de profil de zeff en dehors de W
	    nzp = nerp;
		% normalisation de nzp pour retrouver nimpm
		nzpm    = trapz(x,nzp .* vpr,2) ./ trapz(x,vpr,2);
		nzp     = nzp .* ((nimpm ./ nzpm) * ve);
	else
		zacc = max(zmax,zimp);

		if faccu >= 0
			nzp_max  = (nerp ./ (nerp(:,1) * ve)) .^ zacc .* ((tip(:,1) * ve) ./ max(30,tip)) .^ (zacc /2 -1);
		else
			nzp_max  = (nerp ./ (nerp(:,1) * ve)) .^ zacc;
		end
	    if any(~isfinite(nzp_max(:)))
            nzp_max(~isfinite(nzp_max(:))) = nerp(~isfinite(nzp_max(:)));
        end
	
		nzp_max  = nzp_max ./ max(1,(max(nzp_max,[],2) * ve));
		nzp_imp  = nerp ./ max(1,(max(nerp,[],2) * ve));
	
		%nzp_imp  = (nep ./ (nep(:,1) * ve)) .^ zimp .* ((tip(:,1) * ve) ./ max(30,tip)) .^ (zimp /2 -1);
		%nzp_imp  = nzp_imp ./ max(1,(max(nzp_imp,[],2) * ve));
	
		nzp_no  = nerp ./ max(1,(max(nerp,[],2) * ve));
	
		% calcul des valeurs moyennes
		nzmaxm  = trapz(x,nzp_max .* vpr,2) ./ trapz(x,vpr,2);
		nznom   = trapz(x,nzp_no .* vpr,2) ./ trapz(x,vpr,2);
		nzimpm  = trapz(x,nzp_imp .* vpr,2) ./ trapz(x,vpr,2);
                %figure(21);clf;plot(temps,nzmaxm,'r-',temps,nznom,'b.-',temps,nzimpm,'c:');drawnow
		% normalisation
		nzp_max = nzp_max .* ((nmaxm ./ max(eps,nzmaxm)) * ve); 
		nzp_no  = nzp_no .* ((nmaxm ./  max(eps,nznom)) * ve); 
		nzp_imp = nzp_imp .* ((nimpm ./  max(eps,nzimpm)) * ve); 
		% les rapports sont fixes en valeur moyenne
		nzp_max = abs(faccu) .* nzp_max + (1- abs(faccu)) .* nzp_no;
		nzp     = nzp_max  + nzp_imp;
		% normalisation de nzp pour retrouver nimpm
		nzpm    = trapz(x,nzp .* vpr,2) ./ trapz(x,vpr,2);
                %figure(22);clf;plot(temps,nzpm);drawnow
		nzp     = nzp .* ((nimpm ./ max(1,nzpm)) * ve);
	end
end

% Modif D. Moreau 30/04/2014 pour donn?es externes sur ZEFF.
% completer le 27/08/2014 par J-F Artaud

if isappdata(0,'ZEFF_SHAPE_EXP')
    zeffexp = getappdata(0,'ZEFF_SHAPE_EXP');
    zeffp  = max(1,interp1_ex(zeffexp.temps,zeffexp.zeff,temps,'nearest','extrap'));
    zeffp  = pchip(zeffexp.x,zeffp,x);
    indnok = find(any(~isfinite(zeffp),2));
    indok  = find(all(isfinite(zeffp),2));
    zeffp(indnok,:) = ones(length(indnok),1) * mean(zeffp(indok,:),1);
    % normalisation au zeff moyen lineique de consigne
    zeffp    = zeffp .* ((zeff0./trapz(x,zeffp,2)) * ve);
    % calcul de nzp
    nzp = (nep .* (zeffp -1) - 2 .* (nhep + nhe4p)  + nwz1 - nwz2) ./  ((zimp .^ 2 + rimp * zmax .^ 2) - (zimp + rimp * zmax));   
    % normalisation de nzp pour retrouver nimpm (et le zeffm)
    nzpm    = trapz(x,nzp .* vpr,2) ./ trapz(x,vpr,2);
    nzp     = nzp .* ((nimpm ./ max(1,nzpm)) * ve);
end

nzp  = factor_ave2prof .* nzp;
nerp = max(1,nerp - (zimp + rimp * zmax) .* nzp);
%  boron profile is assumed to be proportional to  hydrogen one.

% normalisation de nip pour retrouver n1m
nipm    = trapz(x,nerp .* vpr,2) ./ trapz(x,vpr,2);
n1p     = nerp .* ((n1m ./ nipm) * ve);
n1p     = factor_ave2prof .* n1p;

% compute nion
nip = n1p + nzp .* (1 + rimp) + nwp + nhep + nhe4p;

%  figure(21)
%  plot(1:length(nwpm),nwz1m,'or',1:length(nwpm),newm,'b.');
%  drawnow

% recalcul final  pour augmenter la precision 
zeffp  = real((n1p + 4 .* (nhep + nhe4p) + (zimp .^ 2 + rimp * zmax .^ 2) .* nzp + nwz2) ./ nep); 
zeffl  = min(zeffmax,max(zeffmin,real(trapz(x,zeffp,2))));
% zeff volumique moyen
zeff = min(zeffmax,max(zeffmin,real(trapz(x,nep .* zeffp.* vpr,2) ./ trapz(x,nep.* vpr,2))));
% rapport zeff moyen volumique sur zeff moyen lineique
zmszl  = max(0.1,min(10,zeff ./ zeffl));




% no plot
return

if length(zu1)> 1	
        zu1  = zu1 * ve;
        zu2  = zu2 * ve;
end


% figure(101);clf;
% subplot(2,2,1);
% plot(x,(nip + (zimp + rimp * zmax) .* nzp + 2 .* nhep + nwz1 - nep)./nep);
% subplot(2,2,2);
% plot(x,(nip + (zimp .^ 2 + rimp * zmax .^ 2) .* nzp + 4 .* nhep + nwz2)./nep - zeffp);
% subplot(2,2,3);
% plot((nip + zu1 .* nzp + 2 .* nhep - nep)./nep);
% subplot(2,2,4);
% plot((nip + zu2 .* nzp + 4 .* nhep)./nep - zeffp);
% drawnow
% 
% test = abs((nip + zu1 .* nzp + 2 .* nhep - nep)./nep);
% if any(test(:) > 1) || any(n1m < 0) || any (n1m > nem) || any(nzp(:) < 0)
% 	disp('in z0zeff');
% 	keyboard
% end
% 
