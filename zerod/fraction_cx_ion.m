% faction of neutral that get charge exchange instead of ionisation
% test :
%   fcxion = fraction_cx_ion(post.profil0d.tep(:,end),post.profil0d.tip(:,end),post.profil0d.nep(:,end),post.profil0d.nip(:,end));

function fcxion = fraction_cx_ion(tex,tix,nex,nix)

% calcul des sections efficaces 
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(tex,tix);

%% equilibre entre 1s et 2s pour le neutres de centre 
alphas = nex .* sss ./ ( nex .* sss + Ass);
% etat d'equilibre 1s/2s
sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
%
fcxion = nix .* svcx ./ (nex .* svcx + nex .* sviss + nix .* sii);
