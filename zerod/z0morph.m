% fonction de morphing des surface de flux
% Rin, Zin : coordonnees de la surface non deformee
% Rsepa,Zsepa : coordonnees de la vrai separatrice
% Rmom,Zmom ; coordonnees de la separtrice donnees par les moments (ou de la dsmf)
% xin : coordonnees radiale de la surface de flux
% expo : exposant de la fonction de morphing
function [Rx,Zx,z0] = z0morph(Rin,Zin,Rsepa,Zsepa,Rmom,Zmom,xin,expo,lkz)

% acceleration de la fonction avec un seul tri
persistent Rsepa_mem Zsepa_mem Rsepa_sort Zsepa_sort 

%gestion de la memorisation
if isempty(Rsepa_mem) && isempty(Zsepa_mem)
    Rsepa_mem = Rsepa;
    Zsepa_mem = Zsepa;
    Rsepa_sort = [];
    Zsepa_sort = [];
    %disp('1')
elseif ~all(size(Rsepa_mem) == size(Rsepa)) || any(Rsepa_mem(:) ~= Rsepa(:)) || any(Zsepa_mem(:) ~= Zsepa(:))
    %size(Rsepa_mem)   
    %size(Rsepa)
    %try
    %    max(Rsepa_mem(:) - Rsepa(:))
    %    max(Zsepa_mem(:) - Zsepa(:))
    %end
    %disp('2')
    Rsepa_mem = Rsepa;
    Zsepa_mem = Zsepa;
    Rsepa_sort = [];
    Zsepa_sort = [];    
%else
%    disp('3')
end

% pas minimal en angle
da = pi ./ 21 ./ size(Rin,2);

% vecteur auxiliaire
vt = ones(size(Rin,1),1);
vu = ones(1,size(Rin,2));
vs = ones(1,size(Rsepa,2));

% calul de l'offeset vertical
Rsepa_max = max(Rsepa,[],2);
%maskrmax  = (Rsepa == Rsepa_max(:,ones(1,size(Rsepa,2))));
%maskrmax  = (Rsepa == (max(Rsepa,[],2) * ones(1,size(Rsepa,2))));
%z0         = sum(Zsepa .* maskrmax,2) ./ sum(maskrmax,2);
%Zsepa      = Zsepa - z0 * vs;
%Zsepa      = Zsepa - z0(:,ones(size(vs)));
Raxe       = 0.5 .* (min(Rin,[],2) + max(Rin,[],2));


% deformation
%cx   = (Rin - Raxe * vu) + sqrt(-1) .* Zin;
cx   = (Rin - Raxe(:,ones(size(vu)))) + sqrt(-1) .* Zin;
thx = angle(cx);
thx(thx<0) = thx(thx<0) + 2 .* pi;
rhox = abs(cx);
[thx,indx] = sort(thx,2);
for k = 1:size(indx,1)
	rhox(k,:)       = rhox(k,indx(k,:));
end

% separtrice a ce temps version analytique
%cl   = (Rmom - Raxe * vu) + sqrt(-1) .* Zmom;
cl   = (Rmom - Raxe(:,ones(size(vu)))) + sqrt(-1) .* Zmom;
thl = angle(cl);
thl(thl<0) = thl(thl<0) + 2 .* pi;
rhol = abs(cl);
[thl,indl] = sort(thl,2);
for k = 1:size(indl,1)
	rhol(k,:)       = rhol(k,indl(k,:));
end
rhol = cat(2,rhol,rhol,rhol);
thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
indnok = find(any(abs(diff(thl,1,2)) <= da,1));
while(~isempty(indnok))
	thl(:,indnok) =[];
	rhol(:,indnok)  = [];
	indnok = find(any(abs(diff(thl,1,2)) <= da,1));
end

% separtrice a ce temps complete
% le probleme de precision vient de ce block de calcul
%cc   = (Rsepa - Raxe * vs) + sqrt(-1) .* Zsepa;
if ~isempty(Rsepa_sort) && ~isempty(Zsepa_sort) 
    cc   = (Rsepa_sort - Raxe(:,ones(size(vs)))) + sqrt(-1) .* Zsepa_sort;
    mem_sort = 0;
else
    cc   = (Rsepa - Raxe(:,ones(size(vs)))) + sqrt(-1) .* Zsepa;
    Rsepa_sort = NaN .* Rsepa;
    Zsepa_sort = NaN .* Zsepa;
    mem_sort = 1;
end
thc = angle(cc);
thc(thc<0) = thc(thc<0) + 2 .* pi;
rhoc = abs(cc);
if mem_sort
    [thc,indc] = sort(thc,2);
    for k = 1:size(indc,1)
        rhoc(k,:)       = rhoc(k,indc(k,:));
        Rsepa_sort(k,:) = Rsepa(k,indc(k,:));
        Zsepa_sort(k,:) = Zsepa(k,indc(k,:));
    end
end
rhoc = cat(2,rhoc,rhoc,rhoc);
thc = cat(2,thc -2.*pi,thc,thc+2.*pi);
indnok = find(any(abs(diff(thc,1,2)) <= da,1));
while(~isempty(indnok))
    thc(:,indnok)   = [];
    rhoc(:,indnok)  = [];
    indnok = find(any(abs(diff(thc,1,2)) <= da,1));
end


% regle de trois, mise a l'echelle
%  vhc      = (1:size(thx,1))' * ones(1,size(thc,2));
%  vhl      = (1:size(thx,1))' * ones(1,size(thl,2));
%  vhx      = (1:size(thx,1))' * ones(1,size(thx,2));
%  rho_sepa = griddata(vhc,thc,rhoc,vhx,thx,'cubic');
%  rho_mom  = griddata(vhl,thl,rhol,vhx,thx,'cubic');
if size(thc,1) == 1
	rho_sepa = spline(thc,rhoc,thx);
	rho_mom  = spline(thl,rhol,thx);
else
	rho_sepa = tsplinet(thc,rhoc,thx);
	rho_mom  = tsplinet(thl,rhol,thx);
end
morf    = (1 - xin .^ expo) + xin .^ expo .* rho_sepa  ./ rho_mom;
rhox = rhox .* morf;
%Rx   = rhox .* cos(thx) + Raxe * vu;
Rx   = rhox .* cos(thx) + Raxe(:,ones(size(vu)));
Zx   = rhox .* sin(thx);

%  save(sprintf('contexte_z0morph_loc_%d',lkz));
%  keyboard

% figure(61);clf;
% plot(thx',morf');
% morf
% drawnow
%  keyboard






