%  Z0SEPARATRIX  courte description  
%------------------------------------------------------------------------------- 
% fichier :  z0separatrix.m  ->  z0separatrix 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   z0dinput = z0separatrix(z0dinput,option) 
%  
% entrees :  
%  z0dinput = 
%  option   = 
%  
% sorties :  
%   z0dinput  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.0  du  08/03/2007  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function z0dinput = z0separatrix_jt60sa_scenario2_create(z0dinput)

% script de recopie de la separatrice dans metis
[R,Z] = separatrice_jt60sa_scenario2_create(z0dinput.cons.ip ./ 1e6,z0dinput.cons.temps);
sepa.R = R;
sepa.Z = Z;

% calcul des moments
% centre pour angle d'integration
rc = mean(sepa.R,2);
zc = mean(sepa.Z,2);
vc = ones(1,size(sepa.R,2));
uc = unwrap(angle((sepa.R-rc*vc) + sqrt(-1) .* (sepa.Z  -zc*vc)));
uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
uc(:,1)   = uc(:,end) + 2 .* pi;
xu    = linspace(0,1,length(vc));
%dudx  = pdederive(xu,uc,2,2,2,1);
%dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
%dudx(:,end) = dudx(:,1);
dRdx  = pdederive(xu,sepa.R,2,2,2,1);
dZdx  = pdederive(xu,sepa.Z,2,2,2,1);
% calcul de R0 et Z0
maskrmax  = (sepa.R == (max(sepa.R,[],2) * vc));
% recalcul des parametres sur le vecteur final
rmin  = min(sepa.R,[],2);
rmax  = max(sepa.R,[],2);
geo.a = 0.5 .* (rmax - rmin);
geo.R = 0.5 .* (rmax + rmin);
zmin  = min(sepa.Z,[],2);
zmax  = max(sepa.Z,[],2);
geo.z0   = (zmax + zmin + sum(sepa.Z .* maskrmax,2) ./ sum(maskrmax,2)) ./ 3;
geo.K    = sgolayfilt((abs(trapz(xu,sepa.Z .*  dRdx,2) ./ pi ./ geo.a) + (zmax - zmin)) ./ 3 ./ geo.a,1,3);

% conservation de la raideur magnetique
sepa.b0 = z0dinput.geo.b0 .* z0dinput.geo.R ./ geo.R; 

rzmax = geo.R;
rzmin = geo.R;
for k = 1:size(sepa.Z,1)
	rzmax(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmax(k))));
	rzmin(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmin(k))));
end
uu   =  angle(rzmax - geo.R + sqrt(-1) .* (zmax - geo.z0));
ul   =  angle(rzmin - geo.R + sqrt(-1) .* (zmin - geo.z0));
tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
tm   =  (tl + tu) ./ 2;
d    =   abs(rzmax + rzmin -  2 .* geo.R) ./ 2 ./ geo.a;
geo.d =  sgolayfilt(0.6 .* d + 0.4  .* sin(tm),1,3);

z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
z0dinput.geo.z0      = geo.z0;     % centre geometrique du plasma en Z (m)
z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
z0dinput.geo.K       = geo.K;     % elongation (b/a)
z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
z0dinput.geo.b0      = sepa.b0;  % champ toroidal a R0

z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)

t  = asin(max(0,min(1,z0dinput.geo.d)));
u  = linspace(0,2.*pi,201);
vu = ones(size(u));
Rtest  = z0dinput.geo.R *vu + (z0dinput.geo.a * vu) .* cos(ones(size(geo.R,1),1) * u + t * sin(u));
Ztest = (z0dinput.geo.a .* z0dinput.geo.K) * sin(u);

h = findobj(0,'type','figure','tag','z0geosepa');
if isempty(h)
h=figure('tag','z0geosepa');
else
figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
zplotprof(gca,z0dinput.cons.temps,Rtest,Ztest+z0dinput.geo.z0 * vu,'color','r','marker','o','linestyle','none');
zplotprof(gca,z0dinput.cons.temps,z0dinput.exp0d.Rsepa,z0dinput.exp0d.Zsepa+z0dinput.geo.z0 * vu,'color','b','marker','none','linestyle','-');
axis('square')
axis('equal')


