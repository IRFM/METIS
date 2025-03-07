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
function z0dinput = z0separatrix(z0dinput,option,plotonoff,from_gui)

if nargin == 1
  option = evalin('base','sepa_option');
end
z0dinput.sepa_option = option;
if nargin < 3
  plotonoff = 1;
end
if nargin > 3
    sepa  = z0dsepanew2(z0dinput.cons.temps,option,1);
else
    sepa  = z0dsepanew2(z0dinput.cons.temps,option);
end
if isempty(sepa)
    disp('creation of new LCFS canceled')
    return;
end
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
geo.K    = (abs(trapz(xu,sepa.Z .*  dRdx,2) ./ pi ./ geo.a) + (zmax - zmin)) ./ 3 ./ geo.a;

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
geo.d =  0.6 .* d + 0.4  .* sin(tm);


if ~isfield(option,'ton') | ~isfield(option,'toff')

	z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
	z0dinput.geo.z0      = geo.z0;     % centre geometrique du plasma en Z (m)
	z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
	z0dinput.geo.K       = geo.K;     % elongation (b/a)
	z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
	switch sepa.update_b0
	case 'on'
	  z0dinput.geo.b0      = sepa.b0;  % champ toroidal a R0
	end

	z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)

elseif isempty(option.ton) | isempty(option.toff)
	z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
	z0dinput.geo.z0      = geo.z0;     % centre geometrique du plasma en Z (m)
	z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
	z0dinput.geo.K       = geo.K;     % elongation (b/a)
	z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
	switch sepa.update_b0
	case 'on'
	  z0dinput.geo.b0      = sepa.b0;  % champ toroidal a R0
        end

	z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)

else
	% initilaisation
	z0dinput.exp0d.Rsepa = NaN .* sepa.R;
	z0dinput.exp0d.Zsepa = NaN .* sepa.Z;

	% on integre su theta directement
	t  = asin(max(0,min(1,z0dinput.geo.d)));
	vt = ones(size(t));
	R  = z0dinput.geo.R *vc + (z0dinput.geo.a * vc) .* cos(uc + (t* vc) .* sin(uc));
	Z  = ((z0dinput.geo.a .* z0dinput.geo.K) * vc)  .* sin(uc);

	fprintf('Geo_mode:')
	for k=1:size(sepa.R,1)
		if (z0dinput.cons.temps(k) >= option.ton) & (z0dinput.cons.temps(k) < option.toff)
			z0dinput.geo.R(k)       = geo.R(k);      % grand rayon du plasma (m)
			z0dinput.geo.z0(k)      = geo.z0(k);     % centre geometrique du plasma en Z (m)
			z0dinput.geo.a(k)       = geo.a(k);      % petit rayon du plasma (m)
			z0dinput.geo.K(k)       = geo.K(k);     % elongation (b/a)
			z0dinput.geo.d(k)       = geo.d(k);    % triangularite haute (definition entree de helena)
			switch sepa.update_b0
			case 'on'
			    z0dinput.geo.b0(k)      = sepa.b0(k);  % champ toroidal a R0
			end
			z0dinput.exp0d.Rsepa(k,:) = sepa.R(k,:);       % vecteur R des points de la separatrice (m)
			z0dinput.exp0d.Zsepa(k,:) = sepa.Z(k,:) - geo.z0(k) * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)
			fprintf('S');
		else
			z0dinput.exp0d.Rsepa(k,:) = R(k,:);
			z0dinput.exp0d.Zsepa(k,:) = Z(k,:);
			fprintf('M');
		end
	end
	fprintf('\n');


end



t  = asin(max(0,min(1,z0dinput.geo.d)));
u  = linspace(0,2.*pi,201);
vu = ones(size(u));
Rtest  = z0dinput.geo.R *vu + (z0dinput.geo.a * vu) .* cos(ones(size(geo.R,1),1) * u + t * sin(u));
Ztest = (z0dinput.geo.a .* z0dinput.geo.K) * sin(u);

if plotonoff == 0
  return
end

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


