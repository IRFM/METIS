% calcul le flux de neutron dd
% calcul de l'interraction faisceau-plasma
% attention vpr tel que trapz(x,vpr) =volume plasma
function neutron_beambeam = z0beambeam(temps,x,vpr,nep,tep,n1p,tip,zeffp,profnbi, ...
                            nDm,nTm,n1m,nhem,nimpm,nem,einj_a,einj_b,pnbi_a,pnbi_b,mu0a,mu0b,meme)

vt 	    = ones(size(temps));
ve          = ones(size(x));
pinj        = max(1,trapz(x,profnbi .* vpr,2));
volume      = max(1,trapz(x,vpr,2));
nD          = n1p .* ((nDm ./ max(1,trapz(x,vpr .* abs(n1p),2)) .* trapz(x,vpr,2)) * ve);
nDi         = max(1e13,trapz(x,nD .* profnbi .* vpr,2) ./ pinj);
tii         = max(30,trapz(x,tip .* profnbi .* vpr,2) ./ pinj);
tei         = max(30,trapz(x,tep .* profnbi .* vpr,2) ./ pinj);
nei         = max(1e13,trapz(x,nep .* profnbi .* vpr,2) ./ pinj);
zeff        = max(1,trapz(x,zeffp .* profnbi .* vpr,2) ./ pinj);

taus_nbi          = 6.27e8 .* 2 .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
fact              = (nDm ./2 + nTm ./ 3 + (n1m - nTm - nDm) + nhem  + nimpm ./2) ./ nem; 
ecrit_nbi_slow    = max(30,14.8 .* tei .* (2 .^ (3/2)  .* fact) .^ (2/3));   % c'est l'energie liee a vc
ecrit_nbi         = max(30,14.8 .* tei .* ( 2.* 2 .^ (1/2) .* zeff) .^ (2/3));   % c'est l'energir liee a % constante

mp          = 1.6726485e-27;            % masse au repos du proton (kg)
e           = 1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
vc          = sqrt(2 .*e .*  ecrit_nbi_slow ./ mp ./ 2);
vg          = sqrt(e .* ecrit_nbi./ mp ./ 2) ; %pour D seulement   
lx          = linspace(0,1,101);

% fonction de distribution A
sn0a         = max(1,pnbi_a ./ (e .* einj_a));
if length(sn0a) == 1
  	sn0a  = sn0a .* vt;
end
v0a          = sqrt(2 .*e .*  einj_a ./ mp ./ 2);

% fonction de distribution B
sn0b         = max(1,pnbi_b ./ (e .* einj_b));
if length(sn0b) == 1
 	sn0b  = sn0b .* vt;
end
v0b          = sqrt(2 .*e .*  einj_b ./ mp ./ 2);

% calcul du flux total de neutron 
% ref Physics reports , Tokamak plasma diagnostics based on measured neutron signals, 
% B. Wholle, volume 312, may 1999 + ref 13 
neutron_beambeam = zeros(size(temps));
% boucle sur les temps
for k=1:length(temps)
	if length(v0a) > 1
		lka = k;
	else
		lka = 1;
	end
	if length(v0b) > 1
		lkb = k;
	else
		lkb = 1;
	end
	if pinj(k) > (0.1 .* max(pinj))
		svmbb = 0;
		[va,mu,fvmua,fna,pn] = z0nbi_fvmu(sn0a(k),taus_nbi(k),v0a(lka),mu0a(k),vc(k),vg(k));
		%fva = 2 .* pi .* trapz(mu,fvmua,2);
		%nba = trapz(va,fva .* va .^ 2,1);
		%fna = fna ./ nba;
		%fvva  = sn0a(k) .* taus_nbi(k) ./ (va .^3 + vc(k) .^ 3) .* (va <= v0a) ./4 ./ pi;
		%nbva = 4 .* pi .* trapz(va,fvva .* va .^ 2,1);
		uva = ones(size(va));
		umu = ones(size(mu));	
		[vb,mu,fvmub,fnb,pn] = z0nbi_fvmu(sn0b(k),taus_nbi(k),v0b(lkb),mu0b(k),vc(k),vg(k));
		%fvb = 2 .* pi .* trapz(mu,fvmub,2);
		%nbb = trapz(va,fvb .* vb .^ 2,1);
		%fnb = fnb ./ nbb;
		%fvvb  = sn0b(k) .* taus_nbi(k) ./ (vb .^3 + vc(k) .^ 3) .* (va <= v0b) ./4 ./ pi;
		%nbvb = 4 .* pi .* trapz(va,fvvb .* vb .^ 2,1);
		inte   = zeros(size(va,1),size(vb,2));
		% boucle sur les polynome de legendre
                vamu = (va * mu);
		for n = 1:size(pn,1)	
			for l = 1:length(vb)
				%vrel  = sqrt((va * umu) .^ 2 + vb(l) .^ 2 - 2 .* (va * mu) .* vb(l));
				vrel  = sqrt(va(:,ones(size(mu))) .^ 2 + vb(l) .^ 2 - 2 .* vamu .* vb(l));
				ecm   = 0.5 .* 2 .* mp .* vrel .^ 2 ./ e; 
				cs = z0csdd(ecm./1e3) .* 1e-3 .* 1e-24 .* 1e-4;
                                pn_n = pn(n,:);
				inte(:,l)  = trapz(mu,pn_n(ones(size(va)),:) .* cs .* vrel,2);
				%inte(:,l)  = trapz(mu,(uva * pn(n,:)) .* cs .* vrel,2);
			end	
			svmbb = svmbb + 8 .* pi .^ 2 .* (1 ./ (2 .* (n - 1) + 1)) .* ...
				trapz(va,trapz(vb,inte .* (fna(:,n) * fnb(:,n)') .* (va * vb') .^ 2,2),1) ;
		end
		neutron_beambeam(k) = svmbb ./ (1 + meme);
		if length(temps) > 2
			fprintf('.');
		end
	else
		if length(temps) > 2
			fprintf('.');
		end
	end
end
if length(temps) > 2
	fprintf('\n');
end
% securite + le volume est compte deux fois
% R = na nb <sv> / (1 + dab)
% nombre = int(x,R*vpr) = Pinj ^ 2 / E ^ 2 / volume ^ 2  * f1a * f1b <sv> * Volume.
% dans les reactivite c'est des densite ! 
neutron_beambeam = (pinj > (0.1 .* max(pinj))) .* neutron_beambeam ./ volume;
