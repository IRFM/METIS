% calcul rapide avec extraction de point pertinent
function [zs,info,profli]= zerodghost(option,cons,geo,exp0d,thyb)

info = zero1t;
noms = fieldnames(info);
t    = cons.temps;
for k =1:length(noms);
	zs.(noms{k}) = NaN .* t;
end
zs.temps = t;

lp = z0dprofinfo;
noms = fieldnames(lp);
vp = NaN .* (t * ones(1,21));
for k =1:length(noms);
	profli.(noms{k}) = vp;
end
profli.temps = t;
profli.xli = linspace(0,1,21);


% pour la puissance de fusion
zs.pin = cons.plh + cons.pnbi + cons.picrh + cons.pecrh;
zs.pnbi = cons.pnbi;
zs.pbni_th = cons.pnbi;
zs.picrh  = cons.picrh;
zs.picrh_th  = cons.picrh;
zs.pecrh = cons.pecrh;
zs.plh = cons.plh;
zs.plh_th = cons.plh;
zs.zeff = cons.zeff;
zs.vp   = 2 .* pi .* geo.R .* pi .* geo.a .^ 2 .* geo.K;
zs.sp    = pi .* geo.a .^ 2 .* geo.K; 
zs.sext  = 2 .* pi .* geo.R .* pi .* geo.a  .* sqrt(1 + geo.K .^2); 
zs.peri = pi .* geo.a  .* sqrt(1 + geo.K .^2); 
zs.ane  = 0.1 .* ones(size(t));
zs.ate  = 1 .* ones(size(t));
zs.ape  = 1.1 .* ones(size(t));
zs.nebord =  5e-21 .* cons.nbar .^2;
zs.nem    =  cons.nbar;
zs.ne0    = (1+ zs.ane) .* zs.nem;
zs.nbar   = cons.nbar;
zs.prad   =  real(option.frad .*  1e6 ./ 4.5  .* (zs.zeff - 1) .* (cons.nbar./1e20) .^ 1.89 .* zs.sext .^ 0.94 ./ option.zmax .^ 0.12);
plossl2h  = 0.042 .* (zs.nbar./1e20) .^ 0.73 .* geo.b0 .^0.74 .* zs.sext .^0.98;
zs.ploss  = zs.pin;
zimp = (1- option.rimp) .* option.zimp + option.rimp .* option.zmax;
zs.nhem  = 0 .* t;
zs.pfus  = 0 .* t;
zs.ip = cons.ip;

% boucle de convergence
for k = 1:10   

   zs.tauthl  = 23e-3  .* (zs.ip./1e6) .^ 0.96 .* geo.b0 .^ 0.03 .* (zs.nbar ./ 1e19) .^ 0.4 .* ...
               (zs.ploss./1e6) .^ -0.73 .* geo.R .^ 1.83 .* geo.K .^ 0.64 .* (geo.a ./ geo.R) .^ -0.06 .* 2.5 .^ 0.2; 
	          
   zs.tauh   = 56.2e-3  .* (zs.ip./1e6) .^ 0.93 .* geo.b0 .^ 0.15 .* (zs.nbar ./ 1e19).^ 0.41 .*  ...
               (zs.ploss ./ 1e6) .^ -0.69 .* geo.R .^ 1.97 .* geo.K .^ 0.78 .* (geo.a ./ geo.R)  .^ 0.58 .* 2.5 .^ 0.19;    % s    

	       
   zs.modeh  = (zs.pin > plossl2h);
   zs.modeh  = max(zs.modeh,cat(1,0,zs.modeh(1:end-1)));       
   zs.taue = zs.tauh .* zs.modeh + zs.tauthl .* (~zs.modeh);
   zs.wth  = zs.taue.* zs.ploss;
   e = 1.602176462e-19;
   % hypothese ni = ne et Ti = Te;
   zs.tem = zs.wth ./ zs.vp ./ e ./ zs.nem ./ 2 ./ 1.5; %eV (pour 8.9 attendu)
   zs.te0 = (1+zs.ate) .* zs.tem;
   zs.n1m    = (zs.nem .* (zimp - zs.zeff) + (4 - 2 .* zimp) .* zs.nhem) ./ (zimp - 1); 
   zs.nDm    = zs.n1m ./ (cons.iso + 1);
   zs.nTm    =  max(0,zs.n1m -zs.nDm);

   
    sfus    = 1.1e-24 .* (zs.tem./1e3) .^ 2 .* zs.nDm .* zs.nTm .* zs.vp;
    zs.pfus = 0.5 .* (zs.pfus + (e .* 3.55e6) .* sfus);
    zs.nhem = 0.5 .* (zs.nhem  + sfus .* zs.taue .* 5 ./ zs.vp);
    zs.pin = cons.plh + cons.pnbi + cons.picrh + cons.pecrh + zs.pfus;
   
    zs.nimpm  = max(0, zs.nem - 2 .* zs.nhem - zs.n1m) ./zimp;
 
    zs.pbrem  = 5.35e-37 .* zs.nem .* sqrt(zs.tem./1e3) .* (zs.n1m + 4 .* zs.nhem + zimp .^ 2 .* zs.nimpm) .* zs.vp;
	
	
   betaf  = 2; % exposant de r/a dans Te(r/a)
   alphaf = zs.ate;
   te0    = zs.te0 ./ 1e3;
   ne0    = zs.ne0  ./ 1e20;
   gamma  = zs.ane;
   ka = (gamma + 3.87 .* alphaf + 1.46 ) .^ (-0.79) .* ...
         (1.98 + alphaf) .^ 1.36 .* betaf .^ 2.14 .* ...
         (betaf .^ 1.53 + 1.87 .* alphaf - 0.16) .^ (-1.33);

   % facteur de rapport d'aspect
   rap = max(1.5,min(15,geo.R ./ geo.a));
   ga  = 0.93 .* (1 + 0.85 .* exp( - 0.82 .* rap));
   % pa0
   pa0   = 6.04e3 .* geo.a .* ne0 ./ geo.b0;
   % puissance totale par rayon
   zs.pcyclo  = real(1e6 .* 3.84e-8 .* (1-abs(option.rw)) .^ 0.62 .* geo.R .* geo.a .^ 1.38 .* geo.K .^ 0.79 .* ...
          geo.b0 .^ 2.62 .* ne0 .^ 0.38 .* te0 .* (16 + te0) .^ 2.61 .* ...
         (1 + 0.12 .* te0 .* ((1-max(0,option.rw)) ./ pa0) .^ 0.41) .^ (-1.51) .* ...
          ka .* ga);
    
    zs.ploss  = zs.pin - option.fprad .*  zs.prad - zs.pcyclo - zs.pbrem;	
    zs.disrup = (zs.pin <= (2/3 .* zs.prad + zs.pcyclo + zs.pbrem));
end
% partie courant
eta = 2.8e-8 * ((zs.tem./1e3) .^ -1.5) % cette formule simple sous estime la resitivite
fact = 2; % facteur d'ajustement sur des simulations faites avec un code complet
zs.RR = fact .*2 .* geo.R ./ geo.a .^ 2 ./ geo.K .* eta .* zs.zeff
gcd = 0.35e19 * (zs.tem ./1e3);
zs.inbicd = gcd ./ geo.R ./ zs.nem .* zs.pnbi;
zs.ilh    = 3e19 ./ geo.R ./ zs.nem .* zs.plh;
zs.ieccd  = 0.5e19 ./ geo.R ./ zs.nem .* zs.pecrh;
zs.ifus  = 0.05e19 ./ geo.R ./ zs.nem .* zs.pfus;

mu0 = 4*pi*1e-7;
zs.betap = 8 * zs.wth ./ mu0 ./ 3 ./ (zs.ip * 1e6) .^ 2 ./ geo.R;
zs.iboot = 0.45 .* zs.ip .* sqrt(geo.a ./ geo.R) .* zs.betap;
zs.ini = zs.inbicd + zs.ilh + zs.ieccd + zs.iboot + zs.ifus;

zs.vloop = zs.RR .* (zs.ip - zs.ini);
zs.pohm = zs.vloop .^2 ./ zs.RR;

