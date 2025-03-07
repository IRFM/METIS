nbth = 35;
nbx = 21;
r0  = 2.4;
a   = 0.72;
rb0 = 2.4 .* 3.8;
ip = 1;
li = 0.5;

%
x = linspace(0,1,nbx);
th = linspace(0,2.*pi,nbth);
spr = 2 .* pi .* a  .* x;
vpr = 4 .* pi .^ 2 .* r0 .* a  .* x;
C2  = 4 .* pi .^ 2 .* a .* x ./ sqrt(r0 .^ 2 - a .^2 .* x .^ 2); 
% calcul de psi initial (cf. maple)
piqj = max(0.1,min(10,(exp(li)-1.65)./0.89));
jini = (1 - x.^2) .^ piqj;
iini = a .* trapz(x,spr .* jini);
jini = jini .* ip ./ iini;
mu0 = 4 .* pi .* 1e-7;
interi = -  a .* cumtrapz(x,mu0 ./r0 .* vpr .*jini);
inter = interi ./ max(eps,C2);
inter(1) = 0;
psi_init =  a .* cumtrapz(x,inter);
psi_init =  psi_init - psi_init(end);

volume  = 2 .* pi .^ 2 .* r0 .* a .^ 2 .* x .^2;
surface = pi .* a .^ 2 .* x .^2;
phi = prof.x .^ 2  .* pi .* a .^ 2 .* rb0 ./ r0;

% donnees de metis 
prof.x    = x;
prof.kx   = ones(size(prof.x));     
prof.dx   = zeros(size(prof.x));   
prof.Raxe = r0 .* ones(size(prof.x));    
prof.psi  = psi_init;
prof.fdia = rb0 .* ones(size(prof.x));
prof.jmoy = jini;
prof.ptot = eps .* ones(size(prof.x));
prof.rmx  = sqrt(prof.x .^ 2  .* pi .* a .^ 2 .* rb0 ./ r0) ./ sqrt(pi .* rb0 / r0);
%
geo_input.a     = a;
geo_input.R     = r0;
geo_input.K     = 1;
geo_input.d     = 0;
geo_input.b0    = rb0 ./ r0; 
geo_input.z0    = 0;
geo_input.sp    = pi .* a .^2;
geo_input.vp    = 2 .* pi .* r0 .* geo_input.sp;
geo_input.sext  = 4 .* pi .^ 2 .* r0 .* a;
geo_input.Rsepa = [];
geo_input.Zsepa = []; 
z0dstruct.z0dinput.option.morphing = 0;
 %
phys.mu0 = 4 .* pi .* 1e-7;
[profil,deuxd,moment,scalaire,facteur] = metis2equi1t(prof,geo_input,phys,nbth,1, ...
                                             z0dstruct.z0dinput.option.morphing,0);
