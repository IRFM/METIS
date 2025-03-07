% script de verification des lois d'echelles
% les W de lois d'echelles
% loi d'echelles pour le calcul des  profil initiaux 
zs    = post.zerod;
cons  = post.z0dinput.cons;
geo   = post.z0dinput.geo;
t     = data.gene.temps;
x     = param.gene.x;
R     = data.equi.raxe(:,end);
a     = data.equi.a(:,end);
ep    = a ./ R;
K     = data.equi.e(:,end);
Ka    = data.equi.rhomax .* trapz(x,data.equi.vpr,2) ./ (2 .* pi .^ 2 .* R .* a .^ 2);
b     = a .* K ;
ip    = data.equi.ip;
Bt    = data.geo.b0;
zeff  = trapz(x,data.prof.zeff,2);
nbar  = trapz(x,data.prof.ne,2);
ne    = 10 .* nbar;
qcyl  = 5 .* Ka .* a .^ 2 .* Bt ./ R ./ (ip./1e6);
Fq    = data.equi.q(:,98) ./ qcyl;
ploss  = data.equi.rhomax .* trapz(x,data.equi.vpr .* (data.source.totale.el + data.source.totale.ion + (2/3) .* data.source.prad),2);
%ploss = data.gene.ploss;
% calcul de meff
meff = zeros(size(t));
for k = 1:length(t)
   nions   = squeeze(data.impur.impur(k,:,:))';  % size(nions) = [nbg,nbrho] !
   meff(k,1) = trapz(x, (nions(1,:) .* param.compo.a(1) + nions(2,:) .* param.compo.a(2) + nions(3,:) .* param.compo.a(3)  + ...
        nions(4,:) .* param.compo.a(4) + nions(5,:) .*param.compo.a(5)) ./  (nions(1,:) .* param.compo.z(1) + nions(2,:) .* param.compo.z(2) + nions(3,:) .* param.compo.z(3)  + ...
        nions(4,:) .* param.compo.z(4) + nions(5,:) .*param.compo.z(5)));  
end

% figure pour la geometrie et les parametres
h = findobj(0,'type','figure','tag','z0plotctrl');
if isempty(h)
       h=figure('tag','z0plotctrl');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(8,1,1)
plot(t,R,'ro',t,a,'bo',t,data.geo.r0,'r+',t,data.geo.a,'b+',t,geo.R,'mx',t,geo.a,'cx');
title('Parametres des lois d''echelles');
ylabel('R0 & a')
subplot(8,1,2);
plot(t,ep,'ro',t,K,'bo',t,Ka,'go',t,data.geo.a./data.geo.r0,'r+',t,data.geo.e1,'b+',t,geo.a./geo.R,'mx',t,geo.K,'cx',t,geo.vp ./ (2*pi^2.*geo.R.*geo.a.^2),'gx');
ylabel('ep, K & Ka')
subplot(8,1,3);
plot(t,zeff,'ro',t,data.cons.zeffm,'b+',t,zs.zeff,'mx');
ylabel('<Zeff>')
subplot(8,1,4);
plot(t,nbar./1e20,'ro',t,data.gene.nbar./1e20,'b+',t,cons.nbar./1e20,'mx');
ylabel('Nbar')
subplot(8,1,5);
plot(t,ploss./1e6,'ro',t,zs.ploss./1e6,'bx');
ylabel('Ploss')
subplot(8,1,6);
plot(t,Bt,'ro',t,geo.b0,'bx');
ylabel('Bt')
subplot(8,1,7);
plot(t,ip./1e6,'ro',t,zs.ip./1e6,'bx');
ylabel('ip')
subplot(8,1,8);
plot(t,meff,'ro',t,zs.meff,'bx');
ylabel('Meff')
xlabel('time (s)');
joint_axes(h,8);

% scaling
ip = ip ./ 1e6;
nbar = nbar ./ 1e20;
ne = nbar .* 10;
ploss = ploss ./ 1e6;

% ITERH-96P(th)        
taul96p  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* ploss .^ -0.73 .* ...
            R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s    
% cordey core  ref : NF 43 (2003), 670-674 core ift 2 equation 8 (sans dependance en beta)
taulcordey8  = 0.151 .* ip .^ 0.68 .* R .^ 2.32 .* ne .^ 0.59 .* Bt .^ 0.13 .*  ...
             Ka .^ -0.34 .* ep .^ 1.96 .* meff .^ 0.34  .* ploss .^ (0.42 -1);

% ITERH-98P(y,2)        
tauh98y2   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
            R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s  

%'Cordey all'
taupied  =  0.000643  .* ip .^ 1.58 .* R .^ 1.08 .* ne .^ (-0.08) .* Bt .^ 0.06 .* ....
               Ka .^ 1.81 .* ep .^ (-2.13) .* meff .^ 0.2 .* Fq .^ 2.09 .* ploss .^ (0.42 - 1);
tauhcordeyall     =  taupied + taulcordey8;
 
%'Cordey type 1'
taupied  =  0.00807  .* ip .^ 1.41 .* R .^ 1.37 .* ne .^ (-0.15) .* Bt .^ 0.32 .* ....
               Ka .^ 1.21 .* ep .^ 0.0 .* meff .^ 0.2 .* Fq .^ 1.61 .* ploss .^ (0.5 - 1);
tauhcordey1     =  taupied + taulcordey8;

% figure pour le contenu en energie
% figure pour la geometrie et les parametres
h = findobj(0,'type','figure','tag','z0plotctrl2');
if isempty(h)
       h=figure('tag','z0plotctrl2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

plot(t,data.gene.wdia./1e6,'or',t,data.gene.wth./1e6,'+r',t,zs.w./1e6,'mo',t,zs.wth./1e6,'m+',t,ploss .* taul96p,'bo',t,ploss .* taulcordey8,'b+', ...
      t,ploss .* tauh98y2,'co',t,ploss .* tauhcordeyall,'c+',t,ploss .* tauhcordey1,'cx');

legend('Wdia cronos','Wth cronos','Wdia 0D','Wth 0D','L96P','Core 8','98Y2','All','Type 1');
xlabel('time (s)')
ylabel('MJ');
title('Contenu en energie');


% verification du parametrage de kiauto
parametre = param.cons.coefa;
switch parametre.sccore
case 'ITER96L'
   disp('ITERH-96P(th)')        
otherwise
   disp('cordey core  ref : NF 43 (2003), 670-674 core ift 2 equation 8 (sans dependance en beta)')
end

switch parametre.scpied
case  'IPB98(y,2)'
   disp('ITERH-98P(y,2)');        
case 'Cordey all'
   disp('Cordey all');
otherwise
   disp('Cordey type 1');
end
