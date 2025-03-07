%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
load testsinbad
option.init=0;
dis=[];
dist=[];
xv=[];
ev=[];
mat=[];
s2d_1=[];
s2d_2=[];
s2d_3=[];
s2d_4=[];
s2d_5=[];
s2d_6=[];
s2d_7=[];
s2d_8=[];
s2d_9=[];
s2d_10=[];
s2t=[];
divers=[];
%
% mis en forme des donnees
%
plas.temps = temps;
plas.R0=R0;
plas.a = a;
plas.Z=Z;
plas.e0=e0;
plas.e1=e1;
plas.d0=d0;
plas.tri=tri;
plas.Ip=Ip;
plas.mass=mass;
plas.charg=charg;
plas.fluxn=fluxn;
plas.rho = rho;
plas.Ti = Ti;
plas.Te = Te;
plas.ne = ne;
plas.n1=n1;
plas.n2=n2;
plas.n3=n3;
plas.n4=n4;
plas.Zeff = Zeff;
fais.pin1=pin1;
fais.en1=en1;
fais.fr1=fr1;
fais.A1=A1;
fais.Z1=Z1;
fais.align1=align1;
fais.pin2=pin2;
fais.en2=en2;
fais.fr2=fr2;
fais.A2=A2;
fais.Z2=Z2;
fais.align2=align2;
fais.type=type;
fais.neu=neu;
fais.rext=rext;
fais.cs=cs;
plas.Bcron=Bcron;
plas.Boutcron = Boutcron;
plas.BPcron=BPcron;
plas.thetacron = thetacron;
plas.xcron = xcron;
plas.Rcron = Rcron;
para.rot=rot;
para.option=option;
plas.equis=equis;
para.genes=genes;
plas.geo = geos;
plas.rext=rext;
sort.dis = dis;
sort.dist = dist;
sort.xv=xv;
sort.ev=ev;
sort.mat = mat;
sort.s2d_1 = s2d_1;
sort.s2d_2 = s2d_2;
sort.s2d_3 = s2d_3;
sort.s2d_4 = s2d_4;
sort.s2d_5 = s2d_5;
sort.s2d_6 = s2d_6;
sort.s2d_7 = s2d_7;
sort.s2d_8 = s2d_8;
sort.s2d_9 = s2d_9;
sort.s2d_10 = s2d_10;
sort.s2t=s2t;
sort.divers = divers;
   sor=nbi_sinbad_2temps_saturne(plas,fais,para,sort);
sor1=sor.DEP1;
sor2=sor.DEP2;
sor3=sor.DEP3;
sor4=sor.DEP4;
sor5=sor.DEP5;
sor6=sor.DEP6;
pion1=sor.pion1;
pion2=sor.pion2;
pion3=sor.pion3;
pion4=sor.pion4;
pelec=sor.pelec;
ploss=sor.ploss;
Iidn=sor.Jidn;
RAD=sor.RAD;
wpar=sor.Wpar;
wtot=sor.W;
wth=sor.Wth;
ang=sor.ang;
dis=sor.dis;
dist=sor.dist;
xv=sor.xv;
ev=sor.ev;
neu1=sor.neu1;
neu2=sor.neu2;
neu3=sor.neu3;
neu4=sor.neu4;
neu5=sor.neu5;
neu6=sor.neu6;
mat=sor.mat;
s2d=sor.s2d;
s2t=sor.s2t;
divers=sor.divers;






output.sor1            =  sor1;
output.sor2            =  sor2;
output.sor3            =  sor3;
output.sor4            =  sor4 ;
output.sor5            =  sor5 ;
output.sor6            =  sor6 ;
output.pion1           = pion1 ;
output.pion2           = pion2 ;
output.pion3           = pion3 ;
output.pion4           = pion4;
output.pelec           = pelec;
output.ploss           = ploss ;
output.Jidn            = Jidn;
output.RAD             = RAD ;
output.wpar            = wpar ;
output.wtot            = wtot;
output.wth             = wth  ;
output.ang             = ang ;
output.dis             = dis ;  
output.dist            = dist ;  
output.xv              = xv ;  
output.neutron.thth    = neu1 ;  
output.neutron.beth    = neu2 ;  
output.neutron.bebe    = neu3 ;  
output.neutron.itot    = neu4 ;  
output.neutron.abeb    = neu5 ;  
output.neutron.atot    = neu6 ;  
output.ev              = ev ;  
output.divers          = divers;
output.mat             = mat;
output.s2d_1           = s2d_1;
output.s2d_2           = s2d_2;
output.s2d_3           = s2d_3;
output.s2d_4           = s2d_4;
output.s2d_5           = s2d_5;
output.s2d_6           = s2d_6;
output.s2d_7           = s2d_7;
output.s2d_8           = s2d_8;
output.s2d_9           = s2d_9;
output.s2d_10          = s2d_10;
output.s2t             = s2t;

memoire.data.output    = output;


  pion   =  pion1(:,end)' + pion2(:,end)' +pion3(:,end)' +pion4(:,end)' + ploss(:,end)';
  pel    =  pelec(:,end)';
  js     =  Jidn(:,end)';
  nes    =  Z1(1) .* (sor1(:,end)'+sor2(:,end)'+sor3(:,end)') + ...
          Z2(1) .* (sor4(:,end)'+sor5(:,end)'+sor6(:,end)');
	  
  par1   =  (sor1(:,end)'+sor2(:,end)'+sor3(:,end)');
  par2   =  (sor4(:,end)'+sor5(:,end)'+sor6(:,end)');	  

  psupra =  wtot(:,end)' .* 1e6 - wpar(:,end)' .* 1e6;
  paniso =  wpar(:,end)' .* 1e6 - (1/2) .* psupra; 

% 03/02/2005 :  correction provisoire  de la puissance de bord pour simulation ITER
% la renormalisation n'est pas utile
%Vpr_loc   = interp1(rho,equi.vpr,RAD);
%Pel_loc   = equi.rhomax .* trapz(RAD,Vpr_loc .* pel,2);
pel(end) = 0;
%Pion_loc  = equi.rhomax .* trapz(RAD,Vpr_loc .* pion,2);
pion(end) = 0;


% reechantillonage du  profil
  sortie.el           = zbornes(interp1(RAD,pel,rho,'linear'),0,inf,0);
  %Pel_fin             = equi.rhomax .* trapz(rho,equi.vpr .* sortie.el,2);    % ajout du 03/02/2005
  %sortie.el           = sortie.el .* Pel_loc ./ Pel_fin;    % ajout du 03/02/2005
  
  sortie.ion          = zbornes(interp1(RAD,pion,rho,'linear'),0,inf,0);
  %Pion_fin            = equi.rhomax .* trapz(rho,equi.vpr .* sortie.ion,2);    % ajout du 03/02/2005
  %sortie.ion          = sortie.el .* Pion_loc ./ Pion_fin;    % ajout du 03/02/2005

  sortie.ne           = zbornes(interp1(RAD,nes,rho,'linear'),0,inf,0);
  sortie.j            = zbornes(interp1(RAD,js,rho,'linear'),0,inf,0);
  sortie.psupra       = zbornes(interp1(RAD,psupra,rho,'linear'),0,inf,0);
  sortie.paniso       = zbornes(interp1(RAD,paniso,rho,'linear'),-inf,inf,0);
%sortie.neutron.thth = zbornes(interp1(RAD,neu1(:,end),rho,'linear'),-inf,inf,0);
%sortie.neutron.beth = zbornes(interp1(RAD,neu2(:,end),rho,'linear'),-inf,inf,0);
%sortie.neutron.bebe = zbornes(interp1(RAD,neu3(:,end),rho,'linear'),-inf,inf,0);
%sortie.neutron.itot = zbornes(interp1(RAD,neu4(:,end),rho,'linear'),-inf,inf,0);
%sortie.neutron.abeb = zbornes(interp1(RAD,neu5(:,end),rho,'linear'),-inf,inf,0);
%sortie.neutron.atot = zbornes(interp1(RAD,neu6(:,end),rho,'linear'),-inf,inf,0);
% puissance de chauffage IDN
  paddidn             = zintvol(sortie.el+sortie.ion,gene.x,equi.vpr,equi.rhomax);
  pmaxidn            = sum(sum(pcons,2))/2;
