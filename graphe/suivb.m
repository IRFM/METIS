% sui.m
% pour suivre paramètres importants entre chocs en 
% salle de contrôle-commande.
% d'après batchdepro.m (Bécoulet) et suichoc.m (Hoang)
% Clarisse Bourdelle, 15/7/2002, poste 61-36

nchoc=list(k);

% lecture données générales
addpath /usr/drfc/martin/Pilotage
% Te(0)

[tesh,tsh]=tsbase(nchoc,'gshte');
tt=[0:(tsh(length(tsh)))./2000:tsh(length(tsh))]';
te0=max(tesh');
if isempty(tesh)==1
   te0 = zeros(1,length(tt))';
else   
   te0=interp1q(tsh,te0',tt);
end

[teth,tlid]=tsbase(nchoc,'gtethom%8');
if isempty(teth)==1
   teth = zeros(1,length(tt));
else   
   teth=interp1q(tlid,teth,tt);
end


% Pohmic

[Pohm,tohm]=tsbase(nchoc,'spohm'); 

Pohm=interp1q(tohm,Pohm,tt);

% Grand rayon

[R,tmg2]=tsbase(nchoc,'srmaj'); 
if isempty(R)==1
   R = zeros(1,length(tt));
else   
   R=interp1q(tmg2,R,tt);
end

% Petit rayon
 
[a,tmg2]=tsbase(nchoc,'samin');  
if isempty(a)==1
   a = zeros(1,length(tt));
else   
   a=interp1q(tmg2,a,tt);
end

% Ip

[Ipla,tmg2]=tsbase(nchoc,'sipmes'); 
if isempty(Ipla)==1
   Ipla = zeros(1,length(tt));
else   
   Ipla=interp1q(tmg2,Ipla,tt);
end

% Btor

[Btor,tmg2]=tsbase(nchoc,'sitor'); 
if isempty(Btor)==1
   Btor = zeros(1,length(tt));
else   
   Btor=interp1q(tmg2,Btor,tt);
end
Btor1=max(Btor)*3.07e-3;

% Wdia

[wdia,tmg2]=tsbase(nchoc,'swdia'); 
if isempty(wdia)==1
   wdia = zeros(1,length(tt));
else   
   wdia=interp1q(tmg2,wdia,tt);
end


% fuelling

[inj,tinj]=tsbase(nchoc,'sdebcal'); 
if isempty(inj)==1
   inj = zeros(1,length(tt));
else   
   inj=interp1q(tinj,inj,tt);
end


% Pfci

[Picrh,ticrh]=tsbase(nchoc,'spuiss'); 
if isempty(Picrh)==1
   Picrh = zeros(1,length(tt));
else   
   Picrh=interp1q(ticrh,Picrh,tt);
end

% Pecrh

[Pecrh1,tecrh]=tsbase(nchoc,'spia1'); 
if isempty(Pecrh1)==1
   Pecrh1 = zeros(1,length(tt));
else   
   Pecrh1=interp1q(tecrh,Pecrh1,tt);
end
[Pecrh2,tecrh]=tsbase(nchoc,'spia2'); 
if isempty(Pecrh2)==1
   Pecrh2 = zeros(1,length(tt));
else   
   Pecrh2=interp1q(tecrh,Pecrh2,tt);
end
Pecrh=Pecrh1+Pecrh2;

% Plh

[Plh,tlh]=tsbase(nchoc,'gphyb%3'); 
if isempty(Plh)==1
   Plh = zeros(1,length(tt));
else   
   Plh=interp1q(tlh,Plh,tt);
end

% fraction P rayonnée

[Prad,trad]=tsbase(nchoc,'sprad');
if isempty(Prad)==1
   Prad = zeros(1,length(tt));
else   
   Prad=interp1q(trad,Prad,tt);
end
if isempty(ticrh)==1 & isempty(tlh)==1,  
   frad=Prad./(Pohm);
elseif isempty(ticrh)==1,
   frad=Prad./(Pohm+Plh);
elseif isempty(tlh)==1, 
   frad=Prad./(Pohm+Picrh);
elseif isempty(trad)==1, 
   frad=zeros(1,length(tt));
else
   frad=Prad./(Pohm+Picrh+Plh);
end

% neutrons 

[Rnt,tnt]=tsbase(nchoc,'gfluntn%2');
if isempty(Rnt)==1
   Rnt = zeros(1,length(tt));
else   
   Rnt=interp1q(tnt,Rnt,tt);
end
NT=tsbase(nchoc,'rneutron');
NT=NT*1e12;

% faisceau de neutres

% beta_p

[betap,tmg2]=tsbase(nchoc,'sdiam'); 
if isempty(betap)==1
   betap = zeros(1,length(tt));
else   
   betap=interp1q(tmg2,betap,tt);
end

% bootstrap current fraction 

fboot=.2*sqrt(R./a).*betap;

% li

[li,tmg2]=tsbase(nchoc,'scoupli'); 
if isempty(li)==1
   li = zeros(1,length(tt));
else   
   li=interp1q(tmg2,li,tt);
end

 
% beta_n 

betan=(4*betap.*Ipla)./(a.*Btor*3.07e-3); 

% facteur H_98

[w1,tscal]=tsbase(nchoc,'gbilan%6');  
if isempty(w1)==1
   w1 = zeros(1,length(tt));
else   
   w1=interp1q(tscal,w1,tt);
end
[w2,tscal]=tsbase(nchoc,'gbilan%8');
if isempty(w1)==1
   w2 = zeros(1,length(tt));
else   
   w2=interp1q(tscal,w2,tt);
end
h98=w1./w2;


% Ti(0)

[ti0,tcx]=tsbase(nchoc,'stibrag'); 
[teb,tcx]=tsbase(nchoc,'stebrag'); 
if isempty(ti0)==1
  ti0 = zeros(length(tt),1);
  teb = zeros(length(tt),1);
else   
   ti0=interp1q(tcx,ti0,tt);
   teb=interp1q(tcx,teb,tt);
end

% Ti/Te (0)

if isempty(ti0)==0,
   tiste=ti0./te0;
   tisteth=ti0./teth;
else 
   tiste=zeros(1,length(tt));
   tisteth=zeros(1,length(tt));
end

% Zeff

[zeff,tzeff]=tsbase(nchoc,'szfbrm'); 
if isempty(zeff)==1
   zeff = zeros(1,length(tt));
else   
   zeff=interp1q(tzeff,zeff,tt);
end

% Dalpha

[elm,telm]=tsbase(nchoc,'ghalpha');
if isempty(elm)==1
   elm = zeros(1,length(tt));
else   
   elm=interp1q(telm,elm,tt);
end 

% densité

[nvol,tnvol]=tsbase(nchoc,'snmoy');
if isempty(nvol)==1
   nvol = zeros(1,length(tt));
else
   nvol=interp1q(tnvol,nvol,tt);
end
[nl,tnl]=tsbase(nchoc,'gnl%4'); 
if isempty(nl)==1
   nl = zeros(1,length(tt));
else
   nl=interp1q(tnl,nl,tt);
end
nvol=nvol/1e19; 
nl=nl/1e19; 

% niveau de Fer

% Vloop

[vloop,tvloop]=tsbase(nchoc,'svsur'); 
if isempty(vloop)==1
   vloop = zeros(1,length(tt));
else
   vloop=interp1q(tvloop,vloop,tt);
end

% injection de gaz
%gaz = tsbase(nchoc,'rgazchoc');
figure

plot(tt,ti0./teb,tt,ti0./smooth(te0,5),tpf,pf)
[xg,yg]=ginput(2);
ind1=find(tt>min(xg)&tt<max(xg));

rap1=mean(ti0(ind1)./teb(ind1))
rap2=mean(ti0(ind1)./smooth(te0(ind1),5))
Tesh = mean(te0(ind1))

addpath /usr/drfc/bourdel/outilsdata/analyz/bragg/
%cd bragg
visu2_bragg(nchoc);

frac=input('fraction de nD/ne , 0.6?, :');
[tN,Ti,Ntn,Wd,nd] = Ti_Ntn2( nchoc,frac);
figure
subplot(131)
hold on 
axis([0 max(tt) 0 3.5])
title(['Ti neutron vs Ti bragg (r) avec nD/ne = ', num2str(frac)])
plot(tN,Ti,tt,ti0,'r');
grid on
ylabel(['#TS' int2str(nchoc)])
xlabel('s')
frac=input('essai avec nouvelle fraction de nD/ne :');
[tN,Ti,Ntn,Wd,nd] = Ti_Ntn2 ( nchoc,frac );
subplot(132)
hold on 
axis([0 max(tt) 0 3.5])
title(['Ti neutron vs Ti bragg (r) avec nD/ne = ', num2str(frac)])
plot(tN,Ti,tt,ti0,'r');
grid on
ylabel(['#TS' int2str(nchoc)])
xlabel('s')
frac=input('essai avec nouvelle fraction de nD/ne :');
[tN,Ti,Ntn,Wd,nd] = Ti_Ntn2 ( nchoc,frac );
subplot(133)
hold on 
axis([0 max(tt) 0 3.5])
title(['Ti neutron vs Ti bragg (r) avec nD/ne = ', num2str(frac)])
plot(tN,Ti,tt,ti0,'r');
grid on
ylabel(['#TS' int2str(nchoc)])
xlabel('s')


filename2=['sui_',int2str(nchoc)];
eval(['save ' filename2 ' nvol tt nl  tiste zeff ti0 te0  vloop li fboot Btor Pohm  Plh  Picrh Pecrh wdia Ipla h98 elm betap betan NT Rnt frad Prad a R ']); 

