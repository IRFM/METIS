function [Tisol,Wd] = Ti_Ntn2(data,param)

% function [tN,Ti,Ntn,Wd,nd] = Ti_Ntn2 ( Nch,nDne,delt,PiqT,PiqN );
%
%     Version CIEL : sans traitement de Densit�
%
%  Nch  : Numero du Choc
%  nDne : Proportion Deuterium par Electrons (Defaut : 0.64)
%  delt : Resolution temporelle en ms (Defaut : 16 ms)  ->  0 = Acquisition
%  PiqT : Piquage de Ti (Defaut 1.6)
%  PiqN : Piquage de la densit�(Defaut 0.6)
%
%  tN   : Vecteur temps du calcul ( <=> Suffisament de Neutrons )
%  Ti   : Temperature Centrale en keV
%  Ntn  : Production de Neutrons par seconde
%  Wd   : Energie contenue dans les ions Deuterium en kJ
%  nd   : Densit�centrale de Deut�ium en 10^19 m-3


Nch    = param.from.shot.num
Tisol  = zeros(size(data.prof.ti));
Wd     = zeros(size(data.gene.wth));

x      = param.gene.x;
dx     = mean(diff(x));
temps  = data.gene.temps;
[tibr,ttibr] = tsbase(Nch,'stibrag');
nPiqT        = 100;
PiqT   = linspace(0.5,8,nPiqT)';
vc     = ones(size(PiqT));
vl     = ones(size(x));

if Nch < 30000
  [N1,tn1] = tsbase(Nch,'gntnmd2'); 
  [N4,tn4] = tsbase(Nch,'gntnmd4');
  N4       = interp1(tn4(:,1),N4,temps);
  N1       = interp1(tn1(:,1),N1,temps);
  tn4      = temps;
  tn1      = temps; 
  Gain     = tsmat(0,'DNEUTRON;Calcul_Flux;GCOUPS'); 
  nN       = find(N4(:,2) > 10); 
  if length(nN)< 16,
    Erreur = 'Pas assez de Neutrons'
  end  
  n1 = min(nN); n2 = max(nN);
  N4 = N4(n1:n2,:); N1 = N1(n1:n2,:); tN = tn4(n1:n2);
  Errt = max(abs(tn1(n1:n2,1)-tN)); 
  if Errt > 1e-3,
    Erreur = 'Mauvaise base de Temps sur les Neutrons'
  end 
  Ntn = N4(:,1)*Gain(6); nN = find(N4(:,1) > 400); 
  Ntn(nN) = N1(nN,1)*Gain(1); nN = find(N1(:,1) > 1000); 
  Ntn(nN) = N4(nN,2)*Gain(7); nN = find(N4(:,2) > 1000); 
  Ntn(nN) = N1(nN,2)*Gain(2); Ntn = Ntn./[diff(tN);1];
  clear tn1 tn4 N1 N4 Gain
else
  [Ntn,tN]=tsbase(Nch,'gfluntn%2');
  Ntn = interp1(tN(:,1),smooth(10.^Ntn,5),temps);
end

for k=1:length(temps)
  if temps(k) > min(ttibr) & temps(k) < max(ttibr)
    indbr  = min(find(ttibr>temps(k)));
  else
    indbr = [];
  end
  
  if ~isempty(indbr)
  a      = data.geo.a(k)*100;
  R      = data.geo.r0(k)*100;
  r      = data.equi.a(k,:)*100;
  dr     = dx*a;
  indd   = find(param.compo.a == 2 & param.compo.z == 1);
  nd     = squeeze(data.impur.impur(k,:,indd))/1e6;
  ne     = data.prof.ne(k,:)/1e6;
  Eta2   = (nd(1)/ne(1))^2;

  nd0    = nd(1);

  nlsan0 = 2*sum(nd)*dx/nd0;
  clear Flux Wion

  Ti0    = tibr(indbr)*1.2;
  Ti     = Ti0*(1-vc*(x.*x)).^(PiqT*vl)+0.013;
  sv     = 2.33e-14*exp(-18.76./Ti.^0.33333)./Ti.^0.666667;
  Ntnvb  = 0.5*(vc*(nd.*nd)).*sv; 
  Flux   = sum(Ntnvb.*(vc*r),2)*4*pi*pi*R*dr; 
  Wi     = 1.5*(vc*nd).*Ti; 
  Wion   = sum(Wi.*(vc*r),2)*4*pi*pi*R*dr; 


  a      = a/100; 
  R      = R/100; 
  nd0    = nd0/1e13; 
  Wion   = Wion*param.phys.e;

  Rn  = data.geo.r0(k); 
  an  = data.geo.a(k,1);
  nli = data.gene.nbar(k)/1e19 *2 .* an;
  n0  = nli/nlsan0./an;
% Calcul du Flux normalise et de Ti
  delt = param.gene.dt*1e3;
  Fnor = Ntn*R./Rn*a*a./an./an*nd0*nd0./n0./n0/Eta2;
  indok  = find(isfinite(Fnor));
  Alpha = 1; 
  if delt > 0.5, Alpha = 1000*mean(diff(tN))/delt; end
  if Alpha < 0.8, Fnor(indok) = glisse(Fnor(indok),Alpha); end

  if isfinite(Fnor(k))
    iT = min(find(Flux-Fnor(k)<0));
    if ~isempty(iT)
      Tisol(k,:) = Ti(iT,:); 
      Wd(k) = Wion(iT); 
    end
  end
  end
end

 
