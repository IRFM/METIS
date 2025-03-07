function [gtecr,gte,fcr] = gcrt(data,param)
%
% recalcul du gradient critique
% Syntaxe : [gtecr,gte,fcr] = gcrt(data,param);
%
% Auteur : V. Basiuk
% version 1.9 : 11 avril 2002
%
% dernieres modifications
%

prof      = data.prof;
equi      = data.equi;
t         = data.gene.temps;
x         = param.gene.x;
phys      = param.phys;
%
% regarde dans coefb puis coefa
%
parametre = param.cons.coefb;
if ~isfield(parametre,'ce_es')
  parametre = param.cons.coefa;
  if ~isfield(parametre,'ce_es')
    disp('pas de fonction avec gradient critique')
	 return
  end
end
for k=1:length(t)
	
  vth      = sqrt(2 .* prof.te(k,:) .* phys.e ./ phys.me);
  wc       = phys.e .* sqrt(equi.b2(k,:)) ./ phys.me;
  rhoe     = vth ./ wc;
  wpe      = sqrt(phys.e .^ 2 .* max(1e13,prof.ne(k,:)) ./ phys.me ./ phys.epsi0);

% profil de s/q : 
  sq       = abs(x ./ max(0.1,prof.q(k,:)) .^2 .* pdederive(x,prof.q(k,:),0,2,2,1));

%  longueur de gradient
%lte      = max(13.6,abs(prof.te)) ./ gte;
  lte       = prof.lte(k,:);

% gradient de Te 
%gte      = max(0.01, abs(prof.gte ./ equi.grho)); 
  gte(k,:) = max(13.6,abs(prof.te(k,:))) ./ lte;

% taux de croissance
  gamma    = vth ./ sqrt(lte .* equi.rmoy(k,:));

% longueur de correlation electrostatique
  lce      = rhoe .* equi.rmoy(k,:) ./ lte;

% longueur de correlation electromagnetique
  lcm      = phys.c ./ wpe;

% kie sans gradient critique
  kie      = prof.ne(k,:) .* gamma .* (parametre.ce_es .* lce .^ 2 + parametre.ce_em .* lcm .^2);

% dependance en q 
  if parametre.expoq ~= 0
	 kie      = kie .* equi.q(k,:) .^ parametre.expoq;
  end

% dependance en gradient critique
  pdiff   = equi.ptot(k,:) -prof.pe(k,:);
  mask    = (pdiff > 1e2);
  rapport = equi.ptot(k,:) ./ abs(pdiff+eps);
  rapport = rapport .* mask;
  gtecr(k,:)   = parametre.offset * 1000 + parametre.cl .*  ...
          (sq .* max(13.6,abs(prof.te(k,:))) ./ equi.rmoy(k,:)) .* rapport;
  fcr(k,:)     = max(0, 1 - gtecr(k,:) ./ gte(k,:)); 
%  if k==28
%	  keyboard
%  end	  
end

h = findobj(0,'type','figure','tag','gradient_critique');
if isempty(h)
       h=figure('tag','gradient_critique');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])

subplot(2,1,1)

plotprof(gca,t,x,gtecr,'linestyle','-','color','r')
plotprof(gca,t,x,gte,'linestyle','--','color','b')
ylabel('1/Lt')
legend ('1/Ltc','1/Lt')

subplot(2,1,2)

plotprof(gca,t,x,fcr,'linestyle','-','color','r')


