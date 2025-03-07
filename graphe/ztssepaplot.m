% plot les separtrice de TS (dpolo,dmag efit et cronos)
function ztssepaplot	    
   
   
% donnees de cronos   
R     = double(evalin('base','data.equi.R'));
Z     = double(evalin('base','data.equi.Z'));
Rext  = double(evalin('base','data.geo.R'));
Zext  = double(evalin('base','data.geo.Z'));
temps = evalin('base','data.gene.temps');	    
try
   paroi = evalin('base','param.from.paroi');
catch
   paroi =[];
end

% le numero  du choc associe
numchoc = evalin('base','param.from.shot.num');

% donnees de la base    
% dpolo
[grho,tgrho,void,cert]   = tsbase(fix(numchoc),'grho');
tgrho          = tgrho(:,1);
% attente d'acces base
zaxe           = 0;
raxe           = 2.42;
alpha          = (0:15:345) ./ 180 .* pi;
vt             = ones(size(grho,1),1);
rr             = raxe + grho .* cos(vt*alpha);
zz             = zaxe + grho .* sin(vt*alpha);

% dmag & tmag
[rmg,tmg,amg,cert]   = tsbase(fix(numchoc),'grplasmg');
[zmg,tmg,amg,cert]   = tsbase(fix(numchoc),'gzplasmg');
[zparoimg,rparoimg]        = tsbase(fix(numchoc),'sparoimg');


% efit
tse = tsefitgeo(fix(numchoc));
if ~isempty(tse)
   alef = linspace(0,2*pi);
   rhoef = ones(size(alef));
   vtef  = ones(size(tse.time));
   [ref,zef] = rhotheta2rz(vtef * rhoef,vtef * alef,tse);
else
   ref = [];
   zef = [];
end


h = findobj(0,'type','figure','tag','tssepa');
if isempty(h)
       h=figure('tag','tssepa');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',0.5,'color',[1 1 1])

leg ={};
zplotprof(gca,temps,Rext,Zext,'color','r');
leg{end+1}='separatrice cronos';
zplotprof(gca,temps,squeeze(R(:,end,:)),squeeze(Z(:,end,:)),'color','b');
leg{end+1}='DSMF cronos (helena)';
if ~isempty(tgrho)
   zplotprof(gca,tgrho,rr,zz,'color','k','linestyle','none','marker','o');
   leg{end+1}='separatrice dpolo';
end
if ~isempty(tmg)
   zplotprof(gca,tmg,rmg,zmg,'color','m','linestyle','none','marker','+');
   leg{end+1}='separatrice tmag';
end
if ~isempty(tse)
   zplotprof(gca,tse.time,ref,zef,'color','c');
   leg{end+1}='separatrice efit';
end
if ~isempty(paroi)
   hold on 
   plot(paroi.R,paroi.Z,'k')
   leg{end+1}='paroi apolo';
end
if ~isempty(rparoimg)
   hold on 
   plot(rparoimg,zparoimg,'g')
   leg{end+1}='paroi tmag';
end
xlabel('R (m)');
ylabel('Z (m)');
title(sprintf('Separtrices et paroi de Tore Supra, choc # %d',fix(numchoc)));
legend(leg);
  

