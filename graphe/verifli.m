% figure pour verification du li
h = findobj(0,'type','figure','tag','li');
if isempty(h)
       h=figure('tag','li');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

numchoc=param.from.shot.num;
if strcmp(param.from.machine,'TS')
  [diamgene,tligene]=tsbase(fix(numchoc),'sdiam');
  [beligene,tligene]=tsbase(fix(numchoc),'sbeli');
  ligene = 2.* (beligene-diamgene);
  [lipolo,tlipolo]=tsbase(fix(numchoc),'gqbeli%3');
  [diammag,tlimag]=tsbase(fix(numchoc),'gdiammg%1');
  [belimag,tlimag]=tsbase(fix(numchoc),'sbelimg');
  limag = 2.* (belimag-diammag);
  [diammagcor,tlimagcor]=tsbase(fix(numchoc),'sdiamgcor');
  [belimagcor,tlimagcor]=tsbase(fix(numchoc),'sbelicor');
  limagcor = 2.* (belimagcor-diammagcor);
  [liefit,tliefit]=tsbase(fix(numchoc),'sefli');
end

if ~exist('jeux1')
  chargedata
end
if isempty(jeux1)
  jeux1.data=[];
end

if isfield(jeux1.data,'equi')
plot(data.gene.temps,data.equi.li,jeux1.data.gene.temps,jeux1.data.equi.li,tlipolo,smooth(lipolo,10),tlimag,limag,tlimagcor,limagcor,tliefit,liefit,tligene,smooth(ligene,10))
xlabel('time (s)')
ylabel('li')
title(int2str(numchoc))
if ~isempty(limagcor)
  legend('CRONOS','CRONOS ref.','DPOLO','TMAG','TMAGCOR','EFIT','GENE');
else
  legend('CRONOS','CRONOS ref.','DPOLO','TMAG','EFIT','GENE');
end
else
plot(data.gene.temps,data.equi.li,tlipolo,smooth(lipolo,10),tlimag,limag,tlimagcor,limagcor,tliefit,liefit,tligene,smooth(ligene,10))
xlabel('time (s)')
ylabel('li')
title(int2str(numchoc))
if ~isempty(limagcor)
   legend('CRONOS','DPOLO','TMAG','TMAGCOR','EFIT','GENE');
else
   legend('CRONOS','DPOLO','TMAG','EFIT','GENE');
end
end
