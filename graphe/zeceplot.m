% script de plot de l'ECE
if isempty(post)
   return
end
if isempty(post.ece)
   return
end
numchoc = fix(param.from.shot.num);
occur   = 10*(param.from.shot.num-numchoc)
% lecture des rayons de l'ece
[recet,trecet,xrecet,crecet] = tsbase(fix(numchoc),'gshr');
[teecet,tteecet,xteecet,cteecet] = tsbase(fix(numchoc),'gshtenv');
if isempty(teecet)
	[teecet,tteecet,xteecet,cteecet] = tsbase(fix(numchoc),'gshte');
end
if isempty(teecet)
   return
end
numstr       = [int2str(numchoc),'.',int2str(occur)];
eval(['[rhoprof,trhoprof,void1,void2]=tsbase(',numstr,',''gprofrhote'');']);
eval(['[teprof,tteprof,void1,void2]=tsbase(',numstr,',''gprofte'');']);
eval(['[tefitprof,ttefitprof,void1,void2]=tsbase(',numstr,',''gproftefit'');']);
eval(['[aprof,taprof,void1,void2]=tsbase(',numstr,',''sprofamin'');']);
eval(['[Rprof,tRprof,void1,void2]=tsbase(',numstr,',''sprofrmaj'');']);
eval(['[dprof,tdprof,void1,void2]=tsbase(',numstr,',''sprofd0'');']);

rhoprof = rhoprof(:,13:end);
[mr,imr]=min(rhoprof,[],2);
tetaprof=0*rhoprof;
for k=1:size(tetaprof,1)
  tetaprof(k,1:imr(k))=pi;
end
[Rprof,Zprof]=rho2rz(rhoprof,tetaprof,aprof,Rprof,0*Rprof,dprof);
teprof = teprof(:,13:end);
%rhoprof = interp1(trhoprof,rhoprof,data.gene.temps);
%teprof = interp1(tteprof,teprof,data.gene.temps);

% calcul de Rext et Rint
rout = squeeze(double(data.equi.R(:,end,:)));
Rext = max(rout,[],2);
Rint = min(rout,[],2);
Raxe = mean(squeeze(double(data.equi.R(:,2,:))),2);

% R_ECE
h1 =findobj(0,'type','figure','tag','zeceplot');
if isempty(h1)
   	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
     	          'times','defaultlinelinewidth',0.5,'tag','zeceplot','name','Rayon ece');
else
   	  figure(h1);
end
clf

cc=get(gca,'colororder');
plot(data.gene.temps,post.ece.Rece,'.');
hold on
set(gca,'colororder',cc);
plot(trecet,recet,'-');
plot(data.gene.temps,Rext,'g','linewidth',3);
plot(data.gene.temps,Rint,'g','linewidth',3);
plot(data.gene.temps,Raxe,'k-.');
hold off
xlabel('temps')
ylabel('R_e_c_e (m)');
title(sprintf('choc %d@%d  . -> Cronos, - -> Tshetero', ...
		              fix(numchoc),fix(10.01*(numchoc- fix(numchoc)))));
set(gca,'ylim',[1.5,3.2])
	
% Te_ECE
h1 =findobj(0,'type','figure','tag','zeceplot2');
if isempty(h1)
   	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
     	          'times','defaultlinelinewidth',0.5,'defaultlinemarkersize',0.5, ...
                'tag','zeceplot2','name','Te ece');
else
   	  figure(h1);
end
clf

cc=get(gca,'colororder');
plot(data.gene.temps,post.ece.tece/1e3,'-');
hold on
set(gca,'colororder',cc);
plot(trecet,teecet,'o');
hold off
xlabel('temps')
ylabel('Te_e_c_e (keV)');
title(sprintf('choc %d@%d  - -> Cronos, o -> Tshetero', ...
		              fix(numchoc),fix(10.01*(numchoc- fix(numchoc)))));
	
h1 =findobj(0,'type','figure','tag','zeceplot3');
if isempty(h1)
   	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
     	          'times','defaultlinelinewidth',2,'defaultlinemarkersize',6, ...
                'tag','zeceplot3','name','Te ece (profil)');
else
   	  figure(h1);
end
clf
zplotprof(gca,data.gene.temps,post.ece.Rece,post.ece.tece/1e3,'color','r');
zplotprof(gca,trecet,recet,teecet,'color','b','linestyle','o');
zplotprof(gca,tteprof,Rprof,teprof,'color','k');
zplotprof(gca,data.gene.temps,post.ece.Rece,post.ece.tece/1e3,'color','r','linestyle','*');
xlabel('R (m)')
ylabel('Te  (keV)');
legend('Cronos','TSbase','TPROF')
title(sprintf('choc %d  ',fix(numchoc)));
