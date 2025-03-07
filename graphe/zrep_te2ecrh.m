% reponse plasma 
if isempty(post)
   return
end
if isempty(post.ece)
   return
end
numchoc = param.from.shot.num;
% lecture des rayons de l'ece
[recet,trecet,xrecet,crecet] = tsbase(fix(numchoc),'gshr');
[teecet,tteecet,xteecet,cteecet] = tsbase(fix(numchoc),'gshtenv');
if isempty(teecet)
	[teecet,tteecet,xteecet,cteecet] = tsbase(fix(numchoc),'gshte');
end
if isempty(teecet)
   return
end
h1 =findobj(0,'type','figure','tag','te2ecrh1');
if isempty(h1)
	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
  	          'times','defaultlinelinewidth',3,'tag','te2ecrh1','name','reponse plasma');
else
	  figure(h1);
end
clf

t= data.gene.temps;
x = param.gene.x;
Raxe = mean(squeeze(double(data.equi.R(:,2,:))),2);
indok = find(isfinite(Raxe));
Raxe = mean(Raxe(indok));
dd   = abs(Raxe - post.ece.Rece);
for k = 1:size(dd,2);
   indok = find(isfinite(dd(:,k)));
   if isempty(indok)
      d(k) = inf;
   else
      d(k) = mean(dd(indok,k));
   end
end 
indvp = find(d == min(d));
indvpr = input(sprintf('Indice pour post ? (%d) ',indvp));
if ~isempty(indvpr)
   indvp = indvpr;
end 


ind  = find(recet <=0);
if ~isempty(ind)
   recet(ind) = inf;
end
dd   = abs(Raxe - recet);
for k = 1:size(dd,2);
   indok = find(~isnan(dd(:,k)));
   if isempty(indok)
      d(k) = inf;
   else
      d(k) = mean(dd(indok,k));
   end
end 
indv = find(d == min(d));
indvpr = input(sprintf('Indice pour ece ? (%d) ',indv));
if ~isempty(indvpr)
   indv = indvpr;
end 

second = 0;
try
   if ~isempty(jeux1)
      if isfield(jeux1,'post')
         if  isfield(jeux1.post,'ece')
            if ~isempty(jeux1.post.ece.tece)
               second = 1;
            end
          end
      end
   end
end

ha1=subplot(2,1,1);
if second == 0
   plot(tteecet,teecet(:,indv),'.b', ...
      t,post.ece.tece(:,indvp)/1e3,'k')
   legend('mesures SH','Cronos BGB_T_S');
else
   plot(tteecet,teecet(:,indv),'.b', ...
      t,post.ece.tece(:,indvp)/1e3,'k', ...
      jeux1.data.gene.temps,jeux1.post.ece.tece(:,indvp)/1e3,'r')
   legend('mesures SH','Cronos BGB_T_S (spia)','Cronos BGB_T_S (sika)');
end
ylabel('Te (keV)')
xlabel('temps (s)')
grid on
title(sprintf('Reponse de Te0(Pecrh) @ %d', fix(numchoc)));

ha2=subplot(2,1,2);
if second == 0
   plot(t,data.gene.paddfce/1e3,'k')
else
   plot(t,data.gene.paddfce/1e3,'k', ...
   jeux1.data.gene.temps,jeux1.data.gene.paddfce/1e3,'r')
end
ylabel('P_E_C_R_H (kW)')
xlabel('temps (s)')
set(ha1,'xlim',get(ha2,'xlim'));
grid on
