% script de comparaison de Te cronos et de Te ece pour TS
if isempty(post.ece.tece)
   disp('Pas de postprocessing pour ece')
   return
end
shts =  cgcgettrait(fix(param.from.shot.num),'tshetero','sh');
if isempty(shts.times)
   disp('Pas de donnees pour ece')
   return
end

ind = find(shts.tenv <0);
if ~isempty(ind)
   shts.tenv(ind) = NaN;
end
ind = find(shts.te <0);
if ~isempty(ind)
   shts.te(ind) = NaN;
end

h = findobj(0,'type','figure','tag','cmpece');
if isempty(h)
       h=figure('tag','cmpece');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',0.5,'color',[1 1 1])

subplot(2,1,1)
if isempty(shts.tenv)
   plot(shts.times,shts.te,'r-',data.gene.temps,post.ece.tece,'b-');
else
   plot(shts.times,shts.tenv,'r-',data.gene.temps,post.ece.tece./1e3,'b-');
end
xlabel('temps (s)')
ylabel('Te (ece = -,cronos = ., keV)');
title(sprintf('Comparaison ECE, choc TS # %g',param.from.shot.num));
zoom on

subplot(2,2,3)
if isempty(shts.tenv)
   zplotprof(gca,shts.times,shts.r,shts.te,'color','r','linestyle','none','marker','o');
else
   zplotprof(gca,shts.times,shts.r,shts.tenv,'color','r','linestyle','none','marker','o');
end
zplotprof(gca,data.gene.temps,post.ece.Rece,post.ece.tece ./ 1e3,'color','b','linestyle','-','marker','none');
xlabel('temps (s)')
ylabel('Te (keV)');
 
subplot(2,2,4)
if isempty(shts.tenv)
   temes = interp1(shts.times,shts.te,data.gene.temps,'nearest');
   
else
   temes = interp1(shts.times,shts.tenv,data.gene.temps,'nearest');
end
ll     = find(any(isfinite(post.ece.tece),1));
tef   = NaN .* temes;
for k =1:length(ll)
   pp = polyfit(post.ece.tece(:,ll(k)),temes(:,k),1);
   tef(:,k) = polyval(pp, post.ece.tece(:,ll(k)));
end
if isempty(shts.tenv)
   plot(shts.times,shts.te,'r-',data.gene.temps,tef,'b-');
else
   plot(shts.times,shts.tenv,'r-',data.gene.temps,tef,'b-');
end
xlabel('temps (s)')
ylabel('Te recale');
zoom on
