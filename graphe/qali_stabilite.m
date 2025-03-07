% plot du diagramme (qa,li)
function stab=qali_stabilite

qdif = evalin('base','data.prof.q');
qeq = evalin('base','data.equi.q');
lidif = evalin('base','data.gene.li');
lieq = evalin('base','data.equi.li');
param =evalin('base','param'); 


ind = find(~ isfinite(lieq));
if ~ isempty(ind)
   lieq(ind) = lidif(ind);
end
ind = find(~ isfinite(qeq));
if ~ isempty(ind)
   qeq(ind) = qdif(ind);
end

h = findobj(0,'type','figure','tag','qali');
if isempty(h)
       h=figure('tag','qali');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',6)
contour_qali = load('qali_stabilite');
plot(qeq(2:(end-1),end),lieq(2:(end-1)),'g',qdif(2:(end-1),end),lidif(2:(end-1)),'b');
stab.qeq     = qeq(2:(end-1),end);
stab.qdif    = qdif(2:(end-1),end);
stab.lieq    = lieq(2:(end-1));
stab.lidif   = lidif(2:(end-1));
stab.contour = contour_qali;
xlabel('q_a');
ylabel('li');
title(sprintf('Stability plot (q_a,li) for the shot %s #%d',param.from.machine,fix(param.from.shot.num)));
hold on
% pour les legend
 set(plot(NaN,NaN,'c+'),'linewidth',1);
sty = strvcat('ms','md','mv','m^','m<','m>','mp','mh');
for k = 1:5
    set(plot(NaN,NaN,sty(k,:)),'linewidth',1);
   hold on
end
 set(plot(NaN,NaN,'r'),'linewidth',1);

legend('equilibre','diffusion','dds','q_0 = 1/1','q_0 = 3/2','q_0 = 2/1','q_0 = 5/2','q_0 = 3/1','q_m_i_n < q_0')

plot(qdif(1,end),lidif(1),'ob',qeq(1,end),lieq(1),'og');
plot(qdif(end,end),lidif(end),'*b',qeq(end,end),lieq(end),'*g');

plot(contour_qali.limob(:,1),contour_qali.limob(:,2),'k',contour_qali.limoh(:,1),contour_qali.limoh(:,2),'k',contour_qali.limde(:,1),contour_qali.limde(:,2),'k:');

% recherche des q0 rationnelles
[rep,drep,rattab,dds] = isratq(qdif,qeq);
qm   = 0.5 .* (qdif(:,end) + qeq(:,end));
lim = 0.5 .* (lidif + lieq);
for k = 1:length(rattab)
    ind = find(rep(:,k));
    if ~isempty(ind)
       set(plot(qm(ind),lim(ind),sty(k,:)),'linewidth',1);
    end
    ind = find(drep(:,k));
    if ~isempty(ind)
         set(plot(qm(ind),lim(ind),strrep(sty(k,:),'m','r')),'linewidth',1);
    end
end

ind = find(dds);
if ~isempty(dds)
  set( plot(qm(ind),lim(ind),'+c'),'markersize',6,'linewidth',1);
end

% fonction de recherche de surface rationnelle
% rep -> rationnel au centre
% drep -> rationnel + profil creux
function [rep,drep,rattab,dds] = isratq(qdif,qeq)

   rattab = [1/1,3/2,2/1,5/2,3/1];
   rep = zeros(size(qdif,1),length(rattab));
   drep = rep;
   for k = 1:length(rattab)
      [rep(:,k),drep(:,k),dds] = isratq1(qdif,qeq,rattab(k));
   end

% pour une surface
function [ok,dok,dds] = isratq1(qdif,qeq,rat)

q   = max(qdif,qeq);   
q0  = mean(q(:,1:3),2);
tol = abs(qdif(:,1)-qeq(:,1));
ind = find(tol <1e-2);
if ~ isempty(ind)
   tol(ind) = 0.1;
end
qmin = min(q(:,5:end),[],2);

d  = abs(q0-rat);
ok = (d < tol);
dok = (d < tol)&(q0 > qmin);    

dds  =(q(:,5) < 1);
