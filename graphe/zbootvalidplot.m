% plot de la validite du clacul par la methode perturbative du boostrap 
% qui est la methode de NClass, DKE et O. Sauter
% 
% ref : small-drift expansion of the dke, Joan Decker 2006, private communication
function zbootvalidplot


phys   = evalin('base','param.phys');
gene   = evalin('base','param.gene');
machine = evalin('base','param.from.machine');
num     = evalin('base','param.from.shot.num');
compo  = evalin('base','param.compo');
te     = evalin('base','data.prof.te');
ti     = evalin('base','data.prof.ti');
ne     = evalin('base','data.prof.ne');
ni     = evalin('base','data.prof.ni');
q      = evalin('base','data.prof.q');
a0     = evalin('base','data.geo.a');
raxe   = evalin('base','data.equi.raxe');
b2i     = evalin('base','data.equi.b2i');
rm     = evalin('base','data.equi.rhomax');
t      = evalin('base','data.gene.temps');

% vecteur
vt = ones(size(t));
ve = ones(size(gene.x));

% rapport d'aspect
epsi = max(eps,(a0 * gene.x) ./ raxe);

% calcul pour les ions
% rayon de larmor  (Wesson)
rhoi = 4.57e-3 .* sqrt(min(compo.a)) .* sqrt(ti ./ 1e3) .* sqrt(b2i);
ilti = abs(pdederive(gene.x,ti,0,2,2,1) ./ (rm * ve) ./ max(13.6,ti));
ilni = abs(pdederive(gene.x,ni,0,2,2,1) ./ (rm * ve) ./ max(1e13,ni));
ci   = 4 .* q ./ epsi .* rhoi .* max(ilti,ilni);
ci(:,1) = 0;


% calcul pour les electrons
% rayon de larmor  (Wesson)
rhoe = 1.07e-4  .* sqrt(te ./ 1e3) .* sqrt(b2i);
ilte = abs(pdederive(gene.x,te,0,2,2,1) ./ (rm * ve) ./ max(13.6,te));
ilne = abs(pdederive(gene.x,ne,0,2,2,1) ./ (rm * ve) ./ max(1e13,ne));
ce   = 4 .* q ./ epsi .* rhoe .* max(ilte,ilne);
ce(:,1) = 0;

% pour visualisation
sigma   =  10;
prob_e  = tanh((ce - 1) .* sigma);
prob_i  = tanh((ci - 1) .* sigma);

% mise en forme
t = t * ve;
x = vt * gene.x;


% borne en temps
if strcmp(gene.filetype,'resultat');
   indok = gene.kmin:gene.k;
   t       = t(indok,:);
   x       = x(indok,:);
   prob_e  = prob_e(indok,:);
   prob_i  = prob_i(indok,:);
   ce      = ce(indok,:);
   ci      = ci(indok,:);
end

h = findobj(0,'type','figure','tag','bootvalid');
if isempty(h)
       h=figure('tag','bootvalid');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',6)


colormap(flipud(hot));

subplot(2,2,1)
contourf(t,x,ce,21); 
title(sprintf('%s #%d',machine,fix(num)))
ylabel('x');
axes(colorbar);
title('electron bootstrap')

subplot(2,2,2)
contourf(t,x,ci,21); 
ylabel('x');
axes(colorbar);
title('ion bootstrap')

subplot(2,2,3)
contourf(t,x,prob_e); 
set(gca,'clim',[-1,1]);
title('Prob(Ce > 1)')
xlabel('temps (s)');
ylabel('x');
axes(colorbar);
title('valid if <0')

subplot(2,2,4)
contourf(t,x,prob_i); 
set(gca,'clim',[-1,1]);
title('Prob(Ci > 1)')
xlabel('temps (s)');
ylabel('x');
axes(colorbar);
title('valid if <0')

