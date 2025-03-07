% plot des barières 
% cette focntion trace la probabilte d'avoir une ITB sur le canal eletronique et le canal ionique
% ref : G. Tresset et al ,Nuclear Fusion, vol 42, pp 520-526, 2002
function zitbplot

rho_itb  = 1.4e-2;

phys   = evalin('base','param.phys');
gene   = evalin('base','param.gene');
machine = evalin('base','param.from.machine');
num     = evalin('base','param.from.shot.num');
compo  = evalin('base','param.compo');
te     = evalin('base','data.prof.te');
ti     = evalin('base','data.prof.ti');
b2     = evalin('base','data.equi.b2');
b0     = evalin('base','data.geo.b0');
ind    = find(~isfinite(b2));
if ~isempty(ind)
   b2(ind) = mean(b0(isfinite(b0)) .^2);
end
a0     = evalin('base','data.equi.rhomax');
t     = evalin('base','data.gene.temps');


cs     = sqrt(phys.e .* (5/3 .* ti + compo.z(1) .* te) ./ compo.a(1) ./ phys.mp);
wci    = compo.z(1) .* phys.e .* sqrt(b2) ./ phys.mp ./ compo.a(1);
rhos   = cs ./ wci;

r     = a0 * gene.x;
t     = t  * ones(1,size(r,2));
x     = ones(size(t,1),1) * gene.x;

lte   = max(te,13.6) ./ (abs(pdederive(r,te,0,2,2,1))+ 13.6 ./ 3 ./ max(a0(isfinite(a0))));
lti   = max(ti,13.6) ./ (abs(pdederive(r,ti,0,2,2,1))+ 13.6 ./ 3 ./ max(a0(isfinite(a0))));

rhos_te  = rhos ./ lte;
%rhos_te(rhos_te <0) = 0;
%rhos_te(pdederive(r,te,0,2,2,1)> 0) = 0;

rhos_ti  = rhos ./ lti;
%rhos_ti(rhos_ti <0) = 0;
%rhos_ti(pdederive(r,ti,0,2,2,1)> 0) = 0;

sigma    = 0.05;

warning off
tez1 = pdederive(te,rhos_te,2,2,2,1) .^ 2;
ind  = find(~isfinite(tez1));
if ~isempty(ind)
   tez1(ind) = 0;
end 
tez2 = pdederive(ti,rhos_te,2,2,2,1) .^ 2;
ind  = find(~isfinite(tez2));
if ~isempty(ind)
   tez2(ind) = 0;
end 
tez3 = pdederive(te,rhos_te,2,2,1,1) .^ 2;
ind  = find(~isfinite(tez3));
if ~isempty(ind)
   tez3(ind) = 0;
end 
tez4 = pdederive(ti,rhos_te,2,2,1,1) .^ 2;
ind  = find(~isfinite(tez4));
if ~isempty(ind)
   tez4(ind) = 0;
end 

tiz1 = pdederive(te,rhos_ti,2,2,2,1) .^ 2;
ind  = find(~isfinite(tiz1));
if ~isempty(ind)
   tiz1(ind) = 0;
end 
tiz2 = pdederive(ti,rhos_ti,2,2,2,1) .^ 2;
ind  = find(~isfinite(tiz2));
if ~isempty(ind)
   tiz2(ind) = 0;
end 
tiz3 = pdederive(te,rhos_ti,2,2,1,1) .^ 2;
ind  = find(~isfinite(tiz3));
if ~isempty(ind)
   tiz3(ind) = 0;
end 
tiz4 = pdederive(ti,rhos_ti,2,2,1,1) .^ 2;
ind  = find(~isfinite(tiz4));
if ~isempty(ind)
   tiz4(ind) = 0;
end 


sigma_te = sigma  .* max(te,13.6)  .* sqrt(tez1 + tez2 + tez3 + tez4);
sigma_ti = sigma  .* max(ti,13.6)  .* sqrt(tiz1 + tiz2 + tiz3 + tiz4);

                    
warning on

indnok = find(~isfinite(sigma_te));
indok = find(isfinite(sigma_te));
if ~isempty(indnok)
   sigma_te(indnok) = mean(sigma_te(indok));
end
indnok = find(sigma_te == 0);
if ~isempty(indnok)
   sigma_te(indnok) = mean(sigma_te(indok));
end


indnok = find(~isfinite(sigma_ti));
indok = find(isfinite(sigma_ti));
if ~isempty(indnok)
   sigma_ti(indnok) = mean(sigma_ti(indok));
end
indnok = find(sigma_ti == 0);
if ~isempty(indnok)
   sigma_ti(indnok) = mean(sigma_ti(indok));
end


prob_te  = 0.5 .* (1+erf((rhos_te - rho_itb)./ sqrt(2) ./ sigma_te));
prob_ti  = 0.5 .* (1+erf((rhos_ti - rho_itb)./ sqrt(2) ./ sigma_ti));

% borne en temps
if strcmp(gene.filetype,'resultat');
   indok = gene.kmin:gene.k;
   t       = t(indok,:);
   r       = r(indok,:);
   x       = x(indok,:);
   prob_te = prob_te(indok,:);
   prob_ti = prob_ti(indok,:);
   rhos_te = rhos_te(indok,:);
   rhos_ti = rhos_ti(indok,:);
end


h = findobj(0,'type','figure','tag','itb');
if isempty(h)
       h=figure('tag','itb');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',6)


colormap(flipud(hot));

subplot(2,2,1)
u = min(0.1,rhos_te);
contourf(t,x,u); 
title(sprintf('%s #%d',machine,fix(num)))
ylabel('x');
axes(colorbar);
title('rho_sT_e')

subplot(2,2,2)
u = min(0.1,rhos_ti);
contourf(t,x,u); 
ylabel('x');
axes(colorbar);
title('rho_sT_i')

subplot(2,2,3)
u = max(49,100.* prob_te);
contourf(t,x,u); 
set(gca,'clim',[50 100]);
title('Prob(rho_sT_e > rho_I_T_B)')
xlabel('temps (s)');
ylabel('x');
axes(colorbar);
title('%')

subplot(2,2,4)
u = max(49,100 .* prob_ti);
contourf(t,x,u); 
set(gca,'clim',[50 100]);
title('Prob(rho_sT_i > rho_I_T_B)')
xlabel('temps (s)');
ylabel('x');
axes(colorbar);
title('%')

