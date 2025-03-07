% script de plot des composante de epar
x           = param.gene.x;
t           = data.gene.temps;
ve          = ones(size(x));
vt          = ones(size(t));
dpsidt_rho  = - (vt * x) ./ (data.equi.rhomax * ve) .* (data.equi.drhomaxdt * ve) .* data.prof.psid1;
epar_stab   = - data.equi.F .* data.equi.r2i .* data.prof.dpsidt ./ (data.geo.b0 * ve);
epar_rho    = - data.equi.F .* data.equi.r2i .* dpsidt_rho ./ (data.geo.b0 * ve);
epar_3p     = - data.equi.F .* data.equi.r2i .* (data.prof.dpsidt3p + dpsidt_rho) ./ (data.geo.b0 * ve);

h   = findobj(0,'type','figure','tag','epar');
if isempty(h)
  h = figure('tag','epar','color',[1 1 1]);
end
clf

zplotprof(gca,t,x,data.prof.epar,'color','r');
zplotprof(gca,t,x,epar_stab,'color','b');
zplotprof(gca,t,x,epar_rho,'color','m','linestyle','o');
zplotprof(gca,t,x,epar_3p,'color','g','linestyle','-.');
legend('E//','Estab','Egeo','E3p');
xlabel('x');
ylabel('V/m');

