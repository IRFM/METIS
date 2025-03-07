% calcul de eioniz en fonctiond de Te et ne pour le recyclage dans le divertor
% d'apres figure Stangeby 3.34 en attendant d'avoir mieux
% ne en m^-3 et te en eV
%[ne,te] = meshgrid(logspace(17,23),logspace(-1,4));
% e0 = z0eioniz_div(te,ne)
function e0 = z0eioniz_div(te,ne)

n_tab = [1e8,1e10,1e12,1e14,1e16]' .* 1e6;
t_tab = [1,10,100,1e4];
tab  = [2.5e3,40,27,32;
         2.5e3,40,27,32;
         1.5e3,35,27,30;
         250,23,20,27;
         18,15,15,15];
% ajout d'element pour etendre les bornes
t_tab = cat(2,0.1,t_tab,1e5);
tab   = cat(2,tab(:,1),tab,tab(:,end));
n_tab = cat(1,1,n_tab,1e24);
tab   = cat(1,tab(1,:),tab,tab(end,:));

% mise enforme
si = size(ne);
ne = ne(:);
te = te(:);

% interpolation sur ne
[n_log,t_log] = meshgrid(log10(n_tab),log10(t_tab));
e0   = exp(interp2(n_log,t_log,log(tab).',log10(max(1,ne)),log10(max(eps,te)),'linear',log(25)));
e0(~isfinite(e0)) = 25;
% mise en forme
e0 = reshape(e0,si);
