function test_psi_lcfs_evol(post)
% test formule Ejima

if nargin == 0
    post = evalin('base','post;');
end

mu0       = 4 * pi * 1e-7;
res       = post.zerod.pohm ./ post.zerod.ip .^ 2;
psi_res   = cumtrapz(post.zerod.temps,res .* post.zerod.ip);
L_ind_vol = mu0 * post.z0dinput.geo.R .* post.zerod.li / 2;
dipdt     = z0dxdt(post.zerod.ip,post.zerod.temps);
psi_ind   = 0.5 * L_ind_vol .* post.zerod.ip + 0.5 .* ...
            cumtrapz(post.zerod.temps,L_ind_vol .* dipdt);
        
h = findobj(0,'type','figure','tag','flux_ejima');
if isempty(h)
       h=figure('tag','flux_ejima');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
    'defaultlinelinewidth',1,'color',[1 1 1])

plot(post.zerod.temps,psi_ind(1) - post.zerod.edgeflux,'k.', ...
    post.zerod.temps,psi_res,'r',post.zerod.temps,psi_ind,'b', ...
    post.zerod.temps,psi_res + psi_ind ,'g', ...
    post.profil0d.temps,2*pi*(post.profil0d.psi(1,1) - post.profil0d.psi(:,1)),'m', ...
    post.profil0d.temps,2*pi *(post.profil0d.psi(:,1) - post.profil0d.psi(:,end)),'c');

xlabel('time (s)');
ylabel('Wb');
legend('Psi_{LCFS} - Psi_{offset}','Psi_{res}','Psi_{ind}','Psi_{res+ind}', ...
    'Psi_{res,axis}','Psi_{axis} - Psi_{LCFS}');


title(sprintf('Ejima flux consumption computation : %s@%d', ...
    post.z0dinput.machine,post.z0dinput.shot));
edition2