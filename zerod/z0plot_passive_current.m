% script pour tracer le courant dans les structures passives
if isfield(post.zerod,'eddy_current_eddy') && isfield(post.zerod,'flux_edge_cor')
    [taue,f,ieddy,nec,tref,prf,flux_bord_cor] = z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons,...
                      post.zerod,post.zerod.taue,post.zerod.eddy_current,post.zerod.flux_edge_cor);
else
    [taue,f,ieddy,nec,tref,prf,flux_bord_cor] = z0taue_burnthrough(post.z0dinput.option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.zerod.taue);
end
h = findobj(0,'type','figure','tag','z0plot_passive_current');
if isempty(h)
       h=figure('tag','z0plot_passive_current');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
k    = 2;
subplot(k,1,1)
plot(post.zerod.temps,post.zerod.ip ./ 1e6,post.zerod.temps,post.zerod.irun ./ 1e6, ...
     post.zerod.temps,ieddy ./ 1e6,post.zerod.temps,post.zerod.ini ./ 1e6, ...
     post.z0dinput.cons.temps,post.z0dinput.cons.ip ./ 1e6,'.');
%xlabel('time (s)')
ylabel('MA')
legend('I_p','I_{runaway}','I_{eddy}','I_{ni} (//B !)','I_{p,reference}');
z0loglin(gca);
title(sprintf('METIS : %s@%d / Eddy current ', ...
          post.z0dinput.machine,post.z0dinput.shot));
subplot(k,1,2)
plot(post.zerod.temps,post.zerod.vloop,post.zerod.temps,post.zerod.RR.*1e6);
legend('V_{loop} (V)','Resistor_{plasma} (\mu\Omega)');
xlabel('time (s)');
set(gca,'ylim',[0,10]);
joint_axes(h,k);
drawnow
edition2
