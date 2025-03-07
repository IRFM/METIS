% script du  plot du shinethrough
h = findobj(0,'type','figure','tag','shine');
if isempty(h)
       h=figure('tag','shine');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(4,1,1)
plot(post.zerod.temps,post.zerod.nem./1e19,'b',post.zerod.temps,post.zerod.ne0./1e19,'r');
title(sprintf('Zerod : %s@%d/Shinethrough', ...
          post.z0dinput.machine,post.z0dinput.shot));
legend('<N_e>','n_e_0');
ylabel('10^{19} m^-^3');
subplot(4,1,2)
plot(post.zerod.temps,post.zerod.tem./1e3,'b',post.zerod.temps,post.zerod.te0./1e3,'r');
legend('<T_e>','T_e_0');
ylabel('keV');
subplot(4,1,3)
plot(post.zerod.temps,real(post.zerod.pnbi) ./1e6,'b');
if any(imag(post.zerod.pnbi))
  hold on 
  plot(post.zerod.temps,imag(post.zerod.pnbi) ./1e6,'r');
  legend('P_{NBI 1, deposition}','P_{NBI 2, deposition}');
else
  legend('P_{NBI, deposition}');
end
ylabel('MW');
subplot(4,1,4)
if isfield(post.zerod,'firstorb_nbi') && ~all(post.zerod.firstorb_nbi == 0) && all(isfinite(post.zerod.firstorb_nbi))
    semilogy(post.zerod.temps,1 - real(post.zerod.frnbi),'b',post.zerod.temps,real(post.zerod.firstorb_nbi),'c');
    if any(imag(post.zerod.frnbi))
      hold on 
      plot(post.zerod.temps,1 - imag(post.zerod.frnbi),'r',post.zerod.temps,imag(post.zerod.firstorb_nbi),'m');
      legend('1','2','Location','NorthWest','orientation','horizontal');
      legend('Shinethrough + first orbit losses 1','first orbit losses 1','Shinethrough + first orbit losses 2','first orbit losses 2');
    else
      legend('Shinethrough + first orbit losses','first orbit losses');
    end
else
    semilogy(post.zerod.temps,1 - real(post.zerod.frnbi),'b');
    if any(imag(post.zerod.frnbi))
      hold on 
      plot(post.zerod.temps,1 - imag(post.zerod.frnbi),'r');
      legend('1','2','Location','NorthWest','orientation','horizontal');
      legend('Shinethrough + first orbit losses 1','Shinethrough + first orbit losses 2');
    else
      legend('Shinethrough + first orbit losses');
    end
end
ylabel('fraction');
xlabel('time (s)');
z0loglin(gca);
edition2
joint_axes(h,4);