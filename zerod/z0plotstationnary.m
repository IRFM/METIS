% plot stationarity of a discharge
% ref : J_ss from raptor normalized
stat = sqrt(trapz(post.profil0d.xli,((post.profil0d.Raxe .* post.profil0d.epar - (post.profil0d.Raxe(:,end) .* post.profil0d.epar(:,end)) *  ...
        ones(size(post.profil0d.xli))) ./ post.profil0d.eta) .^ 2,2) ./ ...
        trapz(post.profil0d.xli,  (post.profil0d.jeff .* post.profil0d.Raxe) .^ 2,2)); 
%        trapz(post.profil0d.xli,  (post.profil0d.Raxe ./ post.profil0d.eta) .^ 2,2));
    
h = findobj(0,'type','figure','tag','z0plotstationnary');
if isempty(h)
       h=figure('tag','z0plotstationnary');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(4,1,1)
% choix des unites
if max(post.zerod.ip) > 1e5
  ipexp = 1e6;
  ipu   = 'MA';
else
  ipexp = 1e3;
  ipu   = 'kA';
end
if max(post.zerod.w) > 1e8
  wexp = 1e9;
  wu   = 'GJ';
elseif  max(post.zerod.w) > 1e5
  wexp = 1e6;
  wu   = 'MJ';
else
  wexp = 1e3;
  wu   = 'kJ';
end

plot(post.zerod.temps,post.zerod.ip ./ ipexp,post.zerod.temps,post.zerod.w ./ wexp)
title('Current diffusion stationarity')
legend(sprintf('I_p (%s)',ipu),sprintf('W_{total} (%s)',wu));
subplot(4,1,2)
plot(post.zerod.temps,post.zerod.vloop,post.zerod.temps,post.zerod.li,post.zerod.temps,zeros(size(post.zerod.temps)),'k:')
legend('V_{loop} (V)','l_i','Zero line');
subplot(4,1,3)
plot(post.zerod.temps,post.zerod.ialign,post.zerod.temps,post.zerod.ini ./ post.zerod.ipar,post.zerod.temps,ones(size(post.zerod.temps)),'k:');
legend('Current alignement (better -> one)','fraction on non inductive current','One line');
set(gca,'ylim',[0 2]);
subplot(4,1,4)
plot(post.profil0d.temps,stat);
z0loglin(gca);
xlabel('time (s)');
legend('Stationary criterium (stationary -> 0)');
joint_axes(h,4);
edition2