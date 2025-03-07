% E. Nardon, 25/01/11
% Calcule le flux consomme en fonction du temps.
% On utilise l'approximation Lsol (self du solenoide) = Msol,p (mutuelle entre le solenoide et le plasma) qui serait valable si le solenoide etait infini.

mu0 = 4*pi*10^(-7);

t    = post.zerod.temps;
geo  = post.z0dinput.geo;
R    = geo.R + post.zerod.d0
% Self du plasma (formule approchee pour un profil de courant plat) :
L_plasma = mu0 * R .* (log(8 *  R ./ ((post.zerod.sp/pi).^0.5))-2);

L_plasma_dot = z0dxdt(L_plasma,t);
ip_dot = z0dxdt(post.zerod.ip,t);

flux_cons_dot = post.zerod.vloop + 0.5 * L_plasma_dot.*post.zerod.ip + L_plasma.*ip_dot;
flux_cons = cumtrapz(t ,flux_cons_dot);

figure
hold on
set(gcf,'color','w')
plot(t,flux_cons,'linewidth',2)
grid
font_s = 18;
set(gca,'fontsize',font_s)
xlabel('Time (s)','fontsize',font_s)
ylabel('Consumed flux (Wb)','fontsize',font_s)