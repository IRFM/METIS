% creation des figures
tag=sprintf('l2h_metis');
h = findobj(0,'type','figure','tag',tag);
if isempty(h)
	h=figure('tag',tag);
else
	figure(h);
end   
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultaxeslinewidth',3,'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',12,'Position',[400 40 1024 960],...
'PaperPositionMode','auto')
cc  = get(gca,'colororder');
ml{1} = 'o';
ml{2} = '*';
ml{3} = 's';
ml{4} = 's';
ml{5} = 'p';
ml{6} = 'h';
ml{7} = '^';
ml{8} = '>';
ml{9} = '<';
ml{10} = 'v';
ml{11} = '+';

switch post.z0dinput.option.l2hscaling
case 0
	scl2h = 'LH99(1)';
case 1
	scl2h = 'LH02noZeff (Ryter,PPCF,2002,A415)';
case 2
	scl2h = 'LH02Zeff (Takizuka,PPCF, 2004)';
case 3
	scl2h = 'Y R Martin, Journal of Physics 2008 Conference series 123 p 012033';
case 4
	sc2lh = 'NLM-7 in A. Murari et al, NF 52 (2012) p 0630016-'; 
case 5 
	sc2lh = 'NLM-11 in A. Murari et al, NF 52 (2012) p 0630016-'; 
case 6
	scl2h = 'J_theta = 0';
case 10
	scl2h = 'A. Murari et al, NF 53 (2013) p 043001-043013, formula 6';
case 28
	scl2h = 'Low density case - Ryter et al, NF 54 (2014) 083003, equation 4';
case 30
	scl2h = 'Fit of metalic tokamaks database for horizontal targets (E . Delabie et al, 2025 ?)';
case 31
	scl2h = ' Fit of metalic tokamaks database for vertical targets and corner configuration (E . Delabie et al, 2025 ?);';
otherwise
	if post.z0dinput.option.l2hscaling < 0
		scl2h = 'Rotation  (D. Testa Nucl. Fusion 46 (2006) 562-579)';
	else
		error('unknown scaling law for L2H transition');
	end
end

if post.z0dinput.option.fpped > 0
     if post.z0dinput.option.hmore_pped == 1
 	ppedmax = post.zerod.ppedmax .* abs(post.z0dinput.option.fpped) .* post.z0dinput.cons.hmore;
     elseif post.z0dinput.option.hmore_pped == 2
 	ppedmax    = post.zerod.ppedmax .* abs(post.z0dinput.option.fpped) .* min(post.z0dinput.cons.hmore,1);
     else
	ppedmax    = post.zerod.ppedmax .* abs(post.z0dinput.option.fpped);
     end
else
    ppedmax = post.zerod.ppedmax;
end

subplot(5,1,1)
plot(post.profil0d.temps,post.profil0d.fdia(:,end),post.zerod.temps,post.zerod.nbar ./ 1e19,post.zerod.temps,post.zerod.ip ./ 1e6);
legend('R B_0 (T m)','n_{bar} (10^{19} m^{-3})','Ip (MA)');
title(sprintf('Zerod : %s@%d/L -> H transition',post.z0dinput.machine,post.z0dinput.shot));
subplot(5,1,2)
plot(post.zerod.temps,post.zerod.zeff,post.zerod.temps,post.zerod.meff);
legend('Z_{eff}','m_{eff}');
set(gca,'ylim',[1,Inf])
subplot(5,1,3)
plot(post.zerod.temps,post.zerod.sext,post.zerod.temps,post.zerod.vp);
legend('S_{ext}','Volume');
ylabel('m^2');
subplot(5,1,4)
if post.z0dinput.option.l2hscaling >= 0
	plot(post.zerod.temps,post.z0dinput.option.l2hmul + post.zerod.plossl2h ./ 1e6,post.zerod.temps,post.zerod.plhthr ./ 1e6, ...
        post.zerod.temps,post.zerod.pin ./ 1e6, post.zerod.temps, ...
        (post.z0dinput.cons.picrh + real(post.z0dinput.cons.pnbi) + imag(post.z0dinput.cons.pnbi)  +  ....
         post.z0dinput.cons.pecrh + post.z0dinput.cons.plh + post.zerod.pohm) ./ 1e6);
	legend('P_{threshold}','P_{loss,LCFS}','P_{in}','P_{add,external}');
	ylabel('MW');
else
		meff     = interp1(post.zerod.temps,post.zerod.meff,post.profil0d.temps,'linear','extrap');
		gitg     = sqrt(1.602176462e-19 .* (post.profil0d.tip + post.profil0d.zeff .* post.profil0d.tep)./1.6726485e-27 ./  ...
			(meff * ones(size(post.profil0d.xli)))) ./ ...
			(post.profil0d.rmx(:,end) * ones(size(post.profil0d.xli)));
		semilogy(post.profil0d.temps,mean(post.profil0d.web(:,end-1:end-1),2),'r',post.profil0d.temps,abs(post.z0dinput.option.l2hscaling) .* ...
                         mean(gitg(:,end-1:end-1),2),'b');
		legend('\Gamma _{ExB}','\Gamma _{ITG}');		 
end

subplot(5,1,5)
%Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
% use fit of pedestal database augmented with prediction from MHD stability studies for ITER and DEMO1
switch post.z0dinput.option.ploss_exp
case {'max_power','max(pel)+max(pion)'}
	p_scaling_pped =  post.zerod.ploss;
otherwise
	p_scaling_pped =  post.zerod.pin;
end

% Pped_{sc,minimal} = 3.06928 * ip ^ 0.960252 * Bt ^ 0.63543
% Pped_{sc} = 4.53138 * (d + 0.0034) ^ 0.435509 * ptot ^ 0.121836 * iob ^ 1.51649
iob = (2/5) .* post.z0dinput.geo.R .* post.zerod.ip ./ 1e6  .* post.z0dinput.geo.K .^ 2 ./  post.z0dinput.geo.a .^ 2 ./ (1 + post.z0dinput.geo.K .^ 2);
pped_sc   = 4.53138 .* (post.z0dinput.geo.d + 0.0034) .^ 0.435509 .* (p_scaling_pped./ 1e6) .^ 0.121836 .* iob .^ 1.51649; 
% Pped_{sc} = 3.53664 * fgr ^ 0.0531514 * ip ^ 0.604636 * bt ^ 1.00281 * K ^ 0.566784 * (d + 0.187) ^ 0.631146 * li ^ -0.647307
plot(post.zerod.temps,post.zerod.pped ./1e3,post.zerod.temps,ppedmax ./1e3,post.zerod.temps,pped_sc);
legend('P_{ped}','P_{ped, maximum allowed}','P_{ped,scaling,extended}');
xlabel('Time (s)');
ylabel('kPa')

joint_axes(h,5);





