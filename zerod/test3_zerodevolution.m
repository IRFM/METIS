% script de test de zerodevolution
%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_TS_33850.mat');
%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_JET67940_full_run.mat');
%ref = load('/cgc_data/JA132999/zineb/zerod_data/metis_iter_itb_genetic_ok.mat');
%metis_load certification/metis/JET53521_itb_test.mat
%metis_load certification/metis/ITER_rampup_ECCD_expo_inte.mat
z0dinput = post.z0dinput;
z0dinput.cons.flux = 0 .* z0dinput.cons.temps;
%z0dinput.option.dwdt_method = 'none';
%z0dinput.option.dwdt_method = 'old';
z0dinput.option.dwdt_method = 'explicit';
%z0dinput.option.dwdt_method = 'mixed';
%z0dinput.option.dwdt_method = 'implicit';
%z0dinput.option.dwdt_method = 'v4.2';
%z0dinput.option.dwdt_method = 'implicitfilter';

z0dinput.option.dwdt_method

[zs,profil,z0dstruct] = zerodevolution([],z0dinput.option,z0dinput.cons.temps(1),z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
temps = zs.temps(end);
for k=2:(10*length(z0dinput.cons.temps))
	temps = temps + 0.1;
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,temps,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
	temps = zs.temps(end);
	
%  	figure(16);clf;
%  	subplot(3,1,1)
%  	plot(z0dstruct.zs.temps,z0dstruct.zs.taue,'b',z0dstruct.zs.temps,z0dstruct.zs.wth ./ z0dstruct.zs.pin,'r');
%  	ylabel('taue')
%  	subplot(3,1,2)
%    	plot(z0dstruct.zs.temps,z0dstruct.zs.wth,'b',z0dstruct.zs.temps,cumtrapz(z0dstruct.zs.temps,z0dstruct.zs.dwthdt) + z0dstruct.zs.wth(1),'r' );
%  	ylabel('wth')
%  	subplot(3,1,3)
%  	plot(z0dstruct.zs.temps,z0dstruct.zs.dwthdt,'b',z0dstruct.zs.temps,z0dxdt(z0dstruct.zs.wth,z0dstruct.zs.temps),'r');
%  	ylabel('dwthdt')
%  	set(gca,'ylim',[min(z0dstruct.zs.dwthdt),max(z0dstruct.zs.dwthdt)]);
%  	drawnow

	figure(17);clf;
	subplot(3,1,1)
	plot(z0dstruct.zs.temps,z0dstruct.zs.taue,'b',post.zerod.temps,post.zerod.taue,'r',zs.temps,zs.taue,'go' );
	ylabel('taue')
	subplot(3,1,2)
  	plot(z0dstruct.zs.temps,z0dstruct.zs.wth,'b',post.zerod.temps,post.zerod.wth,'r',zs.temps,zs.wth,'go' );
	ylabel('wth')
	subplot(3,1,3)
	plot(z0dstruct.zs.temps,z0dstruct.zs.dwthdt,'b',post.zerod.temps,post.zerod.dwthdt,'r',zs.temps,zs.dwthdt,'go');
	ylabel('dwthdt')
	drawnow

%  	figure(186)
%  	subplot(1,2,1)
%  	plot(profil.xli,profil.psi);
%  	hold on
%  	subplot(1,2,2)
%  	plot(profil.xli,profil.rmx);
%  	hold on
	
	figure(163);clf;
	subplot(4,1,1)
	plot(z0dstruct.zs.temps,z0dstruct.zs.ploss,'b',post.zerod.temps,post.zerod.ploss,'r');
	subplot(4,1,2)
	plot(z0dstruct.zs.temps,z0dstruct.zs.ini,'b',post.zerod.temps,post.zerod.ini,'r');
	subplot(4,1,3)
	plot(z0dstruct.zs.temps,z0dstruct.zs.te0,'b',post.zerod.temps,post.zerod.te0,'r');
	subplot(4,1,4)
	plot(z0dstruct.zs.temps,z0dstruct.zs.pfus_nbi,'b',post.zerod.temps,post.zerod.pfus_nbi,'r');
	
end

if 0
% supression des points supplementaires
%mise enplace du dernier temps
noms = fieldnames(post.zerod);
for l=1:length(noms)
	nomc = noms{l};
	if size(post.zerod.(nomc),1) > 1
		post.zerod.(nomc) = post.zerod.(nomc)(7:end);
	else
		post.zerod.(nomc) = post.zerod.(nomc);		
	end
end
noms = fieldnames(post.profil0d );
for l=1:length(noms)
	nomc = noms{l};
	if size(post.profil0d.(nomc),1) > 1
		post.profil0d.(nomc) = post.profil0d.(nomc)(7:end,:);
	else
		post.profil0d.(nomc) = post.profil0d.(nomc);		
	end
end

[ref.zerod,info,ref.profil0d] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d,1e-3,[]);


zplotstruct(post.zerod,ref.zerod,'0D','Metis')

figure
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.qjli,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.qjli,'color','b')

figure
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.tep,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.tep,'color','m')
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.tip,'color','b')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.tip,'color','c')

figure
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jli,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jli,'color','m')
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jni,'color','b')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jni,'color','c')
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.jboot,'color','g')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jboot,'color','k')

figure
zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.eta,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.eta,'color','b')

end
