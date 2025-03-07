% script de test de zerodevolution
%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_TS_33850.mat');
%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_JET67940_full_run.mat');
%ref = load('/cgc_data/JA132999/zineb/zerod_data/metis_iter_itb_genetic_ok.mat');
metis_load certification/metis/JET53521_itb_test.mat
z0dinput = post.z0dinput;
z0dinput.cons.flux = 0 .* z0dinput.cons.temps;
z0dinput.option.dwdt_method = 'none';

[zs,profil,z0dstruct] = zerodevolution([],z0dinput.option,z0dinput.cons.temps(1),z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
for k=2:length(z0dinput.cons.temps)
	[zs,profil,z0dstruct] = zerodevolution(z0dstruct,z0dinput.option,z0dinput.cons.temps(k),z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
	post.zerod = z0dstruct.zs;
	post.profil0d = z0dstruct.profil;
	post.z0dinput =z0dstruct.z0dinput;
	z0plotsc;
	drawnow
	figure(16);clf;
	subplot(3,1,1)
	plot(z0dstruct.zs.temps,z0dstruct.zs.taue,'b',z0dstruct.zs.temps,z0dstruct.zs.wth ./ z0dstruct.zs.pin,'r');
	ylabel('taue')
	subplot(3,1,2)
  	plot(z0dstruct.zs.temps,z0dstruct.zs.wth,'b',z0dstruct.zs.temps,cumtrapz(z0dstruct.zs.temps,z0dstruct.zs.dwthdt) + z0dstruct.zs.wth(1),'r' );
	ylabel('wth')
	subplot(3,1,3)
	plot(z0dstruct.zs.temps,z0dstruct.zs.dwthdt,'b',z0dstruct.zs.temps,z0dxdt(z0dstruct.zs.wth,z0dstruct.zs.temps),'r');
	ylabel('dwthdt')
	set(gca,'ylim',[min(z0dstruct.zs.dwthdt),max(z0dstruct.zs.dwthdt)]);
	drawnow

	figure(17);clf;
	subplot(3,1,1)
	plot(z0dstruct.zs.temps,z0dstruct.zs.taue,'b',post.zerod.temps,post.zerod.taue,'r');
	ylabel('taue')
	subplot(3,1,2)
  	plot(z0dstruct.zs.temps,z0dstruct.zs.wth,'b',post.zerod.temps,post.zerod.wth,'r' );
	ylabel('wth')
	subplot(3,1,3)
	plot(z0dstruct.zs.temps,z0dstruct.zs.dwthdt,'b',post.zerod.temps,post.zerod.dwthdt,'r');
	ylabel('dwthdt')
	drawnow

	figure(186)
	subplot(1,2,1)
	plot(profil.xli,profil.psi);
	hold on
	subplot(1,2,2)
	plot(profil.xli,profil.rmx);
	hold on
	
	figure(163);clf;
	subplot(3,1,1)
	plot(z0dstruct.zs.temps,z0dstruct.zs.pnbi,'b',z0dstruct.zs.temps,z0dstruct.zs.pnbi_th,'r');
	subplot(3,1,2)
	plot(z0dstruct.zs.temps,z0dstruct.zs.inbicd,'b',z0dstruct.zs.temps,z0dstruct.zs.ilh,'r');
	subplot(3,1,3)
	plot(z0dstruct.zs.temps,z0dstruct.zs.picrh,'b',z0dstruct.zs.temps,z0dstruct.zs.picrh_th,'r');
	
end

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


