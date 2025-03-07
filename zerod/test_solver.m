% script de test de zerodevolution
% utilisation du mexfile si disponible	 
if isappdata(0,'MEXSOLVER_IN_METIS')
	mexsolver = getappdata(0,'MEXSOLVER_IN_METIS');
else	
	mexsolver =  [];
end
if isempty(mexsolver)
	repmex=which(strcat('mexpde1dsolver.',mexext));
	if ~isempty(repmex)
		mexsolver = 1;
	else
		mexsolver = 0;	
	end
	setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
end  

%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_TS_33850.mat');
%ref = load('/cgc_data/JA132999/zineb/zerod_data/METIS_JET67940_full_run.mat');
ref = load('/cgc_data/JA132999/zineb/zerod_data/metis_iter_itb_genetic_ok.mat');
z0dinput = ref.post.z0dinput;
ref = ref.post;
mexsolver = 1;
setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
[new.zerod,info,new.profil0d] = zerod(ref.z0dinput.option,ref.z0dinput.cons,ref.z0dinput.geo,ref.z0dinput.exp0d,1e-3,[]);
mexsolver = 0;
setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
[ref.zerod,info,ref.profil0d] = zerod(ref.z0dinput.option,ref.z0dinput.cons,ref.z0dinput.geo,ref.z0dinput.exp0d,1e-3,[]);



zplotstruct(new.zerod,ref.zerod,'0D','Metis')
zcompstruct(new.zerod,ref.zerod,'0D','Metis')

figure
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.qjli,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.qjli,'color','b','marker','o','linestyle','none')

figure
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.tep,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.tep,'color','m','marker','o','linestyle','none')
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.tip,'color','b')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.tip,'color','c','marker','o','linestyle','none')

figure
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.jli,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jli,'color','m','marker','o','linestyle','none')
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.jni,'color','b')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jni,'color','c','marker','o','linestyle','none')
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.jboot,'color','g')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.jboot,'color','k','marker','o','linestyle','none')

figure
zplotprof(gca,new.profil0d.temps,new.profil0d.xli,new.profil0d.eta,'color','r')
zplotprof(gca,ref.profil0d.temps,ref.profil0d.xli,ref.profil0d.eta,'color','b','marker','o','linestyle','none')


