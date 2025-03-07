function z0plotnbi

%  try
%   if nargin('mdsdisconnect') > 0
%        mdsdisconnect( 'mdsplus.jet.efda.org');
%   else 
%        mdsdisconnect;
%   end
%  end


post = evalin('base','post');

switch post.z0dinput.machine
case 'JET'
    % script de verification du depot NBI pour JET	
    if exist('jdata')
	    if jdata.shot  ~= post.z0dinput.shot
		    clear jdata id_transp
	    end
    end
    if ~exist('id_transp')  && ~ isappdata(0,'TRANSP_USER')
        %id_transp = input('id transp ? (I03, ...)','s');
        prompt={'id TRANSP ?'};
        name='Choose a TRANSP run identifier';
        numlines=1;
        defaultanswer={'I03'};
        answer = {''};
        while(isempty(answer{1}))
            answer=inputdlg(prompt,name,numlines,defaultanswer);
        end
        id_transp = answer{1};
        id_transp = num2str(100 .* (abs(lower(id_transp(1))) - abs('a') + 1) + str2num(id_transp(2:end)));
        if length(id_transp) < 4
            id_transp = strcat('0',id_transp);
        end
    else
        id_transp = 'I03';
        id_transp = num2str(100 .* (abs(lower(id_transp(1))) - abs('a') + 1) + str2num(id_transp(2:end)));
        if length(id_transp) < 4
            id_transp = strcat('0',id_transp);
        end        
    end
    if ~exist('jdata')
	    
	    liste          = {};
	    liste{end+1}   = 'ppf/@shot/NBIP/GI';    
	    liste{end+1}   = 'ppf/@shot/NBIP/GE';    
	    liste{end+1}   = 'ppf/@shot/NBIP/GTT';    
	    liste{end+1}   = 'ppf/@shot/NBIP/JBDC';  
	    liste{end+1}   = 'ppf/@shot/NBIP/JBDT';
	    liste{end+1}   = 'ppf/@shot/NBIP/TORP';
	    liste{end+1}   = 'ppf/@shot/TRA0/CB';  
	    liste{end+1}   = 'ppf/@shot/TRA0/CBS';  
	    liste{end+1}   = 'ppf/@shot/TRA0/QBE';  
	    liste{end+1}   = 'ppf/@shot/TRA0/QBI';  
	    liste{end+1}   = 'ppf/@shot/TRA0/TI';  
	    liste{end+1}   = 'ppf/@shot/EFIT/Q';  
	    liste{end+1}   = 'ppf/@shot/CXSM/WTOT';    % energie thermique
	    liste{end+1}   = 'ppf/@shot/CXSM/ANGF';    % rotation
	    liste{end+1}   = 'ppf/@shot/CXSM/TIMX';    % Ti0
	    liste{end+1}   = 'ppf/@shot/CXSM/TI';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/CXSM/TICR';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/CXSM/TIRH';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/CXGM/WTOT';    % energie thermique
	    liste{end+1}   = 'ppf/@shot/CXGM/ANGF';    % rotation
	    liste{end+1}   = 'ppf/@shot/CXGM/TIMX';    % Ti0
	    liste{end+1}   = 'ppf/@shot/CXGM/TI';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/CXGM/TICR';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/CXGM/TIRH';      % profil de Ti
	    %liste{end+1}   = 'ppf/@shot/CXFM/WTOT';    % energie thermique
	    liste{end+1}   = 'ppf/@shot/CXFM/ANGF';    % rotation
	    %liste{end+1}   = 'ppf/@shot/CXFM/TIMX';    % Ti0
	    liste{end+1}   = 'ppf/@shot/CXFM/TI';      % profil de Ti
	    %liste{end+1}   = 'ppf/@shot/CXFM/TICR';      % profil de Ti
	    %liste{end+1}   = 'ppf/@shot/CXFM/TIRH';      % profil de Ti
	    liste{end+1}   = 'ppf/@shot/EQUI/Q';  
	    liste{end+1}   = 'ppf/@shot/XCS/TE';       % Te bragg
	    liste{end+1}   = 'ppf/@shot/XCS/TI';       % Ti bragg
	    
	    jdata           = cgcgetjet(post.z0dinput.shot,liste,'','');
	    NBIP            = jdata.ppf.NBIP;
	    TRA0            = jdata.ppf.TRA0;
	    EFIT            = jdata.ppf.EFIT;
	    EQUI            = jdata.ppf.EQUI;
	    CX              = jdata.ppf.CXFM;
	    CXS             = jdata.ppf.CXSM;
	    CXG             = jdata.ppf.CXGM;
	    BRAG            = jdata.ppf.XCS;
	    jdata.shot      = post.z0dinput.shot;
        
        if (isempty(TRA0.QBE.data) || isstr(TRA0.QBE.data)) && ...
                isappdata(0,'TRANSP_USER') && isappdata(0,'TRANSP_RUN')
                seq = getappdata(0,'TRANSP_RUN');
                [qui,arbre] = strtok(getappdata(0,'TRANSP_USER'),'/');
                arbre = strrep(arbre,'/','');
                liste          = {};
                liste{end+1} = sprintf('ppf/@shot/%s/CB?uid=%s+seq=%d',arbre,qui,seq);
                liste{end+1} = sprintf('ppf/@shot/%s/QBE?uid=%s+seq=%d',arbre,qui,seq);
                liste{end+1} = sprintf('ppf/@shot/%s/QBI?uid=%s+seq=%d',arbre,qui,seq);
                % data/pulse/99971/ppf/signal/mporad/trau/zef
                jdata2          = cgcgetjet(post.z0dinput.shot,liste,'','');
                TRA0            = jdata2.ppf.(arbre);
       end
	    
 	    if isempty(TRA0.QBE.data) || isstr(TRA0.QBE.data)
		    liste          = {};
		    liste{end+1} = ['ppf/@shot/TRA0/CB?uid=','TRANPPF'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBE?uid=','TRANPPF'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBI?uid=','TRANPPF'];
		    jdata2          = cgcgetjet(post.z0dinput.shot,liste,'','');
		    TRA0            = jdata2.ppf.TRA0;
	    end        
	    if isempty(TRA0.QBE.data) || isstr(TRA0.QBE.data)
		    liste          = {};
		    liste{end+1} = ['ppf/@shot/TRA0/CB?uid=','TRANPPF','+seq=','0'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBE?uid=','TRANPPF','+seq=','0'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBI?uid=','TRANPPF','+seq=','0'];
		    jdata2          = cgcgetjet(post.z0dinput.shot,liste,'','');
		    TRA0            = jdata2.ppf.TRA0;
	    end 
	    if isempty(TRA0.QBE.data) || isstr(TRA0.QBE.data)
		    liste          = {};
		    liste{end+1} = ['ppf/@shot/TRA0/CB?uid=','XII','+seq=','0'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBE?uid=','XII','+seq=','0'];
		    liste{end+1} = ['ppf/@shot/TRA0/QBI?uid=','XII','+seq=','0'];
		    jdata2       = cgcgetjet(post.z0dinput.shot,liste,'','');
		    TRA0        = jdata2.ppf.TRA0;
        end	
        if isempty(CX.ANGF.data) || isstr(CX.ANGF.data)
            liste          = {};
            liste{end+1}   = 'ppf/@shot/CXG6/WTOT';    % energie thermique
            liste{end+1}   = 'ppf/@shot/CXG6/ANGF';    % rotation
            liste{end+1}   = 'ppf/@shot/CXG6/TIMX';    % Ti0
            liste{end+1}   = 'ppf/@shot/CXG6/TI';      % profil de Ti
            liste{end+1}   = 'ppf/@shot/CXG6/TICR';      % profil de Ti
            liste{end+1}   = 'ppf/@shot/CXG6/TIRH';      % profil de Ti
            jdata37         = cgcgetjet(post.z0dinput.shot,liste,'','');
            if ~isempty(jdata37.ppf.CXG6.ANGF.data) & ~isstr(jdata37.ppf.CXG6.ANGF.data)
                CX             = jdata37.ppf.CXG6;
            end
        end
	    if isempty(CX.ANGF.data) || isstr(CX.ANGF.data)
		    liste          = {};
		    liste{end+1}   = 'ppf/@shot/CXFM/ANGF?uid=cxsbatch+seq=0';    % rotation
		    liste{end+1}   = 'ppf/@shot/CXFM/TI?uid=cxsbatch+seq=0';      % profil de Ti
		    jdata3         = cgcgetjet(post.z0dinput.shot,liste,'','');
		    if ~isempty(jdata3.ppf.CXFM.ANGF.data) & ~isstr(jdata3.ppf.CXFM.ANGF.data)
			    CX             = jdata3.ppf.CXFM;
		    end
	    end	
	    if isempty(CXS.ANGF.data) || isstr(CXS.ANGF.data)
		    liste          = {};
		    liste{end+1}   = 'ppf/@shot/CXSM/WTOT?uid=cxsbatch+seq=0';    % energie thermique
		    liste{end+1}   = 'ppf/@shot/CXSM/ANGF?uid=cxsbatch+seq=0';    % rotation
		    liste{end+1}   = 'ppf/@shot/CXSM/TIMX?uid=cxsbatch+seq=0';    % Ti0
		    liste{end+1}   = 'ppf/@shot/CXSM/TI?uid=cxsbatch+seq=0';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CXSM/TICR?uid=cxsbatch+seq=0';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CXSM/TIRH?uid=cxsbatch+seq=0';      % profil de Ti
		    jdata4         = cgcgetjet(post.z0dinput.shot,liste,'','');
		    if ~isempty(jdata4.ppf.CXSM.ANGF.data) & ~isstr(jdata4.ppf.CXSM.ANGF.data)
			    CXS             = jdata4.ppf.CXSM;
		    end
	    end	
	    if isempty(CXG.ANGF.data) || isstr(CXG.ANGF.data)
		    liste          = {};
		    liste{end+1}   = 'ppf/@shot/CXGM/WTOT?uid=cxsbatch+seq=0';    % energie thermique
		    liste{end+1}   = 'ppf/@shot/CXGM/ANGF?uid=cxsbatch+seq=0';    % rotation
		    liste{end+1}   = 'ppf/@shot/CXGM/TIMX?uid=cxsbatch+seq=0';    % Ti0
		    liste{end+1}   = 'ppf/@shot/CXGM/TI?uid=cxsbatch+seq=0';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CXGM/TICR?uid=cxsbatch+seq=0';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CXGM/TIRH?uid=cxsbatch+seq=0';      % profil de Ti
		    jdata4         = cgcgetjet(post.z0dinput.shot,liste,'','');
		    if ~isempty(jdata4.ppf.CXGM.ANGF.data) & ~isstr(jdata4.ppf.CXGM.ANGF.data)
			    CXG             = jdata4.ppf.CXGM;
		    end
	    end	
	    if isempty(CX.ANGF.data) || isstr(CX.ANGF.data)
		    liste          = {};
		    liste{end+1}   = 'ppf/@shot/CX/WTOT';    % energie thermique
		    liste{end+1}   = 'ppf/@shot/CX/ANGF';    % rotation
		    liste{end+1}   = 'ppf/@shot/CX/TIMX';    % Ti0
		    liste{end+1}   = 'ppf/@shot/CX/TI';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CX/TICR';      % profil de Ti
		    liste{end+1}   = 'ppf/@shot/CX/TIRH';      % profil de Ti
		    jdata31         = cgcgetjet(post.z0dinput.shot,liste,'','');
		    if ~isempty(jdata31.ppf.CX.ANGF.data) & ~isstr(jdata31.ppf.CX.ANGF.data)
			    CX             = jdata31.ppf.CX;
		    end
	    end	

	    % essai de lecture de chain2
	    liste          = {};
	    liste{end+1}   = 'ppf/@shot/NBP2/GI';    
	    liste{end+1}   = 'ppf/@shot/NBP2/GE';    
	    liste{end+1}   = 'ppf/@shot/NBP2/GTT';    
	    liste{end+1}   = 'ppf/@shot/NBP2/JBDC';  
	    liste{end+1}   = 'ppf/@shot/NBP2/JBDT';
	    liste{end+1}   = 'ppf/@shot/NBP2/TORP';
	    jdata5         = cgcgetjet(post.z0dinput.shot,liste,'','');
	    if ~isempty(jdata5.ppf.NBP2.TORP.data) & ~isstr(jdata5.ppf.NBP2.TORP.data)
			    NBIP.TORP            = jdata5.ppf.NBP2.TORP;
	    end
	    if ~isempty(jdata5.ppf.NBP2.JBDT.data) & ~isstr(jdata5.ppf.NBP2.JBDT.data)
			    NBIP.JBDT            = jdata5.ppf.NBP2.JBDT;
	    end
	    if ~isempty(jdata5.ppf.NBP2.JBDC.data) && ~isstr(jdata5.ppf.NBP2.JBDC.data)
			    NBIP.JBDC            = jdata5.ppf.NBP2.JBDC;
	    end
	    if ~isempty(jdata5.ppf.NBP2.GTT.data) && ~isstr(jdata5.ppf.NBP2.GTT.data)
			    NBIP.GTT            = jdata5.ppf.NBP2.GTT;
	    end
	    if ~isempty(jdata5.ppf.NBP2.GE.data) && ~isstr(jdata5.ppf.NBP2.GE.data)
			    NBIP.GE            = jdata5.ppf.NBP2.GE;
	    end
	    if ~isempty(jdata5.ppf.NBP2.GI.data) && ~isstr(jdata5.ppf.NBP2.GI.data)
			    NBIP.GI            = jdata5.ppf.NBP2.GI;
        end

        
        if isappdata(0,'TRANSP_USER') && isappdata(0,'TRANSP_RUN')
            revision = getappdata(0,'TRANSP_RUN');
            altuser = getappdata(0,'TRANSP_USER');
            rep             = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/curb',post.z0dinput.shot,altuser),[],revision);
            if ~isempty(rep)
                x_transp    = rep.object.dimensions.value.x1.value.data.value.data;  
                t_transp    = rep.object.dimensions.value.x0.value.data.value.data - 40;
                curb        = rep.object.data.value.data' / 1e4;
                curoh       = NaN * curb;
                curbs       = NaN * curb;
                rep         = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/pbe',post.z0dinput.shot,altuser),[],revision);
                pbe         = rep.object.data.value.data' / 1e6;
                rep         = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/pbi',post.z0dinput.shot,altuser),[],revision);
                pbi         = rep.object.data.value.data' / 1e6;
                rep         = sal_get(sprintf('data/pulse/%d/ppf/signal/%s/tqin',post.z0dinput.shot,altuser),[],revision);
                torque      = rep.object.data.value.data' / 1e6;
                torque_jxb  = NaN * torque;
                % missing data to compute  integrals
                jfast_transpt = NaN * ones(size(t_transp));
                jboot_transpt = NaN * ones(size(t_transp));
                johm_transpt = NaN * ones(size(t_transp));
            end
        else
            try
                % autres donneees transp
                try
                    % bypass error appening on some system where the first try return an error
                    try
                        mdsconnect('mdsplus.jet.efda.org');
                    catch
                        mdsconnect('mdsplus.jet.efda.org');
                    end
                    eval(sprintf('mdsopen(''trout_jet'',%s%s);',num2str(post.z0dinput.shot),num2str(id_transp)));
                catch
                    % bypass error appening on some system where the first try return an error
                    try
                        mdsconnect('mdsplus.jet.efda.org');
                    catch
                        mdsconnect('mdsplus.jet.efda.org');
                    end
                    eval(sprintf('mdsopen(''trout'',%s%s);',num2str(post.z0dinput.shot),num2str(id_transp)));
                end
                %% ---- Fast ion current ---- %%
                [curb,stat] = mdsvalue('_sig=\top.transp_out:curb');
                [x_transp,stat] = mdsvalue('dim_of(_sig,0)');
                [t_transp,stat] = mdsvalue('dim_of(_sig,1)');
                
                %% ---- Te and Ti profiles ---- %%
                [te_transp,stat] = mdsvalue('_sig=\top.transp_out:te');
                [ti_transp,stat] = mdsvalue('_sig=\top.transp_out:ti');
                
                %% ---- ne and ni profiles ---- %%
                [ne_transp,stat] = mdsvalue('_sig=\top.transp_out:ne');
                [ni_transp,stat] = mdsvalue('_sig=\top.transp_out:ni');
                
                %% ---- Toroidal flux and volume ---- %%
                [trflx,stat] = mdsvalue('_sig=\top.transp_out:trflx');
                [darea,stat] = mdsvalue('_sig=\top.transp_out:darea');
                
                %% ---- Major radius of magnetic axis ---- %%
                [raxis,stat] = mdsvalue('_sig=\top.transp_out:raxis');
                
                %% ---- Plasma current profile ---- %%
                [cur,stat] = mdsvalue('_sig=\top.transp_out:cur');
                
                %% ---- Ohmic current profile ---- %%
                [curoh,stat] = mdsvalue('_sig=\top.transp_out:curoh');
                
                %% ---- q-profile ---- %%
                [q_transp,stat] = mdsvalue('_sig=\top.transp_out:q');
                %  	idx001 = closest(x_transp,0.01);
                %  	idx005 = closest(x_transp,0.05);
                %  	idx095 = closest(x_transp,0.95);
                %  	for i=1:length(t_transp)
                %  	q001_transp(i) = q_transp(idx001,i);
                %  	q005_transp(i) = q_transp(idx005,i);
                %  	q095_transp(i) = q_transp(idx095,i);
                %  	end
                %  	qmin_transp = min(q_transp);
                %  	for i=1:length(t_transp)
                %  	kx=find(q_transp(1:end-5,i)==min(q_transp(1:end-5,i)));
                %  	if length(kx)>0
                %  	if length(kx) > 1
                %  		kx = kx(end);
                %  	end
                %  		xqmin_transp(i) = x_transp(kx);
                %  	else
                %  	xqmin_transp(i) = NaN;
                %  	end
                %  	end
                
                %% ---- Bootstrap current profile ---- %%
                [curbs,stat] = mdsvalue('_sig=\top.transp_out:curbs');
                
                %% ---- Zeff profile ---- %%
                [zeffi,stat] = mdsvalue('_sig=\top.transp_out:zeffi');
                
                %% ---- Power to electrons and ions (profile) ---- %%
                [pbi,stat] = mdsvalue('_sig=\top.transp_out:pbi');
                [pbe,stat] = mdsvalue('_sig=\top.transp_out:pbe');
                
                %% ---- Power to electrons and ions (vs time) ---- %%
                [bpti,stat] = mdsvalue('_sig=\top.transp_out:bpti');
                [bpte,stat] = mdsvalue('_sig=\top.transp_out:bpte');
                
                %% ---- Inductance Li vs time ---- %%
                [li_3,stat] = mdsvalue('_sig=\top.transp_out:li_3');
                
                %% ---- beta vs time ---- %%
                [betat,stat] = mdsvalue('_sig=\top.transp_out:betat');
                [li2pb,stat] = mdsvalue('_sig=\top.transp_out:li2pb');
                [pcurc,stat] = mdsvalue('_sig=\top.transp_out:pcurc');
                %  	a_transp  = interp1(data.gene.temps,data.equi.a(:,101),t_transp+40);
                %  	b0_transp = interp1(data.gene.temps,data.geo.b0,t_transp+40);
                %  	betan_transp = betat./(pcurc/1.e6).*a_transp.*b0_transp;
                
                %% ---- Plasma current vs time ---- %%
                [pcurc,stat] = mdsvalue('_sig=\top.transp_out:pcurc');
                
                %% ---- Fast ion, bootstrap and ohmic current versus time ---- %%
                jfast_transpt = NaN * ones(size(t_transp));
                jboot_transpt = NaN * ones(size(t_transp));
                johm_transpt = NaN * ones(size(t_transp));
                try
                    for i=1:length(t_transp)
                        jfast_transpt(i) = trapz(curb(:,i).*darea(:,i));
                        jboot_transpt(i) = trapz(curbs(:,i).*darea(:,i));
                        johm_transpt(i)  = trapz(curoh(:,i).*darea(:,i));
                    end
                end
                %% ---- Zeff vs time ---- %%
                [zeffc,stat] = mdsvalue('_sig=\top.transp_out:zeffc');
                
                %% ---- Vloop vs time ---- %%
                [vsurc,stat] = mdsvalue('_sig=\top.transp_out:vsurc');
                
                %% ---- Energy confinement time vs time ---- %%
                [taua1,stat] = mdsvalue('_sig=\top.transp_out:taua1');
                
                %% ---- torque ---- %%
                [torque,stat] = mdsvalue('_sig=\top.transp_out:tqin');
                [torque_jxb,stat] = mdsvalue('_sig=\top.transp_out:tqjxb');
                %
                %eval('mdsclose;');
            end
        end
    end
	    
    if ~isempty(TRA0.QBE.data) & ~isstr(TRA0.QBE.data)
	    dtransp  = interp1(post.profil0d.temps,post.profil0d.Raxe,TRA0.QBE.t);
	    atransp  = interp1(post.z0dinput.cons.temps,post.z0dinput.geo.a,TRA0.QBE.t);	
	    vt       = ones(size(atransp));
	    xx       = cat(2,-post.profil0d.xli(end:-1:2),post.profil0d.xli);
	    dd       = cat(2,dtransp(:,end:-1:2),dtransp);
	    ve       = ones(size(xx));
	    rr       = dd + (atransp * xx);
	    xtransp  = abs(tsplinet(rr,vt * xx,vt * TRA0.QBI.x(:)'));
    end
    h = findobj(0,'type','figure','tag','z0plotnbi');
    if isempty(h)
	  h=figure('tag','z0plotnbi');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    subplot(2,2,1)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,real(post.profil0d.pnbi - post.profil0d.pnbi_ion) + ...
              imag(post.profil0d.pnbi - post.profil0d.pnbi_ion),'color','r');
    leg = {'Metis'};
    try
        zplotprof(gca,NBIP.GE.t,NBIP.GE.x,NBIP.GE.data,'color','b');
        leg{end+1} ='Pencil';
    end
    try
	    zplotprof(gca,TRA0.QBE.t,xtransp,TRA0.QBE.data,'color','g','linestyle','none','marker','o');
	    leg{end+1} ='Transp';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e6 .* pbe','color','k');
	    leg{end+1} =sprintf('Transp (%s)',id_transp);
    end

    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.el,'color','c');
	    leg{end+1} ='Cronos';
    end
    xlabel('X_{Lao}')
    ylabel('P_{el} (Wm^{-3})');
    legend(leg)

    subplot(2,2,2)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,real(post.profil0d.pnbi_ion) + imag(post.profil0d.pnbi_ion),'color','r');
    leg = {'Metis'};
    try
        zplotprof(gca,NBIP.GI.t,NBIP.GI.x,NBIP.GI.data,'color','b');
        leg{end+1} ='Pencil';
    end
    try
	    zplotprof(gca,TRA0.QBI.t,xtransp,TRA0.QBI.data,'color','g','linestyle','none','marker','o');
	    leg{end+1} ='Transp';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e6 .* pbi','color','k');
	    leg{end+1} =sprintf('Transp (%s)',id_transp);
    end
    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.ion,'color','c');
	    leg{end+1} ='Cronos';
	    
	    % calcul de ecrit_nbi_slow pour cronos
	    fact = zeros(size(data.prof.ne));
	    for k=1:size(data.impur.impur,3)
		    fact           =  squeeze(data.impur.impur(:,:,k)) .* param.compo.z(k) .^ 2 ./ param.compo.a(k) + fact;
	    end
	    fact = fact ./ data.prof.ne;
	    ecrit_nbi_slow = 14.8 .* data.prof.te .* (2 .^ (3/2)  .* fact) .^ (2/3);
	    e0b  = mean(param.cons.idn.energie(param.cons.idn.energie > 1e3));
	    fr = zfract0(ecrit_nbi_slow,e0b);
	    zplotprof(gca,data.gene.temps,param.gene.x,(data.source.idn.ion + data.source.idn.el) .* fr,'color','g');
	    
    end
    xlabel('X_{Lao}')
    ylabel('P_{ion} (Wm^{-3})');
    legend(leg)

    subplot(2,2,3)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,real(post.profil0d.jnbicd) + imag(post.profil0d.jnbicd),'color','r');
    leg = {'Metis'};
    try
        zplotprof(gca,NBIP.JBDC.t,NBIP.JBDC.x,-NBIP.JBDC.data,'color','b');
        leg{end+1} ='Pencil';
    end
    try
	    zplotprof(gca,TRA0.CB.t,xtransp,TRA0.CB.data,'color','g','linestyle','none','marker','o');
	    leg{end+1} ='Transp';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e4 .* curb','color','k');
	    leg{end+1} =sprintf('Transp (%s)',id_transp);
    end
    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.j,'color','c');
	    leg{end+1} ='Cronos';
    end
    xlabel('X_{Lao}')
    ylabel('J (Am^{-2})');
    legend(leg)

    subplot(2,2,4)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,real(post.profil0d.pnbi) + imag(post.profil0d.pnbi),'color','r');
    leg = {'Metis'};
    try
        zplotprof(gca,NBIP.GE.t,NBIP.GE.x,NBIP.GE.data + NBIP.GI.data,'color','b');
        leg{end+1} ='Pencil';
    end
    try
	    zplotprof(gca,TRA0.QBE.t,xtransp,TRA0.QBE.data +TRA0.QBI.data ,'color','g','linestyle','none','marker','o');
	    leg{end+1} ='Transp';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e6 .* (pbe+ pbi)','color','k');
	    leg{end+1} =sprintf('Transp (%s)',id_transp);
    end
    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.el + data.source.idn.ion,'color','c');
	    leg{end+1} ='Cronos';
    end
    xlabel('X_{Lao}')
    ylabel('P (Wm^{-3})');
    legend(leg)



    h = findobj(0,'type','figure','tag','z0plotnbi_t');
    if isempty(h)
	  h=figure('tag','z0plotnbi_t');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    plot(post.zerod.temps,real(post.zerod.inbicd) ./1e6 + imag(post.zerod.inbicd) ./1e6,'r', ...
         NBIP.JBDT.t,-NBIP.JBDT.data ./1e6,'m', ...
	  post.zerod.temps,real(post.zerod.pnbi_th)./1e7 + imag(post.zerod.pnbi_th)./1e7,'b', ...
          NBIP.GTT.t,NBIP.GTT.data./1e7,'c')
    leg ={'I Metis', 'I Pencil','P_{th} Metis','P Pencil'};
    try
	hold on 
	plot(data.gene.temps,data.gene.iidn ./1e6,'g',data.gene.temps,data.gene.paddidn./1e7,'k');
	leg{end+1} ='I Cronos';
	leg{end+1} ='P Cronos';
    end
    xlabel('time (s)')
    ylabel('Ptot (10 MW) & I (MA)');
    legend(leg)


	
    h = findobj(0,'type','figure','tag','z0plotnbiq');
    if isempty(h)
	  h=figure('tag','z0plotnbiq');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.qjli,'color','r');
    leg = {'METIS'};
    zplotprof(gca,EFIT.Q.t,EFIT.Q.x,EFIT.Q.data,'color','b','marker','o','linestyle','none');
    leg{end+1} ='EFIT';
    zplotprof(gca,post.zerod.temps* ones(1,2),ones(size(post.zerod.temps,1),2),post.zerod.qeff* ones(1,2),'color','g','marker','*');
    leg{end+1} ='q_{eff}';
    try
	    zplotprof(gca,EQUI.Q.t,EQUI.Q.x,EQUI.Q.data,'color','c','marker','+','linestyle','none');
	    leg{end+1} ='EQUI';
    end 
    legend(leg);
    xlabel('x (normalized radius)');
    ylabel('safety factor');
     

    h = findobj(0,'type','figure','tag','z0plotnbiti');
    if isempty(h)
	  h=figure('tag','z0plotnbiti');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    rext  = post.profil0d.Raxe + interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'linear','extrap') * post.profil0d.xli;
    zplotprof(gca,post.profil0d.temps,rext,post.profil0d.tip ./ 1e3,'color','r');
    leg = {'METIS'};
    try
	    zplotprof(gca,CX.TI.t,CX.TI.x,CX.TI.data ./1e3,'color','b','marker','o','linestyle','none');
	    leg{end+1} ='CXFM';
    end 
    try
	    zplotprof(gca,CXS.TI.t,CXS.TI.x,CXS.TI.data ./1e3,'color','c','marker','+','linestyle','none');
	    leg{end+1} ='CXSM';
    end 
    try
	    zplotprof(gca,CXS.TICR.t,CXS.TICR.x,CXS.TICR.data ./1e3,'color','m','marker','+','linestyle','none');
	    leg{end+1} ='CXSM corrected';
    end 
     try
	    zplotprof(gca,CXG.TI.t,CXG.TI.x,CXG.TI.data ./1e3,'color','c','marker','x','linestyle','none');
	    leg{end+1} ='CXGM';
    end 
    try
	    zplotprof(gca,CXG.TICR.t,CXG.TICR.x,CXG.TICR.data ./1e3,'color','m','marker','x','linestyle','none');
	    leg{end+1} ='CXGM corrected';
    end 
   legend(leg);     
    xlabel('Rext (m)')
    ylabel('Ti (keV)');

    h = findobj(0,'type','figure','tag','z0plotnbiti2');
    if isempty(h)
	  h=figure('tag','z0plotnbiti2');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.tip ./ 1e3,'color','r');
    leg = {'METIS'};
    try
	    zplotprof(gca,CXS.TIRH.t,CXS.TIRH.x,CXS.TIRH.data ./1e3,'color','m','marker','+','linestyle','none');
	    leg{end+1} ='CXSM regularised';
    end 
    try
	    zplotprof(gca,CXG.TIRH.t,CXG.TIRH.x,CXG.TIRH.data ./1e3,'color','m','marker','x','linestyle','none');
	    leg{end+1} ='CXGM regularised';
    end 
    legend(leg);     
    xlabel('x (normalized radius)')
    ylabel('Ti (keV)');

    h = findobj(0,'type','figure','tag','z0plottran0ti');
    if isempty(h)
	  h=figure('tag','z0plottran0ti');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.tip ./ 1e3,'color','r');
    leg = {'METIS'};
    try
	    zplotprof(gca,TRA0.TI.t,TRA0.TI.x,TRA0.TI.data ./1e3,'color','b','marker','o','linestyle','none');
	    leg{end+1} ='TRANSP';
    end 
    legend(leg);     
    xlabel('x ')
    ylabel('Ti (keV)');

    h = findobj(0,'type','figure','tag','z0plotnbiomega');
    if isempty(h)
	  h=figure('tag','z0plotnbiomega');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    rext  = post.profil0d.Raxe + interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'linear','extrap') * post.profil0d.xli;
    zplotprof(gca,post.profil0d.temps,rext,post.profil0d.omega ./ 2 ./ pi ./ 1e3,'color','g');
    leg = {'METIS (bulk plasma)'};
    zplotprof(gca,post.profil0d.temps,rext,post.profil0d.vtor ./ 2 ./ pi ./ rext ./ 1e3,'color','r');
    leg{end+1} = 'METIS (main imurities)';
    try
	    zplotprof(gca,CX.ANGF.t,CX.ANGF.x,CX.ANGF.data ./ 2 ./ pi ./ 1e3,'color','b','marker','o','linestyle','none');
	    leg{end+1} ='CXFM';
    end 
    try
	    zplotprof(gca,CXS.ANGF.t,CXS.ANGF.x,CXS.ANGF.data ./ 2 ./ pi ./ 1e3,'color','c','marker','+','linestyle','none');
	    leg{end+1} ='CXSM';
    end 
    try
	    zplotprof(gca,CXG.ANGF.t,CXG.ANGF.x,CXG.ANGF.data ./ 2 ./ pi ./ 1e3,'color','c','marker','+','linestyle','none');
	    leg{end+1} ='CXGM';
    end 
    legend(leg);     
    xlabel('Rext (m)')
    ylabel('Omega (k turns / s)');
    %set(gca,'ylim',[min(post.profil0d.omega(:)),max(post.profil0d.omega(:))] ./ 2 ./ pi ./ 1e3);



    if ~isempty(BRAG.TE.data) && ~isempty(BRAG.TI.data)

        h = findobj(0,'type','figure','tag','z0plotnbi_brag');
        if isempty(h)
            h=figure('tag','z0plotnbi_brag');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',1,'color',[1 1 1])
        plot(post.profil0d.temps,post.profil0d.tep(:,1) ./ 1e3,'r',post.profil0d.temps,post.profil0d.tip(:,1) ./ 1e3,'b', ...
            BRAG.TE.t,BRAG.TE.data./1e3,'m',BRAG.TE.t,BRAG.TI.data./1e3,'c');
        leg ={'Te0 metis', 'Ti0 metis','Te_brag','Ti_brag'};
        xlabel('time (s)')
        ylabel('keV');
        try
            legend(leg)
        end
    end
    h = findobj(0,'type','figure','tag','z0plotnbi_torque');
    if isempty(h)
	  h=figure('tag','z0plotnbi_torque');
    else
	  figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',1,'color',[1 1 1])
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,real(post.profil0d.rot_nbi) + imag(post.profil0d.rot_nbi),'color','r');
    leg = {'Metis'};
    try
	    zplotprof(gca,NBIP.TORP.t,NBIP.TORP.x,-NBIP.TORP.data,'color','b');
	    leg{end+1} ='Pencil';
    end
    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.w,'color','c');
	    leg{end+1} ='Cronos';
    end
    try
	    zplotprof(gca,data.gene.temps,param.gene.x,data.source.idn.wb,'color','g');
	    leg{end+1} ='test Cronos';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e6 .* torque','color','k');
	    leg{end+1} ='Transp';
    end
    try
	    zplotprof(gca,t_transp+40,x_transp',1e6 .* torque_jxb','color','m');
	    leg{end+1} ='Transp (JxB)';
    end
    xlabel('x')
    ylabel('Torque (N m^-2)');
    legend(leg)


otherwise
    disp('These graphs can be plot only with JET data')
end


