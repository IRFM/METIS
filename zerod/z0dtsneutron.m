recompute_ndd;
other_graph = false;

switch post.z0dinput.machine
    case 'TS'
        % script de comparaison dt taux de neutron pour TS
        % lecture flux neutron
        [fn,tfn,rfn,cfn]= tsbase(fix(post.z0dinput.shot),'gfluntn');
        zs      = post.zerod;
        profli  = post.profil0d;
        t       = zs.temps;
        xfn     = linspace(0,1,21);
        vt      = ones(size(t));
        ve      = ones(size(xfn));
        xxfn    = vt * xfn;
        nep = pchip(post.profil0d.temps',post.profil0d.nep',zs.temps')';
        tep = pchip(post.profil0d.temps',post.profil0d.tep',zs.temps')';
        if post.z0dinput.option.gaz  == 4
            nd  = ((zs.nim ./ zs.nem ./ 2) * ve) .* nep;
        else
            nd  = ((zs.nDm ./ zs.nem) * ve) .* nep;
        end
        tip  = pchip(post.profil0d.temps',post.profil0d.tip',zs.temps')';
        rm  =  pchip(post.profil0d.temps',post.profil0d.Raxe',zs.temps')';
        %fndd = fluxneutron(tip ./ 1e3,nd,xfn,rm,post.z0dinput.geo.a);
        fndd   = zs.ndd;
        
        % total shot neutron
        fprintf('Total neutron for Tore supra shot %d = %g\n',post.z0dinput.shot,trapz(zs.temps,zs.ndd));
        
        % donnees bragg
        [tib,tti] =tsbase(fix(post.z0dinput.shot),'stibrag');
        [etib,tti] =tsbase(fix(post.z0dinput.shot),'setibrag');
        nom=tsbase(fix(post.z0dinput.shot),'SCIONBRAG');
        [deplw,tw]=tsbase(fix(post.z0dinput.shot),'SDEPLW');
        [edeplw,tw]=tsbase(fix(post.z0dinput.shot),'SEDEPLW');
        if strcmp(lower(deblank(nom)),'fer')
            fact =3;
        elseif  strcmp(lower(deblank(nom)),'chrome')
            fact = 2.5;
        else
            fact =NaN;
        end
        % Ip sens inverse trigo & deplw >0 sens trigo , choc >28424
        if fix(post.z0dinput.shot)>28424
            depl0dm =  -zs.wrad .* post.z0dinput.geo.R  ./ 1e8 ./ fact;
            depl0d  =  -interp1(profli.temps,profli.vtor(:,1),t,'nearest')  ./ 1e8 ./ fact;
        else
            depl0dm = zs.wrad .* post.z0dinput.geo.R  ./ 1e8 ./ fact;
            depl0d  = interp1(profli.temps,profli.vtor(:,1),t,'nearest')  ./ 1e8 ./ fact;
        end
        
        hz =findobj(0,'type','figure','tag','0dneutron');
        if isempty(hz)
            hz=figure('tag','0dneutron','name','Neutron TS');
        else
            figure(hz);
        end
        clf
        set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',1,'color',[1 1 1])
        
        if isempty(tw)
            k = 2;
        else
            k = 3;
        end
        
        subplot(k,1,1);
        semilogy(t,fndd,'r',tfn,10 .^ fn,'c',tfn,10 .^ sgolayfilt(fn,1,23),'k');
        legend('0D','measurement 1','measurement 2','filtered measurement')
        %xlabel('time (s)');
        ylabel('Flux de neutron (s^-1)')
        axis([min(t),max(t),1e8,inf]);
        title(sprintf('Zerod : %s@%d/neutron DD ( r -> 0D,b -> mesure)', ...
            post.z0dinput.machine,post.z0dinput.shot));
        z0loglin(gca);
        drawnow
        xl= get(gca,'xlim');
        
        subplot(k,1,2);
        plot(t,tip(:,1)/1e3,'r');
        hold on
        plot(t,tip(:,5)/1e3,':r');
        errorbar(tti,tib,etib,'b');
        hold off
        %xlabel('time (s)');
        ylabel('Ti (keV)')
        set(gca,'xlim',xl);
        drawnow
        yl =get(gca,'ylim');
        set(gca,'ylim',[0 min(6,yl(2))])
        
        if ~isempty(tw)
            subplot(k,1,3);
            
            indin = find((t >min(tw)) & (t <max(tw)));
            delta =  mean(deplw) - mean(depl0d(indin));
            deltam =  mean(deplw) - mean(depl0dm(indin));
            
            plot(t,depl0d+delta,'r');
            hold on
            plot(t,depl0dm+deltam,'g');
            errorbar(tw,deplw,edeplw,'b');
            hold off
            xlabel('time (s)');
            ylabel('deplacement (A)')
            set(gca,'xlim',xl);
            set(gca,'ylim',[min(deplw)-max(edeplw) max(deplw)+max(edeplw)])
        end
        drawnow
        
    case 'JET'
        
        zs      = post.zerod;
        profli  = post.profil0d;
        t       = zs.temps;
        % total shot neutron
        fprintf('Total DD neutron for JET shot %d = %g\n',post.z0dinput.shot,trapz(zs.temps,zs.ndd));
        
        % caclul complet du beam beam
        %  	fprintf('BB pp:');
        %  	% reechantillone sur les profils
        %  	vpr             = interp1(profli.temps,profli.vpr,t,'nearest');
        %  	nep             = interp1(profli.temps,profli.nep,t,'nearest');
        %  	tep             = interp1(profli.temps,profli.tep,t,'nearest');
        %  	n1p             = interp1(profli.temps,profli.n1p,t,'nearest');
        %  	tip             = interp1(profli.temps,profli.tip,t,'nearest');
        %  	zeffp            = interp1(profli.temps,profli.zeff,t,'nearest');
        %  	fpnbi           = interp1(profli.temps,profli.pnbi,t,'nearest');
        %  	if isfield(zs,'einj_nbi_icrh')
        %  		zs.ndd_nbi_nbi_pp  = z0beambeam(zs.temps,profli.xli,vpr,nep,tep,n1p,tip,zeffp,fpnbi, ....
        %                         	     		zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem, ...
        %  			     		zs.einj_nbi_icrh,zs.einj_nbi_icrh,zs.pnbi_icrh + zs.pnbi_th,...
        %  			     		zs.pnbi_icrh + zs.pnbi_th, zs.mu0_nbi,1);
        %  	else
        %  		zs.ndd_nbi_nbi_pp  = z0beambeam(zs.temps,profli.xli,vpr,nep,tep,n1p,tip,zeffp,fpnbi, ....
        %                         	     		zs.nDm,zs.nTm,zs.n1m,zs.nhem,zs.nimpm,zs.nem, ...
        %  			     		post.z0dinput.option.einj,post.z0dinput.option.einj,max(post.z0dinput.cons.pnbi),...
        %  			     		max(post.z0dinput.cons.pnbi), zs.mu0_nbi,1);
        %  	end
        
        h = findobj(0,'type','figure','tag','neutron');
        if isempty(h)
            h=figure('tag','neutron');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',1,'color',[1 1 1])
        
        % 	zs.ndd_total = zs.ndd_th + zs.ndd_nbi_th  + zs.ndd_nbi_nbi_pp;
        
        plot(zs.temps,zs.ndd,'k', ...
            zs.temps,zs.ndd_th,'g', ...
            zs.temps,zs.ndd_nbi_th,'b',zs.temps,zs.ndd_nbi_nbi,'c', ...
            post.z0dinput.exp0d.temps,post.z0dinput.exp0d.ndd,'r',...
            post.z0dinput.exp0d.temps,post.z0dinput.exp0d.ndd_th,'k:',...
            post.z0dinput.exp0d.temps,post.z0dinput.exp0d.ndd_nbi_th,'m');
        %, ...
        %zs.temps,zs.ndd_nbi_nbi_pp,'m');%, ...
        %zs.temps,zs.ndd,'k:' );
        
        xlabel('time (s)');
        ylabel('neutron/s');
        title(sprintf('METIS : %s@%d/neutron DD ', ...
            post.z0dinput.machine,post.z0dinput.shot));
        legend('Total','Thermal','Posifif : beam-plasma','Posifif : beam-beam','Experimental total', ...
            'Pencil thermal','Pencil Beam-Plasma');
        z0loglin(gca);
        
        
        % add graph for other neutrons if tritium in the discharge
        if post.z0dinput.option.gaz == 3
            other_graph = true;
            
            h = findobj(0,'type','figure','tag','neutron2');
            if isempty(h)
                h=figure('tag','neutron2');
            else
                figure(h);
            end
            clf
            set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
                'defaultlinelinewidth',1,'color',[1 1 1])
            
            
        end
    case 'WEST'
        h = findobj(0,'type','figure','tag','neutron');
        if isempty(h)
            h=figure('tag','neutron');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',1,'color',[1 1 1])
        
        zs   = post.zerod;
        % total shot neutron
        fprintf('Total neutron for WEST shot %d = %g\n',post.z0dinput.shot,trapz(zs.temps,zs.ndd));
        plot(zs.temps,zs.ndd,'b',zs.temps,zs.ndd_th,'k:');
        z0loglin(gca);
        xlabel('time (s)');
        ylabel('neutron/s');
        title(sprintf('METIS : %s@%d/neutron DD ', ...
            post.z0dinput.machine,post.z0dinput.shot));
        if ~all(post.z0dinput.exp0d.ndd == 0)
            hold on
            plot(zs.temps,post.z0dinput.exp0d.ndd,'r');
            legend('Total','Thermal','Experimental total');
        else
            legend('Total','Thermal');
            
        end
    otherwise
        h = findobj(0,'type','figure','tag','neutron');
        if isempty(h)
            h=figure('tag','neutron');
        else
            figure(h);
        end
        clf
        set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
            'defaultlinelinewidth',1,'color',[1 1 1])
        
        if post.z0dinput.option.gaz ~= 3
            zs   = post.zerod;
            plot(zs.temps,zs.ndd,'k', ...
                zs.temps,zs.ndd_th,'g', ...
                zs.temps,zs.ndd_nbi_th,'b',zs.temps,zs.ndd_nbi_nbi,'c');
            z0loglin(gca);
            xlabel('time (s)');
            ylabel('neutron/s');
            title(sprintf('METIS : %s@%d/neutron DD ', ...
                post.z0dinput.machine,post.z0dinput.shot));
            if ~all(post.z0dinput.exp0d.ndd == 0)
                hold on
                plot(zs.temps,post.z0dinput.exp0d.ndd,'r');
                legend('Total','Thermal','Beam-Plasma','Beam-Beam','Experimental total');
            else
                legend('Total','Thermal','Beam-Plasma','Beam-Beam');
            end
        else
            other_graph =true
        end
        
end

if other_graph
    
    if length(post.profil0d.temps) == length(post.zerod.temps)
        zs   = post.zerod;
        geo  = post.z0dinput.geo;
        cons = post.z0dinput.cons;
    else
        noms = fieldnames(post.zerod);
        temps = post.zerod.temps;
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.zerod.(nomc);
            if length(var) == length(temps)
                zs.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                zs.(nomc) = var;
            end
        end
        zs.temps = post.profil0d.temps;
        noms = fieldnames(post.z0dinput.cons);
        temps = post.z0dinput.cons.temps;
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.z0dinput.cons.(nomc);
            if length(var) == length(temps)
                cons.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                cons.(nomc) = var;
            end
        end
        cons.temps = post.profil0d.temps;
        noms = fieldnames(post.z0dinput.geo);
        for k=1:length(noms)
            nomc = noms{k};
            var  = post.z0dinput.geo.(nomc);
            if length(var) == length(temps)
                geo.(nomc) = interp1(temps,var,post.profil0d.temps,'nearest','extrap');
            else
                geo.(nomc) = var;
            end
        end
    end
    [~,salpha,~,~,~,~,~,~,~,~,~,~,~,splustd,splusdt,splusff,splusicrh] = ...
        zfus0tae(zs.nDm,zs.nTm,zs.tem,zs.nem,zs.zeff,zs.tite,geo.R,geo.a,geo.K,geo.b0,zs.ane,zs.ate,zs.vp,zs.sp, ...
        zs.pnbi_th,zs.taus_nbi,zs.ecrit_nbi,post.z0dinput.option.einj,post.z0dinput.option.einj2 ,cons.ftnbi, ...
        zs.pion_icrh,zs.taus_icrh,zs.ecrit_nbi,zs.einj_icrh,post.profil0d.temps,cons.pnbi,zs.d0,zs.qa,zs.qmin, ...
        zs.te0,zs.nebord,zs.tebord,zs.pped,zs.nim,zs.wth, ...
        post.z0dinput.option.tae,post.z0dinput.option.nb_nbi,post.z0dinput.option.fspot,post.z0dinput.option.e_shielding,post.profil0d, ...
        post.z0dinput.option.fpolarized,post.z0dinput.option.forced_H_NBI);
    
    [neutron_total_tt,neutron_th_tt,neutron_nbi_th_tt,neutron_nbi_nbi_tt,pttfus,proton_tt,picrh_nbi_tt,einj_tt] = ...
        z0neutron_tt(post.z0dinput.option,cons,zs,post.profil0d);
    
    
    k = 3;
    subplot(k,1,1)
    plot(zs.temps,zs.ndd,'k', ...
        zs.temps,zs.ndd_th,'r', ...
        zs.temps,zs.ndd_nbi_th,'b',zs.temps,zs.ndd_nbi_nbi,'c');
    z0loglin(gca);
    xlabel('time (s)');
    ylabel('neutron/s (DD)');
    title(sprintf('METIS : %s@%d/neutron', ...
        post.z0dinput.machine,post.z0dinput.shot));
    legend('Total','Thermal','Beam-Plasma','Beam-Beam');
    
    % total shot neutron
    fprintf('Total DD neutron for %s shot %d = %g\n',post.z0dinput.machine,post.z0dinput.shot,trapz(zs.temps,zs.ndd));
    
    subplot(k,1,2)
    plot(zs.temps,salpha .* (zs.pfus + zs.pfus_loss) ./ zs.pfus,'k', ...
        zs.temps,salpha .* (zs.pfus + zs.pfus_loss) ./ zs.pfus - (splustd + splusdt + splusff + splusicrh),'r', ...
        zs.temps,splustd + splusdt,'b',zs.temps,splusff,'c',zs.temps,splusicrh,'m');
    z0loglin(gca);
    xlabel('time (s)');
    ylabel('neutron/s (DT)');
    legend('Total','Thermal','Beam-Plasma','Beam-Beam','ICRH induced');
    
    % total shot neutron
    fprintf('Total DT neutron for %s shot %d = %g\n',post.z0dinput.machine,post.z0dinput.shot,trapz(zs.temps,salpha));
    
    subplot(k,1,3)
    plot(zs.temps,neutron_total_tt,'k', ...
        zs.temps,neutron_th_tt,'r', ...
        zs.temps,neutron_nbi_th_tt,'b',zs.temps,neutron_nbi_nbi_tt,'c');
    z0loglin(gca);
    xlabel('time (s)');
    ylabel('neutron/s (TT)');
    legend('Total','Thermal','Beam-Plasma','Beam-Beam');
    
    % total shot neutron
    fprintf('Total TT neutron for %s shot %d = %g\n',post.z0dinput.machine,post.z0dinput.shot,trapz(zs.temps,neutron_total_tt));
    
    joint_axes(h,k);
    
    switch post.z0dinput.machine
        case 'JET'
            h = findobj(0,'type','figure','tag','neutron3');
            if isempty(h)
                h=figure('tag','neutron3');
            else
                figure(h);
            end
            clf
            set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
                'defaultlinelinewidth',1,'color',[1 1 1])
            
            % 	zs.ndd_total = zs.ndd_th + zs.ndd_nbi_th  + zs.ndd_nbi_nbi_pp;
            
            plot(zs.temps,zs.ndd + neutron_total_tt + salpha,'k', ...
                 post.z0dinput.exp0d.temps,post.z0dinput.exp0d.ndd,'r');
            
            xlabel('time (s)');
            ylabel('neutron/s');
            title(sprintf('METIS : %s@%d/Neutrons total ', ...
                post.z0dinput.machine,post.z0dinput.shot));
            legend('Total','Experimental total');
            z0loglin(gca);
            
            
            % try to compare with transp
            try
                jet_neutron_transp(post);
            catch
                % TRANSP ?
            end
            
    end

    
end

