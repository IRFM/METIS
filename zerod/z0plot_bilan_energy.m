noms = fieldnames(getappdata(0));
extonoff = 0;
for k = 1:length(noms)
    if findstr(noms{k},'_EXP')
        extonoff =1;
    end
end
if extonoff == 1
    
    ee     =  0.1602176462e-18;
    wel    = (3/2) .* ee .* trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.nep .* post.profil0d.tep,2)/1e6;
    wion   = (3/2) .* ee .* trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.nip .* post.profil0d.tip,2)/1e6;
    
    h = findobj(0,'type','figure','tag','z0plot_ce');
    if isempty(h)
        h=figure('tag','z0plot_ce');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',1,'color',[1 1 1])
    
    plot(post.zerod.temps,post.zerod.wth/1e6,'k',post.profil0d.temps,wel,'r',post.profil0d.temps,wion,'b',post.profil0d.temps,wel + wion,'g');
    title(sprintf('METIS : %s@%d/bilan energy', ...
        post.z0dinput.machine,post.z0dinput.shot));
    xlabel('time (s)');
    ylabel('MJ')
    legend('W_{th,scaling}','W_{el}','W_{ion}','W_{el} + W_{ion}');
    edition2
    
end