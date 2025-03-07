% test shot: 85425
%  liste          = {};
%  liste{end+1}   = 'ppf/@shot/KY63/ZEFF';    % Zeff Profile Li-beam
%  liste{end+1}   = 'ppf/@shot/KY63/TI';      % TI Li-beam
%  liste{end+1}   = 'ppf/@shot/KY63/TE';      %TE Li-beam
%  liste{end+1}   = 'ppf/@shot/KY63/RMID';    % R_mid_plane Li-beam
%  liste{end+1}   = 'ppf/@shot/KY63/NE';      % TI Li-beam
%  liste{end+1}   = 'ppf/@shot/KY63/NER';      %TE Li-beam
%  kY63         = cgcgetjet(post.z0dinput.shot,liste,'','');
liste          = {};
%  liste{end+1}   = 'ppf/@shot/KY6/ZEFF';    % Zeff Profile Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6/TI';      % TI Li-beam
liste{end+1}   = 'ppf/@shot/KY6/TE';      %TE Li-beam
liste{end+1}   = 'ppf/@shot/KY6/RMID';    % R_mid_plane Li-beam
liste{end+1}   = 'ppf/@shot/KY6/NE';      % TI Li-beam
liste{end+1}   = 'ppf/@shot/KY6/NER';      %TE Li-beam
KY6        = cgcgetjet(post.z0dinput.shot,liste,'','');
%  liste          = {};
%  liste{end+1}   = 'ppf/@shot/KY6S/ZEFF';    % Zeff Profile Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6S/TI';      % TI Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6S/TE';      %TE Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6S/RMID';    % R_mid_plane Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6S/NE';      % TI Li-beam
%  liste{end+1}   = 'ppf/@shot/KY6S/NER';      %TE Li-beam
%  kY6S        = cgcgetjet(post.z0dinput.shot,liste,'','');
KY6 = KY6.ppf.KY6;
h = findobj(0,'type','figure','tag','z0plot_jet_libeam');
if isempty(h)
      h=figure('tag','z0plot_jet_libeam');
else
      figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
subplot(2,2,1)
nelibeam =KY6.NE.data / 1e19;
nelibeam(nelibeam < eps) = NaN;
zplotprof(gca,KY6.RMID.t,KY6.RMID.data,nelibeam,'color','b');
amat = interp1(post.z0dinput.cons.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');
rli   = post.profil0d.Raxe + amat * post.profil0d.xli;
zplotprof(gca,post.profil0d.temps,rli,post.profil0d.nep./1e19,'color','r');
legend('Li beam','METIS');
xlabel('R (m)');
ylabel('n_e (10^{19} m^{-3})');
zplotprof(gca,KY6.RMID.t,KY6.RMID.data,nelibeam + KY6.NER.data / 1e19,'color','c');
zplotprof(gca,KY6.RMID.t,KY6.RMID.data,nelibeam - KY6.NER.data / 1e19,'color','c');
subplot(2,2,2)
zplotprof(gca,KY6.NE.t,KY6.NE.x,nelibeam,'color','b');
legend('Li beam');
xlabel('Z (m)');
ylabel('n_e (10^{19} m^{-3})');
zplotprof(gca,KY6.NE.t,KY6.NE.x,nelibeam + KY6.NER.data / 1e19,'color','c');
zplotprof(gca,KY6.NE.t,KY6.NE.x,nelibeam - KY6.NER.data / 1e19,'color','c');

subplot(2,2,3)
telibeam =KY6.TE.data / 1e3;
telibeam(telibeam < eps) = NaN;
zplotprof(gca,KY6.RMID.t,KY6.RMID.data,telibeam,'color','b');
amat = interp1(post.z0dinput.cons.temps,post.z0dinput.geo.a,post.profil0d.temps,'nearest');
rli   = post.profil0d.Raxe + amat * post.profil0d.xli;
zplotprof(gca,post.profil0d.temps,rli,post.profil0d.tep./1e3,'color','r');
legend('Li beam','METIS');
xlabel('R (m)');
ylabel('T_e (keV)');

subplot(2,2,4)
zplotprof(gca,KY6.TE.t,KY6.TE.x,telibeam,'color','b');
legend('Li beam');
xlabel('Z (m)');
ylabel('T_e (keV)');
