%script to converge operating point with Qualikiz NN 10 D prediction
%It implement  a convergence loop between METIS call in operation point 
%(computig plasma quantities at stationnary state for on time slice)
%and QLK-NN 10D  profiles prediction:
%METIS compute sources, equilibrium and other quantities using QLK-NN 10D prediction for Te, Ti and ne;
%and QLK-NN 10D use updated METIS sources, equilibrium and other quantities to compute Te, Ti and ne.
%

% reset external data
%rm_cs4m;

% set parameters
waitfor( zuicreefunform('gui_qlkz_std','qlkzparams',1,0,'','visible'));
if getappdata(0,'ZGUI_CANCEL')
    return;
end

%qlkzparams.calc_part_transport = 1;
% fisrt call
clear QLK_nn_data
klz_max = qlkzparams.ext_nbloop;
dwait_ = 1 ./ (2 + 2 .* klz_max);
fwait_ = 0;
if isappdata(0,'WORKING_POINT_TIME');
  rmappdata(0,'WORKING_POINT_TIME');
end
z0working_point;
hwaitbar = waitbar(fwait_,'Please wait...');
fwait_ = fwait_ + dwait_;
waitbar(fwait_,hwaitbar,'Please wait...');
setappdata(0,'WORKING_POINT_TIME',twp);
% convergence loop with Qualikiz NN
for klz_ = 1:klz_max
  % memorisation of last result
  jeux2.post = post;
  list = findobj(0, 'type', 'figure');
  main = findobj(0, 'type', 'figure', 'tag', 'zeroda');
  
  for k = 1:length(list)
      % Vérifie si c'est ni la fenêtre principale, ni un waitbar
      is_waitbar = strcmp(get(list(k), 'Tag'), 'TMWWaitbar');
      if list(k) ~= main && ~is_waitbar
          close(list(k));
      end
  end

  % Vérification avant d'utiliser ho (évite erreur si supprimé)
  if exist('ho', 'var') && ishandle(ho)
      set(ho, 'value', 0);
  end
  % Qualikiz prediction 
  if qlkzparams.imp_trans == 4
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop,G_loop,D_loop_n,V_loop_n,G_loop_n] = z0qlkz_imp_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d,2,twp,post.z0dinput);

  elseif qlkzparams.imp_trans == 3
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop,G_loop,D_loop_n,V_loop_n,G_loop_n] = z0qlkz_imp_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d,2,twp,post.z0dinput);

  elseif qlkzparams.imp_trans == 2
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop,G_loop] = z0qlkz_imp_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d,2,twp,post.z0dinput);

  elseif qlkzparams.imp_trans == 1
    qlkzparams.nbloop_max=1;
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop,G_loop] = z0qlkz_imp_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d,2,twp,post.z0dinput);
    return
  elseif qlkzparams.imp_trans == 0
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0qlkz_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d,2);
  end

  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  drawnow
  % use external data for Te, Ti and n_e
  if qlkzparams.calc_heat_transport == 1
      cs4m('Te','nointer');
      cs4m('Ti','nointer');
  end
  if qlkzparams.calc_part_transport == 1
    cs4m('Ne','nointer');
    if qlkzparams.imp_trans ~= 0
      cs4m('N1','nointer');
      cs4m('NW','nointer');
    else
    end
  end

  % prediction for working point
  z0working_point;
  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  % test 
  somme = zcompstruct(jeux2.post.profil0d,post.profil0d,1e-2);
  if real(somme) < 5e-2
    break;
  end
end
rmappdata(0,'WORKING_POINT_TIME');
delete(hwaitbar);
