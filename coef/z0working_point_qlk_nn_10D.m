%script to converge operating point with Qualikiz NN 10 D prediction
%It implement  a convergence loop between METIS call in operation point 
%(computig plasma quantities at stationnary state for on time slice)
%and QLK-NN 10D  profiles prediction:
%METIS compute sources, equilibrium and other quantities using QLK-NN 10D prediction for Te, Ti and ne;
%and QLK-NN 10D use updated METIS sources, equilibrium and other quantities to compute Te, Ti and ne.
%

% not useful as Te Ti and n_e will be replace after the first QLK-NN
% computation. Allows to keep external source in the loop.
%% reset external data
%rm_cs4m;

% set parameters
waitfor( zuicreefunform('gui_qlkANN10D','qlkparams',1,0,'','visible'));
if getappdata(0,'ZGUI_CANCEL')
    return;
end

% fisrt call
clear QLK_nn_data
klz_max = 31;
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
  % Qualikiz prediction 
  [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0qlkANN_10D(qlkparams,post.z0dinput.option,post.zerod,post.profil0d,2);
  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  drawnow
  % use external data for Te, Ti and n_e
  if qlkparams.calc_heat_transport == 1
      cs4m('Te','nointer');
      cs4m('Ti','nointer');
  end
  if qlkparams.calc_part_transport == 1
      cs4m('Ne','nointer');
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
