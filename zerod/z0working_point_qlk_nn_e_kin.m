% script to converge operating point with Qualikiz NN prediction
% question on limit in which QLK is used 
warndlg('This tool use a prototype version of neural network interpolation of QUALIKIZ; please consider results with waryness','Disclamer');
disp('=============================================================')
fprintf('This tool use a prototype version of neural network interpolation of QUALIKIZ;\nplease consider results with waryness\n');
disp('=============================================================')
prompt={'shift_chi_e:','shift_chi_i:','shift_DV:','threshold:','neo_ion_in_D:','neo_ion_in_Chie:','factor_neo_chii:'};
name='QLKNN parameters';
numlines=1;
defaultanswer={'0','0','0','0','0','0','1'};

answer=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(answer)
    return
end
option = post.z0dinput.option;
option.shift_chi_e  = str2num(answer{1});
option.shift_chi_i  = str2num(answer{2});
option.shift_DV     = str2num(answer{3});
option.seuil        = str2num(answer{4});
option.neo_ion_in_D = str2num(answer{5});
option.neo_ion_in_Chie = str2num(answer{6});
option.factor_neo_chii = str2num(answer{7});

                       

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
  z0qlkANN_kin_e_2018(option,post.zerod,post.profil0d); 
  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  drawnow
  % use external data for Te, Ti and n_e
  cs4m('Te','nointer');
  cs4m('Ti','nointer');
  cs4m('Ne','nointer');
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
