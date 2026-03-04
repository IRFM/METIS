%script to converge operating point with TGLF prediction
%It implement  a convergence loop between METIS call in operation point 
%(computig plasma quantities at stationnary state for one time slice)
%and TGLF  profiles prediction:
%METIS compute sources, equilibrium and other quantities using TGLF prediction for Te, Ti and ne;
%and TGLF use updated METIS sources, equilibrium and other quantities to compute Te, Ti and ne.

% reset external data
%rm_cs4m;

% set parameters
waitfor( zuicreefunform('gui_tglf_op','tglfparams',1,0,'','visible'));
if getappdata(0,'ZGUI_CANCEL')
    return;
end

% first call
clear TGLF_data
klz_max = 31;
dwait_ = 1 ./ (2 + 2 .* klz_max);
fwait_ = 0;
if isappdata(0,'WORKING_POINT_TIME');
  rmappdata(0,'WORKING_POINT_TIME');
end
z0working_point
which z0working_point -all
hwaitbar = waitbar(fwait_,'Please wait...');
fwait_ = fwait_ + dwait_;
waitbar(fwait_,hwaitbar,'Please wait...');
setappdata(0,'WORKING_POINT_TIME',twp);
disp(twp)
% convergence loop with TGLF
outer_convergence_time  = 0
for klz_ = 1:klz_max
  t_start_outer = tic;
  % memorisation of last result
  jeux2.post = post;
  % TGLF prediction 
  [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0tglf(tglfparams,post.z0dinput.option,post.zerod,post.profil0d,2,klz_);
  whos chie_loop chii_loop chii_neo D_loop V_loop;
  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  drawnow
  % use external data for Te, Ti and n_e
  if tglfparams.calc_heat_transport == 1
      cs4m('Te','nointer');
      cs4m('Ti','nointer');
  end
  if tglfparams.calc_part_transport == 1
      cs4m('Ne','nointer');
  end

  % prediction for working point
  z0working_point
  which z0working_point -all
  fwait_ = fwait_ + dwait_;
  waitbar(fwait_,hwaitbar,'Please wait...');
  % test 
  somme = zcompstruct(jeux2.post.profil0d,post.profil0d,1e-2);
  if real(somme) < 3e-2
    break;
  end
outer_convergence_time = outer_convergence_time + toc(t_start_outer);
end
rmappdata(0,'WORKING_POINT_TIME');
delete(hwaitbar);
outer_it_count = klz_
fprintf('Outer convergence completed in %d iterations and %.2f seconds \n',outer_it_count,outer_convergence_time);
