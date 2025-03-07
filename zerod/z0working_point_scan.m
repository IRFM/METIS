% script to scan a parameter with working point mode
% 1/ load some reference simulation
% metis_load ...
% time on interrest (s)
twp = 800;   


% fisrt call
if isappdata(0,'WORKING_POINT_TIME');
  rmappdata(0,'WORKING_POINT_TIME');
end
setappdata(0,'WORKING_POINT_TIME',twp);  
% memorisation of last result
jeux2.post = post;
nbval = 3; %30
xiioxie_vect = logspace(-1,1,nbval);
rep = {};
tep = NaN * ones(nbval,length(post.profil0d.xli));
tip = NaN * ones(nbval,length(post.profil0d.xli));
pfus = NaN * ones(nbval,1);
% convergence loop with Qualikiz NN
for kzl_ = 1:nbval
  % scan parameter
  z0dinput = jeux2.post.z0dinput;
  z0dinput.option.xiioxie = xiioxie_vect(kzl_);
  % prediction for working point
  z0working_point;
  rep1.zerod = zerod_get1t(post.zerod,3);
  rep1.profil0d = zerod_get1t(post.profil0d,3);
  tep(kzl_,:) = rep1.profil0d.tep;
  tip(kzl_,:) = rep1.profil0d.tip;
  pfus(kzl_) = rep1.zerod.pfus;
  rep{kzl_} = rep1;
  
end
rmappdata(0,'WORKING_POINT_TIME');


figure;
plot(post.profil0d.xli,tep,'r',post.profil0d.xli,tip,'b')
pfus
