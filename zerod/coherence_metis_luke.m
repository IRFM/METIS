ind = post.lukeinmetis.computation_time_slice
time = post.profil0d.temps(post.lukeinmetis.computation_time_slice);
ind0d = find(post.zerod.temps >= time,1);		     
disp(' ')
disp('ieccd')
trapz(post.profil0d.xli,post.profil0d.jeccd(ind,:).* post.profil0d.spr(ind,:))
post.zerod.ieccd(ind0d)
trapz(post.profil0d.xli,(post.lukeinmetis.j_TOT - post.lukeinmetis.j_OHM).* post.profil0d.spr(ind,:))
post.lukeinmetis.I_TOT - post.lukeinmetis.I_OHM
disp(' ')
disp('Pohm')
trapz(post.profil0d.xli,post.profil0d.pohm(ind,:).* post.profil0d.vpr(ind,:))
trapz(post.profil0d.xli,post.profil0d.ej(ind,:).* post.profil0d.vpr(ind,:))
post.zerod.pohm(ind0d)
trapz(post.profil0d.xli,post.lukeinmetis.p_OHM.* post.profil0d.vpr(ind,:))
post.lukeinmetis.P_OHM

% test volume 
disp(' ')
disp('volume')
trapz(post.profil0d.xli,post.profil0d.vpr(ind,:))
post.zerod.vp(ind0d)
[vps,sps,sexts,peris,geo] = zgeo0(post.z0dinput.geo);
vps(ind0d)
% test surface 
disp(' ')
disp('surface')
trapz(post.profil0d.xli,post.profil0d.spr(ind,:))
post.zerod.sp(ind0d)
[vps,sps,sexts,peris,geo] = zgeo0(post.z0dinput.geo);
sps(ind0d)
