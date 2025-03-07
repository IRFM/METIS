% benchmark avec cronos
% load CRONOS data (zuiload)
% create a METIS data set from CRONOS data
% use external data for CRONOS for ne,Te and Ti.
% run METIS
% then call this script
rm    = data.equi.rhomax * ones(size(param.gene.x));
qe_c  = rm .* cumtrapz(param.gene.x,data.source.totale.el .* data.equi.vpr,2);
qi_c  = rm .* cumtrapz(param.gene.x,data.source.totale.ion .* data.equi.vpr,2);
qei_c = rm .* cumtrapz(param.gene.x,data.source.qei .* data.equi.vpr,2);
ge    = rm .* cumtrapz(param.gene.x,data.source.totale.ne .* data.equi.vpr,2) ./ ...
        data.equi.grho2 ./ data.equi.vpr;
ge(:,1) = 0;
tt    = data.gene.temps * ones(size(param.gene.x));
xx    = ones(size(data.equi.rhomax)) * param.gene.x;
tm    = post.profil0d.temps * ones(size(post.profil0d.xli));
rmx   = post.profil0d.rmx ./ (max(post.profil0d.rmx,[],2) * ones(size(post.profil0d.xli))); 
post.profil0d.qe = griddata(tt,xx,qe_c,tm,rmx,'cubic');
post.profil0d.qi = griddata(tt,xx,qi_c,tm,rmx,'cubic');
post.profil0d.qei = griddata(tt,xx,qei_c,tm,rmx,'cubic');
post.profil0d.ge = griddata(tt,xx,ge,tm,rmx,'cubic');
post.profil0d.tep = griddata(tt,xx,data.prof.te,tm,rmx,'cubic');
post.profil0d.tip = griddata(tt,xx,data.prof.ti,tm,rmx,'cubic');
post.profil0d.nep = griddata(tt,xx,data.prof.ne,tm,rmx,'cubic');
post.profil0d.nip = griddata(tt,xx,data.prof.ni,tm,rmx,'cubic');
post.profil0d.ware = griddata(tt,xx,data.neo.coef.vn,tm,rmx,'cubic');
post.profil0d.chii_neo = griddata(tt,xx,data.neo.coef.ii ./ data.prof.ni,tm,rmx,'cubic');
post.profil0d.xie = griddata(tt,xx,data.coef.ee ./ data.prof.ne + data.neo.coef.ee ./ data.prof.ne,tm,rmx,'cubic');
post.profil0d.xii = griddata(tt,xx,data.coef.ii ./ data.prof.ni + data.neo.coef.ii ./ data.prof.ni,tm,rmx,'cubic');
post.profil0d.dn = griddata(tt,xx,data.coef.nn + data.neo.coef.nn ,tm,rmx,'cubic');
post.profil0d.vn = griddata(tt,xx,data.coef.vn + data.neo.coef.vn,tm,rmx,'cubic');
post.profil0d.qjli = griddata(tt,xx,data.equi.q,tm,rmx,'cubic');
%post.profil0d.grho = post.profil0d.grho2;
post.z0dinput.geo.r0 = interp1(data.gene.temps,data.equi.raxe(:,end),'pchip');
post.z0dinput.geo.a  = interp1(data.gene.temps,data.equi.a(:,end),'pchip');
post.z0dinput.geo.b0 = interp1(data.gene.temps,data.equi.F(:,end) ./data.equi.raxe(:,end),'pchip');


z0qlkANN_kin_e_2018(post.z0dinput.option,post.zerod,post.profil0d); 



