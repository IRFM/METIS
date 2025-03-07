% test de zpelletcoor
tsscenario
disp('select a time');
[t,void] = ginput(1);
ind   = find(data.gene.temps >= t,1);
datak= zget1t(data,ind);
origine.r0 = double(max(datak.equi.R(:)));
origine.r1 = double(min(datak.equi.R(:)));
origine.z0 = double(max(datak.equi.Z(:)));
origine.z1 = double(min(datak.equi.Z(:)));
pos        = [0.1,0.5,1,5] .* datak.equi.a(end);
x          = param.gene.x;
[inout,rho,R,Z,volume,ne,te,nions,tions,gradbob,indice]=zpelletcoor(origine,pos,x,datak.equi,datak.prof,datak.source,datak.impur,datak.neo,1,'cubic')