% test de zsauter0d
rmx = data.equi.rhomax * param.gene.x;
epsi = rmx ./ data.equi.raxe; 

jboot = zsauter0d(param.gene.x,data.prof.te,data.prof.ti,data.prof.ne, ....
                  data.prof.ni,data.prof.q,data.gene.zeffm,rmx,data.equi.raxe,data.equi.ftrap,epsi,data.geo.b0,data.equi.e)
