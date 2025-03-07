% test d'un autre calcul dejmoy
dpsidx  = pdederive(param.gene.x,data.prof.psi,0,2,2,1);
dptotdx  = pdederive(param.gene.x,data.equi.ptot,0,2,2,1);
df2dx  = pdederive(param.gene.x,data.equi.F,0,2,2,1);
jmoy    = (dptotdx ./ dpsidx + 1./2./param.phys.mu0.* df2dx .* data.equi.r2i) ./ data.equi.ri;


