function ntheory=read_neo_out_theory(dirname)

try

file=load([dirname '/out.neo.theory']);

[aa bb] = size(file);

n_radial = aa;
nbspecies = (length(file(1,:))-16)./2;
for jr = 1:n_radial;
ntheory.rova(jr,1) = file(jr,1);
ntheory.HHpflux(jr,1) = file(jr,2); 
ntheory.HHqfluxi(jr,1) = file(jr,3); 
ntheory.HHqfluxe(jr,1) = file(jr,4); 
ntheory.HHBScur(jr,1) = file(jr,5); 
ntheory.HHKcoef(jr,1) = file(jr,6); 
ntheory.HHupara(jr,1) = file(jr,7); 
ntheory.HHupol0(jr,1) = file(jr,8); 
ntheory.CHqfluxi(jr,1) = file(jr,9); 
ntheory.TGqfluxi(jr,1) = file(jr,10); 
ntheory.SBScur(jr,1) = file(jr,11); 
ntheory.SKcoef(jr,1) = file(jr,12); 
ntheory.Supara(jr,1) = file(jr,13); 
ntheory.Supol0(jr,1) = file(jr,14); 
ntheory.HRphi1(jr,1) = file(jr,15); 
for js = 1:nbspecies;
ntheory.HSpflux(jr,js) = file(jr,15+(js-1)*2+1); 
ntheory.HSqflux(jr,js) = file(jr,15+(js-1)*2+2); 
end
ntheory.KBScur(jr,1) = file(jr,bb); 
end;

catch

n_radial = 1; nbspecies = 3;
for jr = 1:n_radial;
ntheory.rova(jr,1) = NaN;
ntheory.HHpflux(jr,1) = NaN;
ntheory.HHqfluxi(jr,1) = NaN;
ntheory.HHqfluxe(jr,1) = NaN;
ntheory.HHBScur(jr,1) = NaN;
ntheory.HHKcoef(jr,1) = NaN;
ntheory.HHupara(jr,1) = NaN;
ntheory.HHupol0(jr,1) = NaN;
ntheory.CHqfluxi(jr,1) = NaN;
ntheory.TGqfluxi(jr,1) = NaN;
ntheory.SBScur(jr,1) = NaN;
ntheory.SKcoef(jr,1) = NaN;
ntheory.Supara(jr,1) = NaN;
ntheory.Supol0(jr,1) = NaN;
ntheory.HRphi1(jr,1) = NaN;
for js = 1:nbspecies;
ntheory.HSpflux(jr,js) = NaN;
ntheory.HSqflux(jr,js) = NaN;
end
ntheory.KBScur(jr,1) = NaN;
end;

end;
