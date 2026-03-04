function nclass=read_neo_out_nclass(dirname)

try

file=load([dirname '/out.neo.theory_nclass']);

[aa bb] = size(file);

n_radial = aa;
nbspecies = (length(file(1,:))-2)./5;
for jr = 1:n_radial;

nclass.rova(jr) = file(jr,1);
nclass.BScur(jr) = file(jr,2); 
for js = 1:nbspecies;
nclass.pflux(jr,js) = file(jr,2+(js-1)*5+1);
nclass.qflux(jr,js) = file(jr,2+(js-1)*5+2);
nclass.upara(jr,js) = file(jr,2+(js-1)*5+3);
nclass.upol(jr,js)  = file(jr,2+(js-1)*5+4);
nclass.utor(jr,js)  = file(jr,2+(js-1)*5+5);
end
end;

catch

n_radial = 1;
nbspecies = 6;
for jr = 1:n_radial;

nclass.rova(jr) = NaN;
nclass.BScur(jr) = NaN;
for js = 1:nbspecies;
nclass.pflux(jr,js) = NaN;
nclass.qflux(jr,js) = NaN;
nclass.upara(jr,js) = NaN;
nclass.upol(jr,js)  = NaN;
nclass.utor(jr,js)  = NaN;
end;
end;

end;
