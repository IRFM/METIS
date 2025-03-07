% version matricielle du calcul des section efficaces
function [si1s,si2s,scx,srec,sii,sss,Ass] = z0sectionh(te,ti)

% securite 
sl = size(te);
te = max(0.1,real(te(:)));
ti = max(0.1,real(ti(:)));

% ionisation de H (1s)
tion = [2     ,  6    ,10   ,30   ,80   ,220  ,1200 ,8500  ,1e5];
sion = [1e-11 , 1e-9  ,5e-9 ,2e-8 ,3e-8 ,3e-8 ,2e-8 ,9e-9  ,3e-9] .* 1e-6; % m^3/s
si1s = exp(pchip(log(tion),log(sion),log(te)));
si1s = reshape(si1s,sl);

% ionisation de H (2s)
tion = [0.4   ,1      ,2     ,6     ,30      ,40     ,700   ,8000  ,20000  ,1e5  ];
sion = [1e-11 ,2.5e-9 ,2e-8  ,1e-7  ,1.7e-7  ,1.7e-7 ,1e-7  ,5e-8  ,3.7e-8 ,2e-8     ] .* 1e-6; % m^3/s
si2s = exp(pchip(log(tion),log(sion),log(te)));
si2s = reshape(si2s,sl);


% transition 1s -> 2s
tss  = [1    , 10      , 80    , 300  , 2000 , 20000];
sss  = [1e-11, 8.5e-9  , 3e-8  , 3e-8 , 2e-8 , 9e-9] .* 1e-6; % m^3/s;  
sss = exp(pchip(log(tss),log(sss),log(te)));
sss = reshape(sss,sl);

% desexitation
Ass = 1.1e9 .* ones(sl); % s^ -1

% echange de charge
tcx  = [1       ,3      ,10      ,23    ,120   ,700  ,2300 ,8000   ,9000   ,1e5];
scx  = [3.3e-8  ,3.5e-8 ,3.7e-8  ,4e-8  ,5e-8  ,7e-8 ,1e-7 ,1.6e-7 ,1.6e-7 ,2e-8] .* 1e-6;% m^3/s
scx  = exp(pchip(log(tcx),log(scx),log(te)));
scx = reshape(scx,sl);

% recombinaison
trec = [0.1   ,1       ,2     ,8     ,20       ,50       ,100     ,300     ,500   ,900   ,1500  ,1e5];
srec = [5e-13 ,1.7e-13 ,1e-13 ,4e-14 ,2.3e-14  ,1.2e-14  ,6e-15   ,1.5e-15 ,7e-16 ,3e-16 ,1e-16 ,1e-25] .* 1e-6;% m^3/s;
srec  = exp(pchip(log(trec),log(srec),log(te)));
srec = reshape(srec,sl);

% ionisation ion  +H  -> ion + e + H+
tii  = [0.1   ,300   ,1e5   ,1e6];
sii  = [1e-25   ,1e-11 ,3e-8  ,1e-9] .* 1e-6;% m^3/s
sii  = exp(pchip(log(tii),log(sii),log(ti)));
sii = reshape(sii,sl);
%il manque la dissociation de H2





