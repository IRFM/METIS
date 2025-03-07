% script de test des section efficaces.
Zeff = 2;
te = 10e3;
ti = 15e3;
ne = [1e12,1e13,1e14,1e15] .* 1e6;
vn  = ones(size(ne));
A  = 1;

E  = logspace(2,4)' .* 1e3;
ve = ones( size(E));
nhe = 0;
nhigh = 0;
zhigh = 56;
zlow  = 6;
nlow  = ne .* (Zeff -1) ./ (zlow .^ 2 -zlow);
ni    = nlow + (ne - zlow .* nlow);
sv0 = z0nbistop(A,E*vn,ve*ne,te,nhe,ve*nlow,nhigh,zlow,zhigh);
sv1 = z0signbi(te,ve * ne,ve*ni,E*vn,A)./ (ve*ne);
sv2  = z0sinbad(A,E*vn,ve*ne,te,nhe,ve*nlow,nhigh,zlow,zhigh)./ (ve*ne);
eref = [1e2        ,5e2           ,1e3        ,2e3        ,5e3       ,1e4] .*1e3;
s1   = [2.3e-16    ,7.5e-17       ,4.4e-17      ,2.8e-17    ,1.3e-17   ,8e-18].*1e-4;    
s2   = [2.5e-16    ,8e-17         ,5e-17      ,3e-17      ,1.5e-17   ,9e-18].*1e-4;    
s3   = [2.9e-16    ,1e-16         ,6.5e-17    ,3.5e-17    ,1.9e-17   ,1e-17].*1e-4;    
s4   = [3.7e-16    ,1.5e-16       ,9e-17      ,5e-17      ,2.7e-17   ,1.5e-17].*1e-4;    

h = findobj(0,'type','figure','tag','svnbi');
if isempty(h)
       h=figure('tag','svnbi');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

loglog(E,sv0,E,mean(sv1,2),'-.',E,sv2,':');
xlabel('E0 (eV)');
ylabel('SV(E,...) (m^2)');
legend('1e18','1e19','1e20','1e21','VB mfile','asymptotic')
hold on
loglog(eref,s1,'o',eref,s2,'o',eref,s3,'o',eref,s4,'o')
edition2

h = findobj(0,'type','figure','tag','svnbi2');
if isempty(h)
       h=figure('tag','svnbi2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

semilogx(E,1./(sv0.*(ve*ne)),E,1./(sv1.*(ve*ne)),'-.',E,1./(sv2.*(ve*ne)),':');
xlabel('E0 (eV)');
ylabel('L (m)');


