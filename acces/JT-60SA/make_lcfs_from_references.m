function z0dinput = make_lcfs_from_references(z0dinput,lcfs_param)

sepa_option = [];
z0dinput.exp0d.Rsepa = [];
z0dinput.exp0d.Zsepa = [];
%
for k = 1:length(z0dinput.cons.temps)
    noms = fieldnames(lcfs_param);
    for l=1:length(noms)
        if length(lcfs_param.time) == length(lcfs_param.(noms{l}))
            sepa_option.(noms{l}) = interp1lcfs(lcfs_param.time,lcfs_param.(noms{l}),z0dinput.cons.temps(k)); 
        else
            sepa_option.(noms{l}) = lcfs_param.(noms{l});
        end
    end    
    sepa_option.b0    = 1;
    sepa_option.delta = 1;
    sepa = z0dsepanew2([],sepa_option);
    
    rb0 = z0dinput.geo.R(k) .* z0dinput.geo.b0(k);
    z0dinput.geo.K(k) = sepa.e1;
    z0dinput.geo.d(k) = (sepa.trl + sepa.trh)/2;
    z0dinput.geo.R(k) = sepa.r0;
    z0dinput.geo.a(k) = sepa.a;
    z0dinput.geo.z0(k) = sepa.z0;
    z0dinput.exp0d.Rsepa(k,:) = sepa.R(1,:);
    z0dinput.exp0d.Zsepa(k,:) = sepa.Z(1,:) - sepa.z0;
    % keep magnetic rigidity
    z0dinput.geo.b0(k) = rb0 ./ z0dinput.geo.R(k);
end


t  = asin(max(0,min(1,z0dinput.geo.d)));
u  = linspace(0,2.*pi,201);
vu = ones(size(u));
Rtest  = z0dinput.geo.R *vu + (z0dinput.geo.a * vu) .* cos(ones(size(z0dinput.geo.R,1),1) * u + t * sin(u));
Ztest = (z0dinput.geo.a .* z0dinput.geo.K) * sin(u);

h = findobj(0,'type','figure','tag','z0geosepa');
if isempty(h)
h=figure('tag','z0geosepa');
else
figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
zplotprof(gca,z0dinput.cons.temps,Rtest,Ztest+z0dinput.geo.z0 * vu,'color','r','marker','o','linestyle','none');
zplotprof(gca,z0dinput.cons.temps,z0dinput.exp0d.Rsepa,z0dinput.exp0d.Zsepa+z0dinput.geo.z0 * vu,'color','b','marker','none','linestyle','-');
axis('square')
axis('equal')




function out = interp1lcfs(tin,in,tout)

out = interp1(tin,in,tout,'linear','extrap');
if any(tout > max(tin))
    out(tout > max(tin)) = in(end);
end
if any(tout< min(tin))
    out(tout < min(tin)) = in(1);
end





