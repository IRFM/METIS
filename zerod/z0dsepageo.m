function [ z0dinput ] = z0dsepageo(z0dinput, Rsepa, Zsepa, t_efit)
%z0dsepageo Calculate the geo struct from R,Z sepa

vt             = ones(size(z0dinput.cons.temps));
teta           = linspace(0,2*pi,201);
ve             = ones(size(teta));
z0dinput.exp0d.Rsepa = NaN .* vt * ve;
z0dinput.exp0d.Zsepa = NaN .* vt * ve;
for k = 1:length(vt)
    indn   = find(t_efit >= z0dinput.cons.temps(k),1);
    if isempty(indn)
        indn  = length(t_efit);
    end
    % check if the separatrix is not nan, use a previous one if possible
    while any(isnan(Rsepa(indn,:))) || any(isnan(Zsepa(indn,:)))
        indn = indn - 1;
    end
    rr     = Rsepa(indn,:);
    zz     = Zsepa(indn,:);
    % first and last LCFS points are indentical in our EFIT
    rr(:,end)    = rr(:,1);
    zz(:,end)    = zz(:,1);
    rmin  = min(rr);
    rmax  = max(rr);
    ra    = max(0.3, 0.5 .* (rmin + rmax));
    a     = max(0.01,0.5 .* (rmax - rmin));
    % SECURITE TEMPS ETRANGE
    ra = max(a .* 1.01 , ra);
    zmin  = min(zz);
    zmax  = max(zz);
    za    = (zmin + zmax) ./ 2;
    b     = 0.5 .* (zmax - zmin);
    elon     = max(0.5,b ./ a);
    mask1 = (zz == zmax);
    mask2 = (zz == zmin);
    rzmax = max(rr .* mask1,[],2);
    rzmin = max(rr .* mask2,[],2);
    cl    = ra - rzmin;
    cu    = ra - rzmax;
    d     = (cl+cu) ./2 ./ a;
    z0dinput.geo.a(k)      = a;
    z0dinput.geo.R(k)      = ra;
    z0dinput.geo.z0(k)     = za;
    z0dinput.geo.K(k)      = elon;
    z0dinput.geo.d(k)      = d;
    
    % iso angle
    cl   = (rr - ra) + sqrt(-1) .* (zz -za);
    thl  = unwrap(angle(cl),[],2);
    thl(thl<0) = thl(thl<0) + 2 .* pi;
    rhol = abs(cl);
    [thl,indl] = sort(thl);
    rhol    =   rhol(indl);
    rhol = cat(2,rhol,rhol,rhol);
    thl = cat(2,thl -2.*pi,thl,thl+2.*pi);
    indnok = find(any(diff(thl,1,2)<=0,1));
    thl(:,indnok) =[];
    rhol(:,indnok)  = [];
    rho    = pchip(thl,rhol,teta);
    Rext   = ra + rho .* cos(teta);
    Zext   = za + rho .* sin(teta);
    z0dinput.exp0d.Rsepa(k,:) = Rext;
    z0dinput.exp0d.Zsepa(k,:) = Zext - za;
end

end

