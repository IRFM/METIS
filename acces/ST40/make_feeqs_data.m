[s,rep] = unix(' cat V2coilsAndPassiveStructures.JSON');
data    = jsondecode(rep);

figure;
plot(data.IVC.r,data.IVC.z,'g');
r_lim = cat(1,data.IVC.r(:),NaN);
z_lim = cat(1,data.IVC.z(:),NaN);

hold on 
plot(data.IVC.r,-data.IVC.z,'g');
r_lim = cat(1,r_lim,data.IVC.r(:),NaN);
z_lim = cat(1,z_lim,-data.IVC.z(:),NaN);
plot(data.first_wall.r,data.first_wall.z,'k');
r_lim = cat(1,r_lim,data.first_wall.r(:),NaN);
z_lim = cat(1,z_lim,data.first_wall.z(:),NaN);

plot(data.passive_rings.lower.r,data.passive_rings.lower.z,'b','linewidth',3);
r_lim = cat(1,r_lim,data.passive_rings.lower.r(:),NaN);
z_lim = cat(1,z_lim,-data.passive_rings.lower.z(:),NaN);
plot(data.passive_rings.upper.r,data.passive_rings.upper.z,'b','linewidth',3);
r_lim = cat(1,r_lim,data.passive_rings.upper.r(:),NaN);
z_lim = cat(1,z_lim,-data.passive_rings.upper.z(:),NaN);
plot(r_lim,z_lim,'.k');

% open file for limiter data
fid = fopen('ST80-HTS_limiter.txt','w');
fprintf(fid,'R (m)\tZ (m)\n');
for k=1:length(r_lim)
   fprintf(fid,'%g\t%g\n', r_lim(k),z_lim(k));
end
fclose(fid);

% open file for FEEQS data
fid = fopen('ST80-HTS_coils_fast_mode.txt','w');
fprintf(fid,'ST80-HTS\tRcentre\tZcentre\tDR_horizontal\tDZ_vertical\tnumber_of_turns\tPower_supply_cabling\n');
fprintf(fid,'Coil name\t(m)\t(m)\t(m)\t(m)\t(-)\t(-)\n');
%
cabling = 1;
noms = fieldnames(data.Solenoids);
for k=1:length(noms)
    sol = data.Solenoids.(noms{k});
    for l=1:length(sol.a1)
        r = [sol.a1(l),sol.a1(l),sol.a2(l),sol.a2(l),sol.a1(l)];
        z = [sol.b1(l),sol.b2(l),sol.b2(l),sol.b1(l),sol.b1(l)];
        switch k
            case 1
                plot(r,z,'r');
            case 2
                plot(r,z,'m');
            case 3
                plot(r,z,'c');                
        end
        Rc = (sol.a1(l) + sol.a2(l)) / 2;
        Zc = (sol.b1(l) + sol.b2(l)) / 2;
        dR = abs(sol.a1(l) - sol.a2(l));
        dZ = abs(sol.b1(l) - sol.b2(l));
        S  = dR*dZ;
        nbturns = ceil(S / 1e-3);
        lname = sprintf('%s_%d',noms{k},l);
        fprintf(fid,'%s\t%g\t%g\t%g\t%g\t%g\t%d\n',lname,Rc,Zc,dR,dZ,nbturns,cabling);
    end
    cabling = cabling + 1;
end

noms = fieldnames(data.PFCoils);
for k=1:length(noms)
    sol = data.PFCoils.(noms{k});
    for l=1:length(sol.a1)
        r = [sol.a1(l),sol.a1(l),sol.a2(l),sol.a2(l),sol.a1(l)];
        z = [sol.b1(l),sol.b2(l),sol.b2(l),sol.b1(l),sol.b1(l)];
        plot(r,z,'b');
        Rc = (sol.a1(l) + sol.a2(l)) / 2;
        Zc = (sol.b1(l) + sol.b2(l)) / 2;
        dR = abs(sol.a1(l) - sol.a2(l));
        dZ = abs(sol.b1(l) - sol.b2(l));
        S  = dR*dZ;
        nbturns = ceil(S / 1e-3);
        if length(sol.a1) > 1
            lname = sprintf('%s_%d',noms{k},l);
        else
            lname = noms{k};
        end
        fprintf(fid,'%s\t%g\t%g\t%g\t%g\t%g\t%d\n',lname,Rc,Zc,dR,dZ,nbturns,cabling);
    end
    cabling = cabling + 1;
end
fclose(fid);

axis('equal')
xlabel('R (m)');
ylabel('Z (m)');

