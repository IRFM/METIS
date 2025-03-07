% script to test equilibrium in metis versu chease for ITER cocos 11
% lecture IDS equilibrium for chease
load('equi_ids_metis_astra.mat')

%
%[output,idss_list] = litidss('equilibrium',55,11);
chease = equilibrium.time_slice{1};
% load of matfile containing IDS for metis
%ref = load('testexportimas.mat');
%ref = ref.exported_data.equilibrium.time_slice{300};
ref = equilibrium2.time_slice{100};

% global quantities
noms = fieldnames(ref.global_quantities);
for k=1:length(noms)
   if isempty(strfind(noms{k},'_error_'))
      var1 = ref.global_quantities.(noms{k});
      var2 = chease.global_quantities.(noms{k});
      if iscell(var1)
	keyboard
      elseif isstruct(var1)
	nn = fieldnames(var1);
	for l=1:length(nn)
	  if isempty(strfind(nn{l},'_error_'))
	    fprintf('%s.%s METIS = %g & %s.%s ASTRA = %g\n', noms{k},nn{l},var1.(nn{l}),noms{k},nn{l},var2.(nn{l}));   
	  end
	end
      else
	fprintf('%s METIS = %g & %s ASTRA = %g\n', noms{k},var1,noms{k},var2);
      end
   end
end

% x_points
for k=1:max(length(ref.boundary.x_point),length(chease.boundary.x_point))
  if k <= length(ref.boundary.x_point)
     if length(ref.boundary.x_point) > 1
        fprintf('METIS x_point @ (%g,%g) &',ref.boundary.x_point{k}.r,ref.boundary.x_point{k}.z);
     elseif (length(ref.boundary.x_point) == k) && isfield(ref.boundary.x_point,'r')
      fprintf('METIS x_point @ (%g,%g) &',ref.boundary.x_point.r,ref.boundary.x_point.z);
     end
  end
  if k <= length(chease.boundary.x_point)
    if length(chease.boundary.x_point) > 1
      fprintf('ASTRA x_point @ (%g,%g)',chease.boundary.x_point{k}.r,chease.boundary.x_point{k}.z);  
   
    elseif (length(chease.boundary.x_point) == k) && isfield(chease.boundary.x_point,'r')
      fprintf('ASTRA x_point @ (%g,%g) &',chease.boundary.x_point.r,chease.boundary.x_point.z);
    end
  end
  fprintf('\n');
end
figure;
plot(ref.boundary.outline.r,ref.boundary.outline.z,'b',chease.boundary.lcfs.r,chease.boundary.lcfs.z,'r');
try 
  load iterwall
   hold on;
   plot(R,Z,'k')
end
title('Boundary : blue = METIS & red = ASTRA');
xlabel('R')
ylabel('Z')
legend('METIS','ASTRA');
ref.profiles_1d.psi = ref.profiles_1d.psi - ref.profiles_1d.psi(1);
noms = fieldnames(ref.profiles_1d);
for k=1:length(noms)
   if isempty(strfind(noms{k},'_error_'))
       figure
       leg = {};
       if length(ref.profiles_1d.volume) == length(ref.profiles_1d.(noms{k}))
	  plot(ref.profiles_1d.volume,ref.profiles_1d.(noms{k}),'b');
	  hold on
	  leg{1} = 'METIS';
       end  
       if length(chease.profiles_1d.volume) == length(chease.profiles_1d.(noms{k}))
	  plot(chease.profiles_1d.volume,chease.profiles_1d.(noms{k}),'r');
	  leg{end + 1} = 'ASTRA';
       end
       ylabel(noms{k});
       xlabel('Volume (m^3)')
      legend(leg);
   end
end
drawnow
% first METIS grid
% iso psi
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.psi,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.psi,101,'color','r')
%axis('equal');
title('PSI: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso PHI
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.phi,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.phi,101,'color','r')
%axis('equal');
title('PHI: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso j_tor
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.j_tor,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.j_tor,101,'color','r')
%axis('equal');
title('j_tor: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso j_parallel
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.j_parallel,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.j_parallel,101,'color','r')
%axis('equal');
title('j_parallel: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso br
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_r,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_r,101,'color','r')
%axis('equal');
title('Br: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso bz
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_z,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_z,101,'color','r')
%axis('equal');
title('Bz: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso btor
figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_tor,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_tor,101,'color','r')
%axis('equal');
title('Btor: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% vector

figure
contour(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.psi,101,'color','b','linestyle',':')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.psi,101,'color','r','linestyle',':')
quiver(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_r,ref.profiles_2d{1}.b_field_z,'color','b');
quiver(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_r,chease.profiles_2d{1}.b_field_z,'color','r');
title('Bpoloidal: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');



% value on boundary
bphi_ref = griddata(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_tor,ref.boundary.outline.r,ref.boundary.outline.z,'cubic');
br_ref = griddata(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_r,ref.boundary.outline.r,ref.boundary.outline.z,'cubic');
bz_ref = griddata(ref.profiles_2d{1}.r,ref.profiles_2d{1}.z,ref.profiles_2d{1}.b_field_z,ref.boundary.outline.r,ref.boundary.outline.z,'cubic');
bphi_chease = griddata(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_tor,chease.boundary.outline.r,chease.boundary.outline.z,'cubic');
br_chease = griddata(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_r,chease.boundary.outline.r,chease.boundary.outline.z,'cubic');
bz_chease = griddata(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_z,chease.boundary.outline.r,chease.boundary.outline.z,'cubic');

figure
plot(linspace(0,1,length(bphi_ref)),bphi_ref,'b')
hold on 
plot(linspace(0,1,length(bphi_chease)),bphi_chease,'r')
title('Btor @ LCFS : blue = METIS & red = ASTRA');

figure
plot(linspace(0,1,length(br_ref)),br_ref,'b')
hold on 
plot(linspace(0,1,length(br_chease)),br_chease,'r')
title('Br @ LCFS : blue = METIS & red = ASTRA');

figure
plot(linspace(0,1,length(bz_ref)),bz_ref,'b')
hold on 
plot(linspace(0,1,length(bz_chease)),bz_chease,'r')
title('Bz @ LCFS : blue = METIS & red = ASTRA');


drawnow
% second METIS grid
% iso psi
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.psi,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.psi,101,'color','r')
%axis('equal');
title('PSI: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso PHI
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.phi,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.phi,101,'color','r')
%axis('equal');
title('PHI: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso j_tor
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.j_tor,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.j_tor,101,'color','r')
%axis('equal');
title('j_tor: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso j_parallel
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.j_parallel,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.j_parallel,101,'color','r')
%axis('equal');
title('j_parallel: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso br
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.b_field_r,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_r,101,'color','r')
%axis('equal');
title('Br: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso bz
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.b_field_z,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_z,101,'color','r')
%axis('equal');
title('Bz: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% iso btor
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.b_field_tor,101,'color','b')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_tor,101,'color','r')
%axis('equal');
title('Btor: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% vector

figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.psi,101,'color','b','linestyle',':')
hold on
contour(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.psi,101,'color','r','linestyle',':')
quiver(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.b_field_r,ref.profiles_2d{2}.b_field_z,'color','b');
quiver(chease.profiles_2d{1}.r,chease.profiles_2d{1}.z,chease.profiles_2d{1}.b_field_r,chease.profiles_2d{1}.b_field_z,'color','r');
title('Bpoloidal: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');

% vector on the same grid
FR = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.b_field_r(:),'natural','linear');
br_chease = FR(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z);
FZ = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.b_field_z(:),'natural','linear');
bz_chease = FZ(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z);
FPSI = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.psi(:),'natural','linear');
psi_chease = FPSI(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z);
figure
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.psi,101,'color','b','linestyle',':')
hold on
contour(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,psi_chease,101,'color','r','linestyle',':')
quiver(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,ref.profiles_2d{2}.b_field_r,ref.profiles_2d{2}.b_field_z,'color','b');
quiver(ref.profiles_2d{2}.r,ref.profiles_2d{2}.z,br_chease,bz_chease,'color','r');
title('Bpoloidal: blue = METIS & red = ASTRA');
xlabel('R (m)');
ylabel('Z (m)');



% value on boundary
F = scatteredInterpolant(ref.profiles_2d{2}.r(:),ref.profiles_2d{2}.z(:),ref.profiles_2d{2}.b_field_tor(:),'natural','linear');
bphi_ref = F(ref.boundary.outline.r,ref.boundary.outline.z);
F = scatteredInterpolant(ref.profiles_2d{2}.r(:),ref.profiles_2d{2}.z(:),ref.profiles_2d{2}.b_field_r(:),'natural','linear');
br_ref = F(ref.boundary.outline.r,ref.boundary.outline.z);
F = scatteredInterpolant(ref.profiles_2d{2}.r(:),ref.profiles_2d{2}.z(:),ref.profiles_2d{2}.b_field_z(:),'natural','linear');
bz_ref = F(ref.boundary.outline.r,ref.boundary.outline.z);
F = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.b_field_tor(:),'natural','linear');
bphi_chease = F(chease.boundary.outline.r,chease.boundary.outline.z);
F = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.b_field_r(:),'natural','linear');
br_chease = F(chease.boundary.outline.r,chease.boundary.outline.z);
 F = scatteredInterpolant(chease.profiles_2d{1}.r(:),chease.profiles_2d{1}.z(:),chease.profiles_2d{1}.b_field_z(:),'natural','linear');
bz_chease = F(chease.boundary.outline.r,chease.boundary.outline.z);

figure
plot(linspace(0,1,length(bphi_ref)),bphi_ref,'b')
hold on 
plot(linspace(0,1,length(bphi_chease)),bphi_chease,'r')
title('Btor @ LCFS : blue = METIS & red = ASTRA');

figure
plot(linspace(0,1,length(br_ref)),br_ref,'b')
hold on 
plot(linspace(0,1,length(br_chease)),br_chease,'r')
title('Br @ LCFS : blue = METIS & red = ASTRA');

figure
plot(linspace(0,1,length(bz_ref)),bz_ref,'b')
hold on 
plot(linspace(0,1,length(bz_chease)),bz_chease,'r')
title('Bz @ LCFS : blue = METIS & red = ASTRA');

drawnow

