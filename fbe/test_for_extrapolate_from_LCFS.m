% test for extrapolate_from_LCFS
eq = imas_west_get(55550,'equilibrium');
itime = 25;
R_LCFS = eq.boundary.outline.r{itime};
Z_LCFS = eq.boundary.outline.z{itime};
%FPSI   = scatteredInterpolant(eq.RNodes(:),eq.ZNodes(:),eq.Psi_val(itime,:)','natural','linear');
FBR    = scatteredInterpolant(eq.RNodes(:),eq.ZNodes(:),eq.b_r_val(itime,:)','natural','linear');
FBZ    = scatteredInterpolant(eq.RNodes(:),eq.ZNodes(:),eq.b_z_val(itime,:)','natural','linear');
BR_LCFS = FBR(R_LCFS,Z_LCFS);
BZ_LCFS = FBZ(R_LCFS,Z_LCFS);
PSI_LCFS = ones(size(R_LCFS)) * eq.global_quantities.psi_boundary(itime) ./ 2 ./ pi;
%
[FPSI,FBR,FBZ] = extrapolate_from_LCFS(R_LCFS,Z_LCFS,PSI_LCFS,BR_LCFS,BZ_LCFS,1);

figure;
hold on
plot(R_LCFS,Z_LCFS,'k');
psi_list = linspace( eq.global_quantities.psi_boundary(itime), 3 .* eq.global_quantities.psi_boundary(itime) - 2 .* eq.global_quantities.psi_axis(itime),31);
contour(eq.RInterp,eq.ZInterp,squeeze(eq.Psi_interp(itime,:,:)),psi_list,'color','r');
contour(eq.RInterp,eq.ZInterp,2 .* pi .* FPSI(eq.RInterp,eq.ZInterp),psi_list,'color','b');




