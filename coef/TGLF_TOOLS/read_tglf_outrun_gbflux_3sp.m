function outrun = read_tglf_outrun_gbflux_3sp(simname)

% This function reads fluxes from the TGLF output namelist file out.tglf.run as provided by the "simname" input string (given with full path and name) into the matlab structure "outrun"

%simname = '/Home/RK278477/codes/metis/coef/TGLF_TOOLS/runs/tglf_tmp_radial2';
clear outrun;
outrun = struct();
aa = load([simname, '/out.tglf.gbflux']);

%if ~isfield(outrun, 'ion2_flx')
%	outrun.ion2_flx = outrun.ion1_flx*0.0;
%end

outrun.Ge_gB = aa(1);
outrun.G1_gB = aa(2);
outrun.G2_gB = aa(3);
outrun.Qe_gB = aa(4);
outrun.Q1_gB = aa(5);
outrun.Q2_gB = aa(6);

end

