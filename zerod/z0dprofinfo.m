%  Z0DPROFINFO  courte description  
%------------------------------------------------------------------------------- 
% fichier :  z0dprofinfo.m  ->  z0dprofinfo ,    
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   profinfo = z0dprofinfo 
%  
% entrees :  
%  
%  
% sorties :  
%   profinfo  = 
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  3.1  du  18/11/2005  
%  
%@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
function profinfo = z0dprofinfo

profinfo.xli = 'radial normalized coordinate (Lao coordinate ~ r/a, normalized to 1 at the edge).';
profinfo.Rsepa='radial coordinate of the LCMS points (m)';
profinfo.Zsepa='vertical coordinate of the LCMS points (m)';
profinfo.nbishape_el = 'NBI power deposition on electrons shape; reserved to an internal used only';
profinfo.nbishape_ion = 'NBI power deposition on ions shape; reserved to an internal used only';
profinfo.jnbishape= 'NBICD current density shape; reserved to an internal used only';
profinfo.pitch = 'NBI pitch angle profile = cos(beam,B)';
profinfo.tep= 'electron temperature  (eV)';
profinfo.tip= 'ion temperature  (eV)';
profinfo.jli= 'averaged current density, equivalent at jmoy in CRONOS (A/m^2)';
profinfo.jeff= '<J.B>/B0 (A/m^2)';
profinfo.qjli= 'safety factor';
profinfo.jboot= ' bootstrap current density (A/m^2)';
profinfo.eta= 'neoclassical resistivity (ohm * m)';
profinfo.jnbicd= 'parallel current density source  due to NBICD (j = <J.B>/Bo, A/m^2)';
profinfo.jlh= 'parallel current density source  due to LHCD (j = <J.B>/Bo, A/m^2)';
profinfo.jeccd= 'parallel current density source  due to ECCD (j = <J.B>/Bo, A/m^2)';
profinfo.jfwcd= 'parallel current density source  due to FWCD (j = <J.B>/Bo, A/m^2)';
profinfo.jfus= 'fast alpha bootstrap current density (A/m^2)';
profinfo.jrun= 'average runaway  current density (A/m^2)';
profinfo.plh= 'power density  due to LH (W/m^3)';
profinfo.pnbi= 'power density  due to NBI (W/m^3)';
profinfo.pnbi_ion= 'power density  due to NBI, coupled to ions (W/m^3)';
profinfo.pecrh= 'power density  due to ECRH (W/m^3)';
profinfo.pfweh= 'power density  due to FWEH (W/m^3)';
profinfo.picrh= 'power density  due to ICRH minority scheme (W/m^3)';
profinfo.picrh_ion= 'power density  due to ICRH minority scheme, coupled to ions (W/m^3)';
profinfo.pfus= 'power density heating the plasma due to fusion reactions (W/m^3)';
profinfo.pfus_ion= 'power density heating the plasma due to fusion reactions, coupled to ions (W/m^3)';
profinfo.pbrem= 'bremsstrahlung power  sink (W/m^3)';
profinfo.prad= 'radiated power  sink  (W/m^3)';
profinfo.pcyclo= 'cyclotron radiation power  source (W/m^3)';
profinfo.pohm= 'ohmic power deposition  (W/m^3)';
profinfo.nep= 'electron density  (m^-3)';
profinfo.nip= 'ions density  (m^-3)';
profinfo.vpr= 'volume element (m^3, int(vpr,x= 0..1) = plasma volume)';
profinfo.vpr_tor= 'volume element, reserved for internal used (m^2, rhomax * int(vpr,x= 0..1) = plasma volume) ';
profinfo.spr= 'surface element (m^2, int(spr,x= 0..1) = plasma poloidal section surface) ';
profinfo.grho2r2= '<|gradient(rho)|^2/R^2> (see CRONOS technical document)';
profinfo.r2i= '<1/R^2> (see CRONOS technical document)';
profinfo.ri= '<1/R> (see CRONOS technical document)';
profinfo.C2= 'C2 geometrical coefficient (see CRONOS technical document)';
profinfo.C3= 'C3 geometrical coefficient (see CRONOS technical document)';
profinfo.grho= '<|gradient(rho)|> (see CRONOS technical document)';
profinfo.grho2= '<|gradient(rho)|^2> (see CRONOS technical document)';
profinfo.ej= 'ohmic power density (W/m^3)';
profinfo.kx= 'flux surface elongation ';
profinfo.dx= 'flux surface geometrical triangularity ';
profinfo.Raxe= ' of major radius of the centre of each flux surface (m)';
profinfo.epsi= ' of aspect ratio (m)';
profinfo.rmx= ' of average radius of each flux surface: rmx   = sqrt(phi ./ pi ./ b0) (m)';
profinfo.bpol= 'average poloidal magnetic field  (T)';
profinfo.fdia= 'diamagnetic function  (T*m)';
profinfo.psi= 'poloidal flux  (Wb)';
profinfo.dpsidt= 'time derivative of poloidal flux  (V)';
profinfo.epar= 'parallel electric field  (V/m) ';
profinfo.zeff= 'effective charge';
profinfo.n1p= 'density of HDT ions (m^-3) or H+D if option.gza = 1';
profinfo.nhep= 'density of helium (m^-3) or helium 3 density if option.gza = 5';
profinfo.nzp= 'density of main impurity (m^-3)';
profinfo.xieshape= 'heat transport coefficient (Ke) shape without ITB effect ';
profinfo.xieshape_itb= 'heat transport coefficient (Ke) shape with ITB effect';
profinfo.source_ion= 'total heat power density coupled to ions (W/m^3)';
profinfo.source_el= 'total heat power density coupled to electrons (W/m^3)';
profinfo.jni= 'total current density source (A/m^2)';
profinfo.ftrap= 'effective trap fraction profile ';
profinfo.ptot= 'total pressure profile (Pa)';
profinfo.jfusshape= 'alpha particles bootstrap current density shape; reserved to an internal used only ';
profinfo.salf= 'alpha particles source (m^-3 s^-1)';
profinfo.palf= 'alpha power source (W/m^-3)';
profinfo.fprad= 'line radiative power shape; reserved to an internal used only';
profinfo.temps='vectors of time associate to the s  (only time at which the s are computed are stored)';
profinfo.xie= 'electron diffusivity estimation (m^2/s)';
profinfo.xii= 'ion diffusivity estimation (m^2/s)';
profinfo.n0 = 'neutral density coming from edge, hot neutral (m^-3)';
profinfo.s0 = 'ionisation sources coming from edge, hot neutral (m^-3/s)';
profinfo.n0m = 'neutral density coming from edge, cold neutral (m^-3)';
profinfo.s0m = 'ionisation sources coming from edge, cold neutral (m^-3/s)';
profinfo.ware = 'Ware pinch estimation (m/s)';
profinfo.dn= 'density diffusivity estimation (m^2/s)';
profinfo.vn= 'anormal density convection velocity estimation (m^2/s)';
profinfo.omega ='plasma solid rotation frequency in toroidal direction (rad/s) { sum(nk*mk*<Vk,phi * R>) /  sum(nk*mk)}';
profinfo.vtheta = 'fluid velocity, theta component, at R = Rmax of each flux surface, for main impurity (m/s)';
profinfo.utheta = 'neoclassical poloidal rotation speed ,for main impurity  (<V_k . theta> / <B . theta> , m/s/T)';
profinfo.vtor   ='toroidal rotation speed (m/s), at R = Rmax of each flux surface, for  main impurity ';
profinfo.er     = 'neoclassical radial electric field (V/m) = Er / gradient(rho)';
profinfo.spellet = 'equivalent continue source due to pellet injection (m^-3/s)';
profinfo.ge      = 'electrons flux  (m^-2/s)';
profinfo.pioniz    = 'loss power due to cold neutral ionization and charge exchange between ions and cold neutrals(W m^-3)';
profinfo.pioniz_i  = 'loss power due charge exchange between ions and cold neutrals(W m^-3)';
profinfo.phi  = 'toroidal flux (Wb)';
profinfo.dphidx  = 'toroidal flux space derivative (Wb)';
profinfo.nbinesource = 'electron source due to fast neutral ionisation (Number s^-1 m^-3)';  
profinfo.web         = 'rotation shear (NClass definition)';
profinfo.qe          = 'electron heat flux (W)';
profinfo.qi          = 'ion heat flux (W)';
profinfo.qei         = 'electron to ion heat flux (W)';
profinfo.rot_nbi     = 'toroidal torque source from NBI (N m^-2)';
profinfo.rot_lh      = 'toroidal torque source from LHCD and ECCD (N m^-2)';
profinfo.rot_n0      = 'toroidal torque sink from edge neutral friction (N m^-2)';
profinfo.frot        = 'toroidal rotation moment flux (N m^-1)';
profinfo.drot        = 'toroidal rotation diffusion coefficient (m^2/s)';
profinfo.vrot        = 'toroidal rotation pinch coefficient (m/s)';
profinfo.rtor        = 'toroidal rotation moment density (kg m^-1 s^-1)';
profinfo.nwp         = sprintf('%s\n%s\n%s','density of W when tungsten effects are taking into account (m^-3);', ...
                               'if option.Sn_fraction > 0, this profile become the sum of the density of W and Sn:', ...
                               'n_W = (1 - option.Sn_fraction) * nwp & n_Sn =  option.Sn_fraction * nwp.');

% profile for coupling to FREEBIE
profinfo.jgs         = 'averaged current density, raw computation for FREEBIE (A/m^2)';
profinfo.df2dpsi     = 'd(F^2)/dPsi for FREEBIE (4*pi*T)';
profinfo.dptotdpsi   = 'dPtotal/dPsi for FREEBIE (2*pi*A/m^3)';

% hard x ray constraint
profinfo.fxdurplh    = 'Hard X ray inversed profile';

profinfo = orderfields(profinfo);

