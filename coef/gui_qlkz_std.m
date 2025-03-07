% This function provide the GUI creation and handeling of QLK-NN 10D tool for METIS
% use :
%
% evalin('base','zuicreefunform(''gui_qlkANN10D'',''qlkparams'',1,0,''[chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkANN10D(qlkparams,post);'',''visible'');');
% 

function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkz_std(qlkzparams,post)

% parameters declaration
if nargin < 1
    chie_loop = declare_QLKZ_std_parameters;   
    return
end

% call of compute_energy_west
if nargout == 0
    z0qlkz_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d);
else
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0qlkz_std(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d);
end
