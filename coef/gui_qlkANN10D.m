% This function provide the GUI creation and handeling of QLK-NN 10D tool for METIS
% use :
%
% evalin('base','zuicreefunform(''gui_qlkANN10D'',''qlkparams'',1,0,''[chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkANN10D(qlkparams,post);'',''visible'');');
% 

function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkANN10D(qlkparams,post)

% test if QLK-NN 10D mexfile exsit
if exist('qlknn_mex') ~= 3
    errordlg(sprintf('Qualikiz neural network 10D do not seem be installed on this version of METIS (qlknn_mex mexfile is not in the path):\n Please, consider to install it using the mfile install_QLKNN10D.'),'Error: QLK-NN 10D not available');
    chie_loop = [];
    return
end

% parameters declaration
if nargin < 1
    chie_loop = declare_QLKNN10D_parameters;   
    return
end

% call of compute_energy_west
if nargout == 0
    z0qlkANN_10D(qlkparams,post.z0dinput.option,post.zerod,post.profil0d);
else
    [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0qlkANN_10D(qlkparams,post.z0dinput.option,post.zerod,post.profil0d);
end
