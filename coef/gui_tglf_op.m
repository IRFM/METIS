% This function provide the GUI creation and handeling of TGLF for METIS
% use :
%
% evalin('base','zuicreefunform(''gui_qlkANN10D'',''qlkparams'',1,0,''[chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_qlkANN10D(qlkparams,post);'',''visible'');');
% 

function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_tglf_op(tglfparams,post)

% parameters declaration
if nargin < 1
    chie_loop = declare_TGLF_op_parameters;   
    return
end

% call of compute tglf
%if nargout == 0
%    z0tglf(qlkzparams,post.z0dinput.option,post.zerod,post.profil0d);
%else
%    [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = run_tglf(tglfparams,post.z0dinput.option,post.zerod,post.profil0d);
%end
