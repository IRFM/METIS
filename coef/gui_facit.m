function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] =  gui_facit(facitparams,post)

% parameters declaration
if nargin < 1
    chie_loop = declare_FACIT_parameters;   
    return
end

