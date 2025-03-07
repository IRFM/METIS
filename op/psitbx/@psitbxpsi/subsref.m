function y = subsref(x,s)

% PSITBXPSI/SUBSREF	Field reference of PsiTbx-Poloidal-Flux objects
try
 y = builtin('subsref',x,s);
catch
 y = builtin('subsref',x.psitbxfun,s);
end
