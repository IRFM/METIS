function y = subsref(x,s)

% PSITBXFUN/SUBSREF	Field reference of PsiTbx-Function objects
y = builtin('subsref',x,s);
