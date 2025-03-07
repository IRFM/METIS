function y = subsref(dcd,s)

% PSITBXDCD/SUBSREF

if length(s) == 0
 y = dcd;
elseif strcmp(s(1).type,'()')
 rd   = builtin('subsref',dcd.rd,  s);
 zd   = builtin('subsref',dcd.zd,  s);
 phid = builtin('subsref',dcd.phid,s);
 pvd  = builtin('subsref',dcd.pvd, s);
 tvd  = builtin('subsref',dcd.tvd, s);
 y = subsref(psitbxdcd(rd,zd,phid,pvd,tvd,dcd.nd,dcd.t),s(2:end));
else
 y = builtin('subsref',dcd,s);
end
