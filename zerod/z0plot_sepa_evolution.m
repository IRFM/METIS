% graph for separatrix 
if ~isfield(z0dinput.exp0d,'Rsepa')
    return
end
h = findobj(0,'type','figure','tag','z0geosepa');
if isempty(h)
    h=figure('tag','z0geosepa');
else
    figure(h);
end
clf
if isfield(z0dinput.exp0d,'Rsepa')
  plot(z0dinput.exp0d.Rsepa',(z0dinput.exp0d.Zsepa + z0dinput.geo.z0 * ones(1,size(z0dinput.exp0d.Zsepa,2)))');
end
axis('square')
axis('equal')
xlabel('R (m)');
xlabel('Z (m)');
try 
  title(sprintf('METIS : %s@%d/LCFS evolution',post.z0dinput.option.machine,post.z0dinput.option.shot));
catch
  title(sprintf('METIS : %s@%d/LCFS evolution',post.z0dinput.machine,post.z0dinput.shot));
end
edition2


