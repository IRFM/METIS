function tracemhd(post,param);
%
%
%
h = findobj(0,'type','figure','tag','mhd');
if isempty(h)
       h=figure('tag','mhd');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1],'defaultlinemarkersize',4,...
	'defaultaxesytick',[-2 -1 0 1 2],'defaultaxesyticklabel','unstable|unstable?|unknown|stable?|stable')

mode = abs(param.cons.mhd.stab.tor);

EW   = post.mhd.n1.EW;
t    = post.mhd.n1.t;
DEW  = post.mhd.n1.DEW;
NVI  = post.mhd.n1.NVI;

val  = zeros(size(EW))+2;

ver1      = EW > 1e-5 & DEW <= 1e-6 & NVI < 21;
if sum(ver1) > 0
  val(ver1) = -2*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW <= 1e-6 & NVI < 21;
if sum(ver1) > 0
  val(ver1) = 2*ones(size(val(ver1)));
end
ver1      = EW > 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI < 21;
if sum(ver1) > 0
  val(ver1) = -1*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI < 21;
if sum(ver1) > 0
  val(ver1) = 1*ones(size(val(ver1)));
end

subplot(3,1,1)

plot(t,val,'ro')
grid
text('units','normalized','position',[1.01 0.5],'string',sprintf('n=%d',mode(1)))
title(sprintf('%s, #%d, MHD Stability from MISHKA',param.from.machine,param.from.shot.num))   
EW   = post.mhd.n2.EW;
t    = post.mhd.n2.t;
DEW  = post.mhd.n2.DEW;
NVI  = post.mhd.n2.NVI;

val  = zeros(size(EW))+2;

ver1      = EW > 1e-5 & DEW <= 1e-6 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = -2*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW <= 1e-6 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = 2*ones(size(val(ver1)));
end
ver1      = EW > 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = -1*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = 1*ones(size(val(ver1)));
end

subplot(3,1,2)

plot(t,val,'ro')
grid
text('units','normalized','position',[1.01 0.5],'string',sprintf('n=%d',mode(2)))

EW   = post.mhd.n3.EW;
t    = post.mhd.n3.t;
DEW  = post.mhd.n3.DEW;
NVI  = post.mhd.n3.NVI;

val  = zeros(size(EW))+2;

ver1      = EW > 1e-5 & DEW <= 1e-6 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = -2*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW <= 1e-6 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = 2*ones(size(val(ver1)));
end
ver1      = EW > 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = -1*ones(size(val(ver1)));
end
ver1      = EW <= 1e-5 & DEW > 1e-6 & DEW <= 1e-5 & NVI <= 21;
if sum(ver1) > 0
  val(ver1) = 1*ones(size(val(ver1)));
end


subplot(3,1,3)

plot(t,val,'ro')
grid
xlabel('time (s)')
text('units','normalized','position',[1.01 0.5],'string',sprintf('n=%d',mode(3)))


