function y = psitbxtcv(shot,varargin)

%*PSITBXTCV	TCV LIUQE poloidal flux
% PSI = PSITBXTCV(SHOT[,T][,FORMAT][,SOURCE])

tw = NaN; form = '*0';
if shot == -1, from = 'FBTE';
else,          from = 'LIUQE'; end
for k = 1:nargin-1
 switch class(varargin{k})
  case 'char'
   varargin{k} = upper(varargin{k});
   switch varargin{k}
    case {'*0' '-0' '+0' '01' 'FS' 'JPHI'}, form = varargin{k};
    otherwise, from = varargin{k};
   end
  case 'double', tw = varargin{k};
 end
end

if ~isempty(mdsopen(shot)), return, end

switch form
 case {'*0' '-0','+0','01'}
  x = tdi('TCV_EQ("PSI",$1)',from);
  sip = tdi('TCV_EQ("I_P",$1)',from); sip = sign(mean(sip.data));
  if sip == 1, fip = '+0'; else, fip = '-0'; end
 case 'FS'
  x = tdi('\results::psitbx:as');
  rmag = mdsdata('\results::psitbx:rmag');
  zmag = mdsdata('\results::psitbx:zmag');
  pmag = mdsdata('\results::psitbx:psimag');
 case 'JPHI'
  x = tdi('TCV_EQ("J_TOR",$1)',from);
end
z = x.data;
t = x.dim{3};
x = x.dim(1:2);

if ~isnan(tw) & ~isempty(z)
 k = ifloor(round(t*1e6),round(tw(tw <= max(t))*1e6),1);
 k = k(k > 0); k(find(diff(k) == 0)) = [];
 t = t(k);
 z = z(:,:,k);
 if strcmp(form,'FS')
  rmag = rmag(k); zmag = zmag(k); pmag = pmag(k);
 end
end

switch form
 case '*0'
  y = psitbxpsi(z,psitbxgrid('C','G',x),t,fip);
 case '+0'
  y = psitbxpsi(sip*z,psitbxgrid('C','G',x),t,'+0');
 case '-0'
  y = psitbxpsi(-sip*z,psitbxgrid('C','G',x),t,'-0');
 case '01'
  y = psitbxp2p(psitbxpsi(sip*z,psitbxgrid('C','G',x),t,'+0'),'01');
 case 'FS'
  y = psitbxpsi(z,psitbxgrid('F','G',x),t,'FS',[rmag(:)';zmag(:)';pmag(:)']);
 case 'JPHI'
  y = psitbxfun(z,psitbxgrid('C','G',x),t);
end
