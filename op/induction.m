function [lieq,Bpol]=induction(BPR,BPZ,dVhel,Rhel,Zhel,rhohel)
%
%[lieq,Bpol]=induction(BPR,BPZ,dVhel,Rhel,Zhel,rhohel);
%
  [n,m]        = size(Rhel);
  thet         = atan2(-BPR,BPZ);
  thet(thet<0) = thet(thet<0)+2*pi;
  thet(:,m)    = zeros(size(thet(:,m)))+2*pi;
  thet(:,1)    = zeros(size(thet(:,m)));  
  dVdx         = dVhel*max(rhohel);
  xhel         = rhohel./max(rhohel);
  modB         = sqrt(BPR.^2+BPZ.^2);
  for k=2:n
  
    Bpol(k)    = trapz(thet(k,:),modB(k,:));
  
  end  
  Bpol(1)      = 0;      
  lieq         = trapz(xhel,Bpol.^2 .* dVdx)./Bpol(end).^2./trapz(xhel,dVdx);
