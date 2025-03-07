function [lieq,Bpol,surf]=inductionJET(BPR,BPZ,Rhel,Zhel,equil,indR,indZ,mask)
%
%[lieq,Bpol,surf]=inductionJET(BPR,BPZ,Rhel,Zhel,equil,indR,indZ,mask);
%
% equil = 1 -> equilibre efit
% equil = 0 -> equilibre HELENA
%
if nargin < 6
  indR = [];
  indZ = [];
  mask = [];
end

if equil == 1
  warning off
  [n,m]       = size(Rhel);
  modB        = sqrt(BPR.^2+BPZ.^2);
  sur              = zeros(size(Rhel));
  nBpol            = zeros(size(Rhel));
  nk               = 1;
  for k=1:n-1
    for l=1:m-1
      if mask(k,l)==1
        R11        = Rhel(k+1,l+1);
        R10        = Rhel(k+1,l);
        R01        = Rhel(k,l+1);
        R00        = Rhel(k,l);
        Z11        = Zhel(k+1,l+1);
        Z10        = Zhel(k+1,l);
        Z01        = Zhel(k,l+1);
        Z00        = Zhel(k,l);
      
        sur(k,l)   = polyarea([R11 R10 R00 R01],[Z11 Z10 Z00 Z01]);
        B00        = modB(k,l);
        B01        = modB(k,l+1);
        B10        = modB(k+1,l);      
        B11        = modB(k+1,l+1);
      
        nBpol(k,l) = ((B11+B01)/2+(B00+B10)/2)/2;
        indf       = find(indR==l&indZ==k);
        if ~isempty(indf)
         
          frontB(nk)   = ((B11+B01)/2+(B00+B10)/2)/2;
          frontsur(nk) = sur(k,l);
          nk           = nk+1;
         
        end
      end        
    end     
  end  
  Bpol2           = sum(sur.*nBpol.^2);
  Bpol            = sum(sur.*nBpol);
  surfront        = sum(frontsur);  
  Bpol2front      = sum(frontsur.*(frontB.^2));
  surf            = sum(sur,2);
  ind             = surf~=0;
  surf            = surf(ind);
  Bpol2           = Bpol2(ind);      
  lieq            = surfront*sum(Bpol2)./Bpol2front./sum(surf);
  surf            = sum(surf);
  warning on
end


if equil == 0
  warning off
  [n,m]            = size(Rhel);
  modB             = sqrt(BPR.^2+BPZ.^2);
  sur              = zeros(size(Rhel));
  nBpol            = zeros(size(Rhel));
  nk               = 1;
  for k=1:n-1
    for l=1:m-1
        R11        = Rhel(k+1,l+1);
        R10        = Rhel(k+1,l);
        R01        = Rhel(k,l+1);
        R00        = Rhel(k,l);
        Z11        = Zhel(k+1,l+1);
        Z10        = Zhel(k+1,l);
        Z01        = Zhel(k,l+1);
        Z00        = Zhel(k,l);
      
        sur(k,l)   = polyarea([R11 R10 R00 R01],[Z11 Z10 Z00 Z01]);
        B00        = modB(k,l);
        B01        = modB(k,l+1);
        B10        = modB(k+1,l);      
        B11        = modB(k+1,l+1);
      
        nBpol(k,l) = ((B11+B01)/2+(B00+B10)/2)/2;
       
    end     
  end  
  Bpol2           = sum(sur.*nBpol.^2);
  Bpol            = sum(sur.*nBpol);
  surfront        = sum(sur(n-1,:));  
  Bpol2front      = sum(sur(n-1,:).*(nBpol(n-1,:).^2));
  surf            = sum(sur,2);
      
  lieq            = surfront*sum(Bpol2)./Bpol2front./sum(surf);
  surf            = sum(surf);
  warning on
end
