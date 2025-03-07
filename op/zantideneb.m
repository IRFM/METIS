function f = zantideneb(f)
%
% remplacement /usr/deneb/groupe/utilisateur par /usr/drfc/utilisateur
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 09/10/2003.
%
% liste des modifications : 
%

part = {};
r    = f;
l    = 0;
while(~isempty(r))
   [p,r]= strtok(r,filesep);
   part{end+1} = p;
  l = l+1;
end

if ~strcmp(part{2},'drfc') &  strcmp(part{1},'usr')
   s = part{3};
   if strcmp(s(1),'g')
      nf = '';
      for m=4:l-1  
       nf = [nf,part{m},filesep];
      end
      nf = [nf,part{l}];
      nf = [filesep,'usr',filesep,'drfc',filesep,nf];
      f  = nf;
   end 
end
