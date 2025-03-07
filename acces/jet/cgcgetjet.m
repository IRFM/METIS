
% liste{1} = 'PPF/@shot/EFIT/RMAG'
%                  | imperratif
% numchoc = 53521
% user,passwd -> le votre

% Update record:
% YYMMDD Who Comments
% 080521 AJC Removed disconnect and retry, complete waste of time.

function data=cgcgetjet(numchoc,liste,void1,void2)

if nargin < 2
	error('il faut 2 arguments en entree ...');
end
if isempty(numchoc)
	error('il faut donner le numero de choc')
end
if isempty(liste)
	error('pas de donnee demandee')
end

data = [];
data.void = [];
connect = 1;
while ~isempty(liste)
	% nouveau signal
	nomc  =liste{end};
	nomc  = strrep(nomc,'@shot/','');
	liste(end) =[];
	% recherche uid
   [noms,other]=strtok(nomc,'?');
	if ~isempty(other)
       % recherche de la sequence     
       [uidstr,other]=strtok(other,'+');
		 if ~isempty(other)
       		[void,seqstr]= strtok(other,'=');
 		      seqstr  = strrep(seqstr,'=','');
				seq          = str2num(seqstr);
		 else
		    seq = [];
		 end 
		 % extraction de uid
       [void,uid]=strtok(uidstr,'=');
 		 uid  = strrep(uid,'=','');
	else
       % recherche de la sequence     
       [noms,other]=strtok(nomc,'+');
		 uid = '';
		 if ~isempty(other)
       		[void,seqstr]= strtok(other,'=');
 		      seqstr  = strrep(seqstr,'=','');
				seq          = str2num(seqstr);
		 else
		    seq = [];
		 end 
   end
%  	fprintf('name = %s\n',noms)
%  	if ~isempty(uid)
%  		fprintf('uid = %s\n',uid)
%     end
%  	if ~isempty(seq)
%  		fprintf('seq = %d\n',seq)
%  	end
	% acces a la donnee
    
   [d,t,x,h,err,flag_lec] = zmdsplusjet(numchoc,noms,seq,uid,connect);
%  if flag_lec & (isempty(d) | isempty(t))
%     disp('inconstancy error -> retry')
%     cr = mdsdisconnect;
%     [d,t,x,h,err] = zmdsplusjet(numchoc,noms,seq,uid,0);
%     connect = 1;
%  elseif length(liste) > 1
	if length(liste) > 1
	   connect = 3;
   else
	   connect = 0; 
    end
   
   if ischar(d)
       dc.data = [];
       dc.t    = [];
       dc.x    = x;
       dc.h    = h;
       dc.err  = d;       
   else
       dc.data = double(d);
       dc.t    = double(t);
       dc.x    = x;
       dc.h    = h;
       dc.err  = err;
   end
   % ecriture dans la structure
   nomout  = strrep(noms,'/','.');
   %data=zsetfield(data,nomout,dc);    Instruction plus recente dans la v2.2 de Deneb, mais ne marche pas sur JET ? Je remets la ligne en dessous
   eval(sprintf('data.%s = dc;',nomout));
   % fin de la boucle
end
