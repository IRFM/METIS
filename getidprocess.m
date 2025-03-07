% GETPIDPROCESS recupere le pid du process matlab
%-------------------------------------------------------
%
% syntaxe :
%	[pid,pid_sess,user]=getidprocess
%
%	pid = pid du process matlab
%  pid_sess = pid du process racine de la session
%  user = nom du user
%
% fonction ecrite par J-F Artaud, poste 46-78
% version 1, derniere mise �jour le 08/08/96
%-------------------------------------------------------
%
function [pid,pid_sess,user]=getidprocess


    if ~strcmp(computer,'GLNX86') & ~strcmp(computer,'GLNXA64')

        %
	    % recupere le pid
	    %
	    [s,t]=unix('ps $$ -j');
	    if s ~= 0
		    [s,t]=unix('ps -j $$');
	    end
	    if s~=0
		    error('Probleme d''execution d''une commande unix');
	    end
	    ind=find(t==10);
	    if length(ind)<2,
		    ind(2)=length(t);
	    end
	    t=t(ind(1):ind(2));
	    ind=find(([t,' ']==' ')&([' ',t]~=' '));
	    user=t(1:(ind(1)-1));
	    pid=str2num(t(ind(2):ind(3)));
	    val=t(ind(4):ind(5));
	    if all(abs(val) > 47 & abs(val) < 58)
	      pid_sess=str2num(val);
	    else
              pid_sess=-1;
	    end


    else
            [s,user]=unix('whoami');
	    if s~=0
		    error('Probleme d''execution d''une commande unix');
        else
            user=tseparec(user);
            user=user(end,:);
	    end
            user(user <= ' ') =[];
        
            [s,t]=unix('ps $$ -j');
	    if s ~= 0
		    [s,t]=unix('ps -j $$');
	    end
	    if s~=0
		    error('Probleme d''execution d''une commande unix');
        end
	t(t==13) = [];  % securite pour certaines configurations linux
	t = tseparec(t);
        t = t(end,:);
	
	    %[r,t]=strtok(t,sprintf('\n'));
           [pid,t]=strtok(t);
 	   pid    = str2num(pid);
           [pid,t]=strtok(t);
 	   pid    = str2num(pid);
          [pid_sess,t]=strtok(t);
	   pid_sess    = str2num(pid_sess);
    end
