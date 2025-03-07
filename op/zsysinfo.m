% ZSYSINFO retourne des information sur la memoire et le temps cpu utilise
%---------------------------------------------------------------------------
% fichier zsysinfo.m ->  zsysinfo,zcvps
%
%
% fonction Matlab 5 :
%
% Cette retourne le temps cpu utilise, la memoire utilisee et le
% temps ecoule depuis le lancement
% 
%
% syntaxe  :
%  
%       [flops_matlab,cpu,mem,datation] = zsysinfo;
%       
% entrees : aucune
%
% sorties :
% 
%       flops_matlab  = nombre d'operation flottantes dans matlab
%       cpu           = temps cpu totale (s) + i * temps cpu utilisateur (s)
%                       ( la difference est le temps systeme)
%       mem           = memoire reel utilisee (Mo) + i * memoire virtuelle (Mo)
%       datation      = temps ecoule depuis le lancement (s)
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.8, du 20/02/2002.
% 
% 
% liste des modifications : 
%
%   * 11/07/2001 -> correction bug en cas de plantage de l'instruction unix ps
%   * 20/02/2002 -> securite si probleme de format
%   * 28/08/2002 -> patch pour hercule
%   * 28/05/2005 -> mise a niveau pour saturne
%
%--------------------------------------------------------------
%
function [flops_matlab,cpu,mem,datation] = zsysinfo

%flops_matlab      =  nombre de flops (Matlab)
%cpu         =  temps de calcul du cpu (s)
%mem    =  memoire utilisee par le process (o)
%datation   =  temps ecoule depuis le demarrage du process
% premiere estimation grossiere
flops_matlab=NaN;
cpu=cputime;
mem=NaN;
datation=NaN;

% cas de la session trait
% contournement du probleme systeme d'hercule
% les donnees ne sont pas utilisee dans les traitements
if strcmp(getenv('USER'),'trait') | strcmp(getenv('USER'),'devtrait')
	return
end

if isappdata(0,'pid')
	pid = getappdata(0,'pid');
else
	pid = [];
end
if isempty(pid)
	[pid,pid_sess,user]=getidprocess;
	setappdata(0,'pid',pid);
end

try
    [s,t]=unix(sprintf('ps -ocputime,etime,usertime,rssize,vsize -p%d',pid));
    if s~=0
        disp('probleme avec la commande ps :')
	     disp(t)
	return
    end
    t = tseparec(t);
    t = t(end,:);

    [cputime_str,tt] = strtok(t);
    [date_str,tt] = strtok(tt);
    [user_str,tt] = strtok(tt);
    [mem_res_str,tt] = strtok(tt);
    [mem_vir_str,tt] = strtok(tt);

    % converion
    cpu      = zcvps(cputime_str) + i .* zcvps(user_str);
    mem      = zcvps(mem_res_str) + i .* zcvps(mem_vir_str);
    datation = zcvps(date_str);

catch
    disp('Probleme dans zsysinfo :')
    disp(lasterr)
    return
end
    

% fonction de conversion
function n = zcvps(s)

n=NaN;

if s(end) == 'H'
   % temps en heures + minutes decimales 
   ind = min(find(s ==':'));
   s1 =  s(1:(ind-1));
   s3 =  s((ind+1):(end-1));
   n = str2num(s1) .* 3600 + str2num(s3) .* 60;
elseif s(end) == 'D'     	 
   % temps en jours + heures  decimales 
   ind = min(find(s ==':'));
   s1 =  s(1:(ind-1));
   s3 =  s((ind+1):(end-1));
   n = str2num(s1) .* 3600 * 24 + str2num(s3) .* 3600;
elseif s(end) == 'K'     	 
   % memoire en ko
   n = str2num(s(1:(end-1))) ./1024;
elseif s(end) == 'M'     	 
   % memoire en Mo
   n = str2num(s(1:(end-1))) ;
elseif s(end) == 'G'     	 
   % memoire en Go
   n = str2num(s(1:(end-1))) .* 1000;
elseif length(find(s ==':')) ==2
   % temps en  h m s	 
   ind = find(s ==':');
   s1 =  s(1:(ind(1)-1));
   s2 =  s((ind(1)+1):(ind(2)-1));
   s3 =  s((ind(2)+1):end);
   n = str2num(s1) .* 3600 + str2num(s2) .* 60 + str2num(s3);

else
   % temps en minute + s	 
   ind = min(find(s ==':'));
   s1 =  s(1:(ind-1));
   s2 =  s((ind+1):end);
   n = str2num(s1) .* 60 + str2num(s2);
end	 

if isempty(n)
    n = 0;
end    
