%  NOM DE LA FONCTION  courte description  
%------------------------------------------------------------------------------- 
% fichier : nom_du_fichier -> nom_fonction_principale, nom_fonction_locale 1,... 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function .... 
%  
% entrees :  
%  
%  
% sorties :  
%  
%  
% fonction ecrite par xxxxxxx , poste XX-XX  
% version  4.1  du  08/04/2008  
%  
%*@auto@   Help genere automatiquement  
%  
% liste des modifications :  
%   * XX/XX/20XX -> commentaires  
%  
%-------------------------------------------------------------------------------  
%  
chemin = getappdata(0,'root');

if exist(fullfile(chemin,'WAVE'))
    if exist(fullfile(chemin,'EQUIL'))
        [reseq,fileq] = unix(['ls ',chemin,'/EQUIL/*.mat']);
    end
    
    chemin = fullfile(chemin,'WAVE','C3PO');
    
    [res,file] = unix(['ls ',chemin,'/WAVE*EAST*.mat']);
    posmat     = strfind(file,'.mat');
    numfile = length(posmat);
    if numfile > 1
       
        disp(file)
        choixfile = 1;
        choixfile = inputd('more than one wave file,choose which file you want',choixfile);
        if choixfile == 1
          ndeb      = 1;
          nfin      = posmat(1)+3;
        else
          ndeb      = posmat(choixfile-1)+4;
          nfin      = posmat(choixfile)+3;
        end
        file=file(ndeb:nfin)
        
    end
    if ~exist(file)
        indd = strfind(file,'/');
        if indd(1) ~= 1
           file(1:indd(1)-1)=[]; 
        end
        indf = strfind(file,'.mat');
	file(indf+4:end)=[];
    end
    val=load(file);
    %
    % indice
    %
    indk = strfind(file,'_k');
    indk = indk+2;
    indfk = strfind(file(indk:end),'_');
    indfk = indfk(1)-1+indk;
    kluke = sscanf(file(indk:indfk),'%f');
    datak = zget1t(data,kluke);
    geo = datak.geo;  
    ray = val.wave.rays{1};
    
    h = findobj(0,'type','figure','tag','rayon_luke');
    if isempty(h)
       h=figure('tag','rayon_luke');
    else
       figure(h);
    end   
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

    
    subplot(2,2,1)
    
    plot(ray.sx,ray.sy,geo.R-geo.r0,geo.Z)
    axis('equal')
    xlabel('x')
    ylabel('y')
    title(['ray tracing t=',num2str(data.gene.temps(kluke),3),' s'])
    
    subplot(2,2,2)
    
    plot(ray.ss,abs(ray.sNpar),ray.ss,6.5 ./ sqrt(ray.sTe))
    title('N//')
    xlabel('s')
    axis([-inf inf 0 max(abs(ray.sNpar))*2])
    
    subplot(2,2,3)
    
    plot(ray.ss,ray.srho)
    
    title('minor radius')
    xlabel('s')
    
    subplot(2,2,4)
    
    plot(ray.ss,ray.sTe)
    title('Te')
    xlabel('s')
    
end


