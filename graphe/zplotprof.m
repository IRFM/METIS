% PLOTPROF affichage de profil, avec fenetre de defilement du temps
%-----------------------------------------------------------------------
%
% fonctions plotprof.m
%
% 	Cette fonction affiche un profil ou plusieurs profils a un temps donne. 
%	Le temps peut etre balayer  l'aide d'un slider ou saisie dans la fenetre. 
%	
% syntaxe :
%
%	plotprof(action),
%
%			-> gestion des evenement (ne pas utiliser directement)
%
% ou 
%
%	 [ha,hc,hl]=plotprof(handle_axes,temps,espace,valeur,prop1,val1, ...)
%
%
% entree
%			
%	handle 		= 	handle des axes (optionnel). Permet de superposer 
%					plusieur profils dans un meme axe (cf exemple) .
%
%	temps 		= 	vecteur des temps (size(temps,2)=1 !),
%
%	espace 		= 	vecteur ou matrice des coordonnees d'espace du profil
%
%	valeur		= 	matrice des valeurs valeurs(temps,espace)
%                               plot un ligne si reel ou un patch si complexe
%                               a la maniere d'une barre d'erreur de largeur
%                               [real(valeur) imag(valeur)]
%	
%	propX, valX	=	couple (nom_propriete,valeur) appliquer a la ligne qui sert a tracer le profil
%					6 couples peuvent etre utilisees.
%	
% sortie :
%			
%	ha			= handles des axes du profil
%	hc			= handle de la fenetre de controle
%	hl			= handle de la ligne dans le cas d'une nouvelle ligne
%
% exemple
%
%	[ne,t,y]=tsbase(16421,'GTETHOM');
%	ha=plotprof(gca,t,y,ne,'color',[1 0 0],'linestyle','o');
%
% fonction ecrite par J-F Artaud, poste 62-15
% version 3, derniere mise a jour le 18/07/2005
%	* 18/07/2005  -> adaptation matlab 6/7
%
%-----------------------------------------------------------------------
%
function [ha,hl]=plotprof(ha,temps,espace,valeur,n1,p1,n2,p2,n3,p3,n4,p4,n5,p5,n6,p6)

	action='new';
	nbpaire	= 0;
	hf =[];
	%
	% gestion des entrees
	%
	if nargin ==1,
		action =ha;
	elseif nargin == 2
		action =ha;
		hf = temps;		
	elseif nargin < 4
		error('nombre de parametres incorrect');
	elseif isempty(ha)
		ha=gca;
	elseif ~ishandle(ha)
		ha=gca;
	end
	if strcmp(action,'new')
		type_ha=get(ha,'type');
		if strcmp(type_ha,'figure')
			ha=get(ha,'currentaxes');
		elseif ~strcmp(type_ha,'axes')
			ha=get(ha,'parent');
		end
	end
	%
        % figure qui porte les axes
	%
	if isempty(hf)
		hf=gcf;
	end
	%
	% recherche des uicontrols 
	%
	hcurs=findobj_local(hf,'tag','curseur');
	htemps=findobj_local(hf,'tag','temps');
	hmin=findobj_local(hf,'tag','minimum');
	hmax=findobj_local(hf,'tag','maximum');
	%
	hl=[];
	ttmax  = -inf;
	ttmin  = inf;
	%
	% mise en place si necessaire
	%
	if strcmp(action,'new')
		%
		% test du nombre de paires de parametres
		%
		if nargin >4
			nbpaire=(nargin-4)/2;
			if nbpaire ~=fix(nbpaire)
				error('il faut des proprietes de la forme <nom>,<valeur>');
			end
		end
		%
		% test des dimensions
		%
		if isempty(temps)
			temps=(1:size(valeur,1))';
		elseif size(temps,1) ==1
			temps = temps';
		end
		if size(temps,1)~=size(valeur,1)
			error('erreur de dimension 1 temps-valeur');
		elseif (size(temps,2) >1)&(size(temps,2)~=size(valeur,2))
			error('erreur de dimension 2 temps-valeur');
		elseif size(temps,2) == 1
			temps=temps*ones(1,size(valeur,2));
		end
		if isempty(espace),
			espace=1:size(valeur,2);
		elseif size(espace,2)==1
			espace=espace';
		end
		if size(espace,2)~=size(valeur,2)
			error('erreur de dimension 2 espace-valeur');
		elseif (size(espace,1) >1)&(size(espace,1)~=size(valeur,1))
			error('erreur de dimension 1 espace-valeur');
		elseif size(espace,1) == 1
			espace=ones(size(valeur,1),1)*espace;
		end
		%
		% mise en place de la figure portant le curseur
		%
		if isempty(hcurs)
                        set(hf,'resizefcn','zplotprof(''resize'');');
			%
			% creation des uicontrols
		        %
     		 	 uicontrol(hf,'Style','frame','Units','normalized' , ...
    		 		          'Position',[0 0 0.999 0.05], ...
                                    'Backgroundcolor',[0.6 0.6 0.6]); 

   		 	 htemps = uicontrol(hf,'Style','edit','Units','normalized' , ...
    		 	                   'Position',[0.1 0.003 0.2 0.04], ...
                                           'String','0', ...
                                           'HorizontalAlignment','left', ...
                                           'tag','temps', ...
                                           'userdata',0, ...
                                           'Callback','zplotprof(''temps'');'); 

    		 	 uicontrol(hf,'Style','text','Units','normalized' , ...
    		 		          'Position',[0 0.003 0.1 0.04], ...
                           	          'String','t = ','HorizontalAlignment','left'); 

     		 	 uicontrol(hf,'Style','text','Units','normalized' , ...
     		 		  'Position',[0.3 0.003 0.05 0.04], ...
                           	  'String',' s','HorizontalAlignment','left'); 

   		 	 hcurs = uicontrol(hf,'Style','slider','Units','normalized' , ...
    		 	                     'Position',[0.46 0.003 0.439 0.04], ...
                                             'Min',-1e308,'Max',1e308, ...
                                             'tag','curseur', ...
                                             'Value',0,'Callback','zplotprof(''curseur'');');

     	     		 hmin =  uicontrol(hf,'Style','text','Units','normalized' , ...
     	     		                      'Position',[0.36 0.003 0.1 0.04], ...
                                               'String',['0',' s'], ...
                                               'tag','minimum', ...
                                              'HorizontalAlignment','left');

     		 	 hmax = uicontrol(hf,'Style','text','Units','normalized' , ...
     		 	                   'Position',[0.899 0.003 0.1 0.04], ...
                                           'String',['0', ' s'], ...
                                           'tag','maximum', ...
                                           'HorizontalAlignment','right');
		end
		%
		% mise en place du plot des profils
		%
		set(0,'currentfigure',hf);
		set(hf,'currentaxes',ha);
                %
                % selon valeur
                %
                if iscomplex(valeur)
          	  %
		  % creation de lu patch du profil
	 	  %
		  % 
		  % donnees pour le patch
		  %
		  xx=[espace(1,:),espace(1,size(espace,2):-1:1)];
		  yp=real(valeur(1,:));
		  ym=imag(valeur(1,:));
		  yy=[yp,ym(length(ym):-1:1)];
		  %
		  % tracer
		  %
		  hl=patch(xx,yy,'y','erasemode','xor','edgecolor','none');
%                   hlp(1)=line(espace(1,:),real(valeur(1,:)),'color',[1 1 0],'linestyle','-', ...
%                          'tag','size','userdata',hl,'erasemode','normal');
%                   hlp(2)=line(espace(1,:),imag(valeur(1,:)),'color',[1 1 0],'linestyle','-', ...
%                          'tag','size','userdata',hl,'erasemode','normal');
                else
         	  %
		  % creation de la ligne du profil
	 	  %
		  hl=line(espace(1,:),valeur(1,:));
                end
		%
		% propriete graphique de la ligne
		%
		if nbpaire > 0
			set(hl,n1,p1);
		end
		if nbpaire > 1
			set(hl,n2,p2);
		end
		if nbpaire > 2
			set(hl,n3,p3);
		end
		if nbpaire > 3
			set(hl,n4,p4);
		end
		if nbpaire > 4
			set(hl,n5,p5);
		end
		if nbpaire > 5
			set(hl,n6,p6);
		end
		if nbpaire > 6
			set(hl,n7,p7);
		end
%                 %
%                 % bord meme couleur que centre pur le patch
%                 %
%                 if strcmp('patch',get(hl,'type'))
%                      set(hlp,'color',get(hl,'facecolor'));
%                 end
		%
		% sauvegarde des autres profils
		%
		set(hl,'UserData',[temps,espace,valeur],'tag','zplotprof');
                %
                % mise a jour des champs
                % 
		ind=isfinite(temps);
		tmin=min(min(temps(ind)));
		tmax=max(max(temps(ind)));
		if (get(hcurs,'min')>tmin) | (get(hcurs,'min') <= -1e308)
			set(hcurs,'min',tmin);
		end
		if (get(hcurs,'max') < tmax) | (get(hcurs,'max') >= 1e308)
			set(hcurs,'max',tmax);
		end
		if get(hcurs,'max') < get(hcurs,'value')
			set(hcurs,'value',get(hcurs,'max'));
		end
		if get(hcurs,'min') > get(hcurs,'value')
			set(hcurs,'value',get(hcurs,'min'));
		end
		if get(hcurs,'min') > get(hcurs,'max')
			set(hcurs,'max',get(hcurs,'min') + eps);
		end
		set(htemps,'string',num2str(get(hcurs,'value')));
		set(hmin,'string',[num2str(get(hcurs,'min')),' s']);
		set(hmax,'string',[num2str(get(hcurs,'max')),' s']);


      %drawnow
      %zplotprof('resize');

	elseif strcmp(action,'temps')
		%
		% mise a jour des champs  de controle
		%
		tc=str2num(get(htemps,'string'));
		if isempty(tc)
			set(htemps,'string',num2str(get(hcurs,'value')));
			return
		elseif ~isfinite(tc)
			set(htemps,'string',num2str(get(hcurs,'value')));
			return
		end
		if get(hcurs,'min')>tc
			tc=get(hcurs,'min');
			set(htemps,'string',num2str(get(hcurs,'min')));
		end
		if get(hcurs,'max')<tc
			tc=get(hcurs,'max');
			set(htemps,'string',num2str(get(hcurs,'max')));
		end
		set(hcurs,'value',tc);
	elseif strcmp(action,'curseur')
		%
		% mise a jour des champs de la fenetre de control
		%
		set(htemps,'string',num2str(get(hcurs,'value')));
	end
	%
	% mise a jour des profils
	%
	if isempty(hcurs)
		return
	end
	val=get(hcurs,'value');
	if isempty(val)
		return
	elseif ~isfinite(val)
		return
	end
	%
	% recherche des objets de la fenetre de type ligne profil
	%
	hli=findobj_local(hf,'type','line','tag','zplotprof')';
	if ~isempty(hli)
           hl=hli;
	   for k=hli
		%
		% recherche de la nouvelle valeur
		%
		data_line=get(k,'UserData');				
		%
		% separation des donnees
		%
		sil=size(data_line,2)/3;
		if fix(sil)~=sil,
			error('erreur dans les dimension des donnees du profil');
		end
		temps=data_line(:,1:sil);
		ttmin = min(ttmin,min(temps(isfinite(temps))));
		ttmax = max(ttmax,max(temps(isfinite(temps))));
		
                if (( val>=min(min(temps(isfinite(temps)))))&(val<=max(max(temps(isfinite(temps))))))| ...
	           ((size(temps,1)==1)&(abs(mean(temps)-val)<0.2))
		  ind=(1:size(temps,1))'*ones(1,size(temps,2));
		  espace=data_line(:,(sil+1):(2*sil));
		  valeur=data_line(:,(2*sil+1):size(data_line,2));
		  %
		  % recherche du nouveau temps
		  %
                  delta_time = (max(max(temps(:)),max(val(:))) - min(min(temps(:)),min(val(:)))) ./ max(ind(:));
		  diff=abs(temps-val)+sqrt(eps) .* delta_time * ind;
		  indt=(diff==(ones(size(temps,1),1)*min(diff)));
		  set(k,'xdata',espace(indt),'ydata',valeur(indt),'visible','on');
                else
                  set(k,'visible','off');
                end
	    end	
        end	
 	%
	% recherche des objets de la fenetre de type patch profil
	%
	hli=findobj_local(hf,'type','patch','tag','zplotprof')';
	if ~isempty(hli)
	   hl=hli;
	   for k=hli
		%
		% recherche de la nouvelle valeur
		%
		data_line=get(k,'UserData');				
		%
		% separation des donnees
		%
		sil=size(data_line,2)/3;
		if fix(sil)~=sil,
			error('erreur dans les dimension des donnees du profil');
		end
		temps=data_line(:,1:sil);
		ttmin = min(ttmin,min(temps(isfinite(temps))));
		ttmax = max(ttmax,max(temps(isfinite(temps))));
                if (( val>=min(min(temps(isfinite(temps)))))&(val<=max(max(temps(isfinite(temps))))))| ...
	           ((size(temps,1)==1)&(abs(mean(temps)-val)<0.2))
		  ind=(1:size(temps,1))'*ones(1,size(temps,2));
		  espace=data_line(:,(sil+1):(2*sil));
		  valeur=data_line(:,(2*sil+1):size(data_line,2));
		  %
		  % recherche du nouveau temps
		  %
		  diff=abs(temps-val)+1e-100*ind;
		  indt=min(find((diff==(ones(size(temps,1),1)*min(diff)))));
	          % 
		  % donnees pour le patch
	          %
		  xx=[espace(indt,:),espace(indt,size(espace,2):-1:1)];
		  yp=real(valeur(indt,:));
		  ym=imag(valeur(indt,:));
		  yy=[yp,ym(length(ym):-1:1)];
 		  set(k,'xdata',xx,'ydata',yy,'visible','on');
                  hll=findobj_local(hf,'type','line','tag','size','userdata',k);
                  if ~isempty(hll)
                    set(hll(1),'xdata',espace(indt,:),'ydata',real(valeur(indt,:)),'visible','on');
                    set(hll(2),'xdata',espace(indt,:),'ydata',imag(valeur(indt,:)),'visible','on');
                  end
              else
                  set(k,'visible','off');
                  hll=findobj_local(hf,'type','line','tag','size','userdata',k);
                  set(hll,'visible','off');
               end
	   end	
        end	
        %
	% mise a jour de limites
	%
	if isfinite(ttmin)
		set(hcurs,'min',ttmin);
	end
	if isfinite(ttmax)
		set(hcurs,'max',ttmax);
	end
	if get(hcurs,'max')<get(hcurs,'value')
		set(hcurs,'value',get(hcurs,'max'));
	end
	if get(hcurs,'min')>get(hcurs,'value')
		set(hcurs,'value',get(hcurs,'min'));
	end
	set(htemps,'string',num2str(get(hcurs,'value')));
	set(hmin,'string',[num2str(get(hcurs,'min')),' s']);
	set(hmax,'string',[num2str(get(hcurs,'max')),' s']);

        %
        % nomalisation
        %
	if strcmp(action,'resize')
            huic=findobj_local(hf,'type','uicontrol');
            hfr=findobj_local(hf,'type','uicontrol','style','frame');
	    set(huic,'units','pixels');		     
	    for k=1:length(huic)
                  pos=get(huic(k),'position');
                  if (huic(k)==hfr)
                     pos(4)=20;
                  else
                     pos(4)=16;
                  end
                  set(huic(k),'position',pos);
            end
            %drawnow
	    set(huic,'units','normalized');
         end
 
