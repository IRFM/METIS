%  ZUIEDIT_SOURCE_ACTIVE  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zuiedit_source_active.m  ->  zuiedit_source_active 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   [col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active(source,type,popupmode,valeurmode,nb1,nb2,fonction,tag) 
%  
% entrees :  
%  source     = 
%  type       = 
%  popupmode  = 
%  valeurmode = 
%  nb1        = 
%  nb2        = 
%  fonction   = 
%  tag        = 
%  
% sorties :  
%  col1 = 
%  col2 = 
%  col3 = 
%  col4 = 
%  col5 = 
%  col6 = 
%  col7 = 
%  col8 = 
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
function [col1,col2,col3,col4,col5,col6,col7,col8]=zuiedit_source_active(source,type,popupmode,valeurmode,nb1,nb2,fonction,tag)

     tag1=source;
     tag2=type;
     tag3=strcat('info.param.fonction.',tag1); 
     tag4=strcat('data.mode.',tag1);
     col1 = {tag1,'text',tag2,nb1,tag3,''} ;
     fct = getfield(fonction,source) ; if isempty(fct) fct=' '; end
     val = evalin('base',tag4) ;
     ind = fctval(val) + 1 ;
     tag4=strcat('mode_',tag1);
     tag5=strcat('fct_',tag1);
     tag6=strcat('chng_',tag1);
     tag7=strcat('para_',tag1);
     tag8=strcat('edit_',tag1);
     tag9=strcat('refval_',tag1);
     tag10=strcat('presprof_',tag1);
     col2 = {tag4, 'popup',popupmode            ,ind,''              ,valeurmode,'','Enable','on'};
     %all the column
     if tag == 1 
     switch lower(ind)
          case {1}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','off'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','off'};
          case {2}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','off'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','on'};
          case {3}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','on'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','off'};
          case {4}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','on'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','on'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','on'};
     end
     elseif tag == 2 %all column except col7 (reference value)
     col7 = {'sources','text','',[],'',''} ;
     switch lower(ind)
          case {1}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','off'};
          case {2}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','on'};
          case {3}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','off'};
          case {4}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','on'};
              col8 = {tag10,'radio','prescribed profiles',0  ,''              ,[]        ,'','Enable','on'};
     end
     elseif tag == 3 %all column except col8 (prescribed profiles)
     col8 = {'sources','text','',[],'',''} ;
     switch lower(ind)
          case {1}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','off'};
          case {2}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','off'};
          case {3}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','on'};
          case {4}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','on'};
              col7 = {tag9, 'radio','reference values'   ,0  ,''              ,[]        ,'','Enable','on'};
     end
     elseif tag == 4
     col7 = {'sources','text','',[],'',''} ;
     col8 = {'sources','text','',[],'',''} ;
     switch lower(ind)
          case {1}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
          case {2}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','off'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','off'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','off'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
          case {3}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','off'};
          case {4}
              col3 = {tag5, 'text' ,fct                  ,nb2,''              ,[]        ,'','Enable','on'};
              col4 = {tag6, 'radio','change module'      ,0  ,''              ,[]        ,'','Enable','on'};
              col5 = {tag7, 'radio','parameters'         ,0  ,''              ,[]        ,'','Enable','on'};
              col6 = {tag8, 'radio','edition mode'       ,0  ,''              ,[]        ,'','Enable','on'};
     end
     end

