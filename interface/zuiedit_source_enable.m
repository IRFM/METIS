%  ZUIEDIT_SOURCE_ENABLE  courte description  
%------------------------------------------------------------------------------- 
% fichier :  zuiedit_source_enable.m  ->  zuiedit_source_enable 
% 
% 
% fonction Matlab 5 : 
% 
% description de la fonction sur quelques lignes  
%  
% syntaxe :  
%   function   zuiedit_source_enable(h,source,value,tag) 
%  
% entrees :  
%  h      = 
%  source = 
%  value  = 
%  tag    = 
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
function zuiedit_source_enable(h,source,value,tag)

     fctd0=strcat('zuidisable(h.','fct_',source,')');
     fctd1=strcat('zuidisable(h.','edit_',source,')');
     fctd2=strcat('zuidisable(h.','chng_',source,')');
     fctd3=strcat('zuidisable(h.','para_',source,')');
     fctd4=strcat('zuidisable(h.','edit_',source,')');
     fctd5=strcat('zuidisable(h.','refval_',source,')');
     fctd6=strcat('zuidisable(h.','presprof_',source,')');

     fcte0=strcat('zuienable(h.','fct_',source,')');
     fcte1=strcat('zuienable(h.','edit_',source,')');
     fcte2=strcat('zuienable(h.','chng_',source,')');
     fcte3=strcat('zuienable(h.','para_',source,')');
     fcte4=strcat('zuienable(h.','edit_',source,')');
     fcte5=strcat('zuienable(h.','refval_',source,')');
     fcte6=strcat('zuienable(h.','presprof_',source,')');

     if tag == 1
     switch lower(value)
          case {0}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fctd5);eval(fctd6);
          case {1}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fctd5);eval(fcte6);
          case {2}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fctd4);eval(fcte5);eval(fctd6);
          case {3}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fcte4);eval(fcte5);eval(fcte6);
     end
     elseif tag == 2
     switch lower(value)
          case {0}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fctd6);
          case {1}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fcte6);
          case {2}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fctd4);eval(fctd6);
          case {3}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fcte4);eval(fcte6);
     end
     elseif tag == 3
     switch lower(value)
          case {0}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fctd5);
          case {1}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);eval(fctd5);
          case {2}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fctd4);eval(fcte5);
          case {3}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fcte4);eval(fcte5);
     end
     elseif tag == 4
     switch lower(value)
          case {0}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);
          case {1}
              eval(fctd0);eval(fctd1);eval(fctd2);eval(fctd3);eval(fctd4);
          case {2}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fctd4);
          case {3}
              eval(fcte0);eval(fcte1);eval(fcte2);eval(fcte3);eval(fcte4);
     end
     end



   
