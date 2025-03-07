% ZMAIL envoi un mail a l'aide de la fonction unix Mail
%-------------------------------------------------------------------------------------------
%
% nom du fichier contenant la fonction : zmail.m
%
% descriptif :
%
%   Cette fonction permet d'envoyer un mail a partir de Matlab a une liste de 
% user en donnant le sujet et le texte. Elle interface la fonction unix Mail
% avec Matlab.
% remarque :
%   si un seul user n'est pas connu du mailer, le mail ne sera pas correctement 
% envoye. Le texte du sujet ne doit pas depasser 80 caracteres. Le texte 
% du message doit vaoir moins de 80 caracetres par lignes.
%
% syntaxe :
%
%	*	cr=zmail(adresse,subject,message);
%
% variables d'entree :
%
%	* adresse	-> liste des user a qui le mail doit etre envoye		 
%						(format : user1 user2 user3 ...)
%  * subject   -> sujet du mail (chaine de caractere d'une ligne de moins 
%                 de 80 caracteres qui peut comporter des mots clefs pour
%                 le classement automatique des mails)
%
%	* message  -> texte du message (c'est une chaine de caracteres pouvant etre 
%                composee de plusieurs lignes de moins de 80 caracteres)
%
% variables de sortie :  
%
%  * cr -> compte rendu d'execution (cr=0 => ok)
%
%
% fonction ecrite par J-F Artaud, poste 46-78
% version 2, derniere mise a jour le 04/10/2001
% 
%-------------------------------------------------------------------------------------------
%
function cr = zmail(adresse,subject,message)

		%
		% testd es entrees
		%
		cr=0;
		if nargin <3
			disp('syntaxe : cr=zmail(adresse,subject,message)')
			cr = -1;
			return
		elseif isempty(adresse)
			disp('il faut donner l''adresse');
			cr = -2
			return
		elseif isempty(message)
			cr =-3;
			disp('et le message ?');
			return
		end
		%
		if isempty(subject)
			subject='pas de sujet';
		end
		subject=subject(:)';
		message=message(:)';
		adresse=adresse(:)';
		%
		% elaboration de la commande
		%
		commande=['Mail -s "',subject,'" ',adresse,' <<@@',10, ...
				message, 10,'.',10,'@@',10 ];
		%
		% envoi
		%	
		[cr,t]=zunix(commande);
		if cr~=0
			disp('erreur lors de l''envoi du mail ...')
			disp(t)
			disp(commande)
		else
				disp(['envoi d''un mail : ',subject])
		end
		%
		% fin de la fonction
		%
