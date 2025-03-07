% appel de lhacc_sym pour plot
if post.z0dinput.option.lhmode == 5
    warndlg('There is no LH power is this scenario but a second EC launcher', 'No LH system');
end
% largeur toroidal du coupleur
if ~isfield(post.z0dinput.option,'upshiftmode')
	post.z0dinput.option.upshiftmode = 'linear';
end
if ~isfield(post.z0dinput.option,'fupshift')
	post.z0dinput.option.fupshift = 1;
end


defaultanswer={num2str(post.z0dinput.option.wlh),num2str(post.z0dinput.option.freqlh),num2str(post.z0dinput.option.npar0), ...
               post.z0dinput.option.upshiftmode,num2str(post.z0dinput.option.fupshift)};
prompt={'toroidal width of LH launcher','LH frequency (GHz)','N_/_/_0','upshiftmode','fupshift'};
name='LH deposition';
numlines=1;
answer=inputdlg(prompt,name,numlines,defaultanswer);
pause(1);


if ~isempty(answer)
   
    switch post.z0dinput.option.gaz
        case 1
            agaz    = 1;
            zgaz    = 1;
        case 2
            agaz    = 2;
            zgaz    = 1;
        case 3
            agaz   = 2.5;
            zgaz    = 1;
        case 5
            agaz   = 2.5;
            zgaz   = 1.5;
         case 11
            agaz   = 1 + 11/5;
            zgaz   = (1 + 25 /5) ./ (1 + 5/5);
       otherwise
            agaz    = 4;
            zgaz    = 2;
    end
	wlh   = str2num(answer{1});
	flh   = str2num(answer{2});
	npar0 = str2num(answer{3});
        upshiftmode = answer{4};
	fupshift = str2num(answer{5});

	plh = interp1(post.zerod.temps,post.zerod.plh,post.profil0d.temps,'linear'); 
	Bt = interp1(post.zerod.temps,post.z0dinput.geo.b0,post.profil0d.temps,'linear'); 
	a = interp1(post.zerod.temps,post.z0dinput.geo.a,post.profil0d.temps,'linear'); 
	ip = interp1(post.zerod.temps,post.zerod.ip,post.profil0d.temps,'linear');
   	agaz = agaz .* ones(size(plh));
   	zgaz = zgaz .* ones(size(plh));
  	 % facteur de propagation LHCD
   	qcyl          =  5 .* a .^ 2 .* Bt ./ (ip./1e6) ./ post.profil0d.Raxe(:,end);
 	qbord         = post.profil0d.qjli(:,end);
    switch upshiftmode
	case {'newmodel','newmodel + tail'}
	  upshift       = fupshift .* ones(size(post.profil0d.temps));
	otherwise
	  upshift       = fupshift .* max(eps,qbord ./ qcyl - 1);
	end

	%upshift       = fupshift .* max(eps,qbord ./ qcyl - 1);

%  	[x,fpout,xlh,dlh,xmem,lc,hc,acc,landau] = ...
%  	           z0lhacc_sym(1e9 .* flh,npar0,wlh,agaz,zgaz,post.profil0d.temps, ...
%  		   post.profil0d.xli,post.profil0d.nep,post.profil0d.tep, ...
%  		   post.profil0d.qjli,post.profil0d.Raxe,post.profil0d.rmx, ...
%  		   post.profil0d.spr,post.profil0d.vpr,post.profil0d.fdia, ...
%  		   plh,0.5.* ones(size(plh)),0.25.* ones(size(plh)));
	
	[x,fpout,xlh,dlh,lc,hc,acc,landau,efficiency,effacc] = ...
                  z0lhacc(flh*1e9,npar0,wlh,agaz,zgaz,post.profil0d.temps, ...
                  post.profil0d.xli,post.profil0d.nep,post.profil0d.tep,post.profil0d.qjli, ...
                  post.profil0d.Raxe,post.profil0d.rmx,post.profil0d.spr,post.profil0d.vpr, ...
                  Bt,plh,post.z0dinput.option.xlh .* ones(size(plh)), ...
                  post.z0dinput.option.dlh .* ones(size(plh)),1,upshift,2,upshiftmode);

end	
	
	
	
	
	
	
	