% must be change for steady state
%f_ni =  1; % for steady state
%f_ni =  0.7; % must be set to value lowe rthan one if shot is not steady state
mode_ni        = input('Type of optimisation: 0 -> Standard; 1 -> Hybrid; 2 -> Steady state ? ');
%radial_ext     = input('Tangencial radius flexibility in unit of minor radius [0.1 0.8] ? ');
%vertical_shift = input('Normalized vertical shift flexibility [0.1 0.8] ? ');
% number of loop
kmax = 31; % from 1 to 31
% memorize original NBI power
pnbi_mem = z0dinput.cons.pnbi;
%fpnbi_fact = 2;
fpnbi_fact = sqrt(2);
% flag fast /full
ffflag = 0;
% tableau des resultats
res_tab = [];
% loop for optimisation
for k=1:kmax
    fprintf('=====> start of optimization epoch %d (over maximum of 31) <=====\n',k);
    % attenuation de la variation de  lapuissance nbi
    if k > 21
	fpnbi_fact = fpnbi_fact ./ sqrt(2);
    end
    % shine_through_limit
    shine_through_limit = 0.05 .* max(real(post.zerod.pnbi)) +  0.05 .* max(imag(post.zerod.pnbi));
    fprintf('Shine through limit set at %g MW\n',shine_through_limit./1e6);
    % NBI optimisation
    disp('Optimisation of NBI geometry:')
    % X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
    options_optim = optimset('display','iter','Algorithm','active-set','MaxFunEvals',1e4);
    switch mode_ni
    case 0
      res = cat(2,post.z0dinput.option.rtang,post.z0dinput.option.rtang2,post.z0dinput.option.zext,post.z0dinput.option.zext2);
      lb  = [max(z0dinput.geo.R) - 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) -  0.8 .* max(z0dinput.geo.a),0,0];
      ub  = [max(z0dinput.geo.R) + 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) +  0.8 .* max(z0dinput.geo.a),0.7,0.7];
    otherwise
      res = cat(2,post.z0dinput.option.rtang,post.z0dinput.option.rtang2,post.z0dinput.option.zext,post.z0dinput.option.zext2,1,1);
      lb  = [max(z0dinput.geo.R) - 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) -  0.8 .* max(z0dinput.geo.a),0,0,1./fpnbi_fact,1./fpnbi_fact];
      ub  = [max(z0dinput.geo.R) + 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) +  0.8 .* max(z0dinput.geo.a),0.7,0.7,fpnbi_fact,fpnbi_fact];    
    end
    if k == 1  
	res_tab = res(:)';
    end
    [res,objective_end,exit_flag] = fmincon(@(x) z0optimnbigeo(x,post,shine_through_limit,mode_ni),res,[],[],[],[],lb,ub,'',options_optim);
    [objective,rtang,rtang2,zext,zext2,ialign,fini,soq,shine,qnot,qmin] = z0optimnbigeo(res,post,shine_through_limit,mode_ni);
    if k > 11
	z0dinput.option.rtang  = 0.3 .* res(1) + 0.7 .* z0dinput.option.rtang;
	z0dinput.option.rtang2 = 0.3 .* res(2) + 0.7 .* z0dinput.option.rtang2;
	z0dinput.option.zext   = 0.3 .* res(3) + 0.7 .* z0dinput.option.zext;
	z0dinput.option.zext2  = 0.3 .* res(4) + 0.7 .* z0dinput.option.zext2;
    else
	z0dinput.option.rtang  = res(1);
	z0dinput.option.rtang2 = res(2);
	z0dinput.option.zext   = res(3);
	z0dinput.option.zext2  = res(4);   
    end
    fprintf('objective (to be minimized): %g\n',objective);
    fprintf('non inductive fraction: %g\n',fini);
    fprintf('current alignement: %g\n',ialign);
    fprintf('qnot = %g & qmin = %g\n',qnot,qmin);
    fprintf('shine through (and first orbit losses) & maximum allowed shine through (MW): %g <? %g\n',shine ./ 1e6 ,shine_through_limit ./ 1e6 );
    fprintf('confinement criterium (<|s|/q>): %g\n',soq);
    fprintf('NBI1: Rtang = %g m & vertical shift = %g (normalized)\n',res(1),res(3));
    fprintf('NBI2: Rtang = %g m & vertical shift = %g (normalized)\n',res(2),res(4));
    switch mode_ni
    case 0
	% rien
    otherwise
	if k > 11
	      z0dinput.cons.pnbi = (0.3 .* res(5) + 0.7) .* real(z0dinput.cons.pnbi) + sqrt(-1) .* (0.3 .* res(6) + 0.7) .* imag(z0dinput.cons.pnbi);
	else
	      z0dinput.cons.pnbi = res(5) .* real(z0dinput.cons.pnbi) + sqrt(-1) .* res(6) .* imag(z0dinput.cons.pnbi);
	end
    end
    fprintf('NBI power multiplied by: NBI_1 = %g (%g MW) & NBI_2 = %g (%g MW)\n', ...
            max(real(z0dinput.cons.pnbi)) ./ max(real(pnbi_mem)),max(real(z0dinput.cons.pnbi)) ./ 1e6 ,  ...
            max(imag(z0dinput.cons.pnbi)) ./ max(imag(pnbi_mem)),max(imag(z0dinput.cons.pnbi)) ./ 1e6);	  
            
    % graphe of results         
    res_tab(end+1,:) = res(:)';
    figure(17);clf;plot(1:size(res_tab,1),res_tab);drawnow;
    % break condition
    if k == 1
	res_mem = res;
    else
	dres = sqrt(sum((res - res_mem) .^ 2))./ length(res); 
	res_mem = res;
	if dres < sqrt(eps) 
	    if ffflag == 0
	      ffflag = 1;
	    else
		break;
	    end
	end
    end
    % call of METIS
    switch ffflag
    case 0
	  [zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
    otherwise
	  [zs,infovoid,profli] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);    
    end
    post.z0dinput = z0dinput;
    post.zerod    = zs;
    post.profil0d =profli;
    z0plotsc;drawnow;
end
