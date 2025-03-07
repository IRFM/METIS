% script de securite pour la configuration des mode de calcul des coef de transport
% warning pour le mode de calcul
if any(data.mode.pe == 2)
	if any(data.mode.ee ~= 2)
      		fprintf('Warning : data.mode.ee is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for electron heat transport \n');
		fprintf('=> values corrected !\n')
		data.mode.ee(:) = 2;
	end	
	if any(data.mode.ve ~= 2)
      		fprintf('Warning : data.mode.ve is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for electron heat transport \n');
		fprintf('=> values corrected !\n')
		data.mode.ve(:) = 2;
	end
	if  param.cons.neomulti.ee ~= 1
      		fprintf('Warning : param.cons.neomulti.ee is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for electron heat transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.ee = 1;
	end
	if  param.cons.neomulti.ve ~= 1
      		fprintf('Warning : param.cons.neomulti.ve is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for electron heat transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.ve = 1;
	end		
end
if any(data.mode.pion == 2)
	if any(data.mode.ii ~= 2)
      		fprintf('Warning : data.mode.ee is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for ion heat transport \n');
		fprintf('=> values corrected !\n')
		data.mode.ii(:) = 2;
	end	
	if any(data.mode.vi ~= 2)
      		fprintf('Warning : data.mode.ve is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for ion heat transport \n');
		fprintf('=> values corrected !\n')
		data.mode.vi(:) = 2;
	end
	if  param.cons.neomulti.ii ~= 1
      		fprintf('Warning : param.cons.neomulti.ii is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for ion heat transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.ii = 1;
	end
	if  param.cons.neomulti.vi ~= 1
      		fprintf('Warning : param.cons.neomulti.vi is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for ion heat transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.vi = 1;
	end		
end
if any(data.mode.nel == 2)
	if any(data.mode.nn ~= 2)
      		fprintf('Warning : data.mode.nn is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for density transport \n');
		fprintf('=> values corrected !\n')
		data.mode.nn(:) = 2;
	end	
	if any(data.mode.vn ~= 2)
      		fprintf('Warning : data.mode.vn is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for density transport \n');
		fprintf('=> values corrected !\n')
		data.mode.vn(:) = 2;
	end
	if  param.cons.neomulti.nn ~= 1
      		fprintf('Warning : param.cons.neomulti.nn is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for density transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.nn = 1;
	end
	if  param.cons.neomulti.vn ~= 1
      		fprintf('Warning : param.cons.neomulti.vn is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for density transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.vn = 1;
	end		
end
if any(data.mode.rot == 2)
	if any(data.mode.rotc ~= 2)
      		fprintf('Warning : data.mode.rotc is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for rotation transport \n');
		fprintf('=> values corrected !\n')
		data.mode.rotc(:) = 2;
	end	
	if any(data.mode.rotv ~= 2)
      		fprintf('Warning : data.mode.rotv is not set to 2\n');
     		fprintf('that can generated, in some case, wrong results for rotation transport \n');
		fprintf('=> values corrected !\n')
		data.mode.rotv(:) = 2;
	end
	if  param.cons.neomulti.rot ~= 1
      		fprintf('Warning : param.cons.neomulti.rot is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for rotation transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.rot = 1;
	end
	if  param.cons.neomulti.rotv ~= 1
      		fprintf('Warning : param.cons.neomulti.rotv is not set to 1\n');
     		fprintf('that can generated, in some case, wrong results for rotation transport \n');
		fprintf('=> value corrected !\n')
		param.cons.neomulti.rotv = 1;
	end		
end
zuisavenonok;