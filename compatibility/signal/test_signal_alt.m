% test alternative functions
for k = 1:1000
	l = fix(rand(1)* 1001);
	o = 1 + fix(rand(1) * 5);
	w = o +1 + fix(rand(1) * 25);
	ws = w;
	if fix(ws ./ 2) == (ws ./ 2)
		ws = ws + 1;
	end
	x = linspace(0,abs(randn(1) * pi),max(l,ws));
	y = cos(x) + randn(size(x)) .* sin(x) + randn(size(x));
	y1 = medfilt1(y,w);
	y1_alt = medfilt1_alt(y,w);
	d1(k) = max(abs(y1 - y1_alt));
	if d1(k) > sqrt(eps)
		keyboard
	end
	y1 = sgolayfilt(y,o,ws);
	y1_alt = sgolayfilt_alt(y,o,ws);
	sgo(k) = max(abs(y1 - y1_alt))	
	if sgo(k) > sqrt(eps)
		keyboard
	end
end
max(d1)
max(sgo)
