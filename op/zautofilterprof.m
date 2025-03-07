function y_out = zautofilterprof(t,x,y,mode)

[u,k,v] = svd(y,0) ;
k       = diag(k);
fprintf('k1 = %g, k2 = %g, k3 = %g, sum(k(4:end)) = %g\n',k(1),k(2),k(3),sum(k(4:end)));
e       = k.^2;
e       = sqrt(sum(e(1:3))/sum(e));
fprintf('Energy kept = %4g %%\n',e*100);

% mise en forme
% valeur entre -1 et 1
for ik = 1:3
    up = max(u(:,ik));
    um = min(u(:,ik));
    if (up > 0) & (um >= 0)
        mu  = 1 ./ up;
    elseif (up <= 0) & (um < 0)
        mu  = 1 ./ um;
    elseif abs(um) > up
        mu  = 1 ./ um;
    else
        mu  = 1 ./ up;
    end
    vp = max(v(:,ik));
    vm = min(v(:,ik));
    if (vp > 0) & (vm >= 0)
        mv  = 1 ./ vp;
    elseif (vp <= 0) & (vm < 0)
        mv  = 1 ./ vm;
    elseif abs(vm) > vp
        mv  = 1 ./ vm;
    else
        mv  = 1 ./ vp;
    end
    
    k(ik)   = k(ik) ./ mv ./ mu;
    u(:,ik) = u(:,ik) .* mu;
    v(:,ik) = v(:,ik) .* mv;
    
end

switch mode
    case 'bezier'
        % bezier fit
        pp1 = fitbezier(t,u(:,1));
        pp2 = fitbezier(t,u(:,2));
        pp3 = fitbezier(t,u(:,3));
        
        % reconstruction
        un = u;
        un(:,1) = bezierval(pp1,t);
        un(:,2) = bezierval(pp2,t);
        un(:,3) = bezierval(pp3,t);
        
        % bezier fit
        pp1 = fitbezier(cat(2,-x(end:-1:1),x),cat(1,v(end:-1:1,1),v(:,1)));
        pp2 = fitbezier(cat(2,-x(end:-1:1),x),cat(1,v(end:-1:1,2),v(:,2)));
        pp3 = fitbezier(cat(2,-x(end:-1:1),x),cat(1,v(end:-1:1,3),v(:,3)));
        
        vn =v;
        %  	vn(:,1) = bezierval(pp1,x);
        %  	vn(:,2) = bezierval(pp2,x);
        %  	vn(:,3) = bezierval(pp3,x);
    otherwise
        % reconstruction
        un = u;
        un(:,1) = sgolayfilt(u(:,1),1,11);
        un(:,2) = sgolayfilt(u(:,2),1,11);
        un(:,3) = sgolayfilt(u(:,3),1,11);
        vn =v;
end

figure
subplot(2,3,1)
plot(t,un(:,1),'b',t,u(:,1),'or');
title(sprintf('k = %g',k(1)));
subplot(2,3,2)
plot(t,un(:,2),'b',t,u(:,2),'or');
title(sprintf('k = %g',k(2)));
subplot(2,3,3)
plot(t,un(:,3),'b',t,u(:,3),'or');
title(sprintf('k = %g',k(3)));
subplot(2,3,4)
plot(x,vn(:,1),'b',x,v(:,1),'or');
subplot(2,3,5)
plot(x,vn(:,2),'b',x,v(:,2),'or');
subplot(2,3,6)
plot(x,vn(:,3),'b',x,v(:,3),'or');



% recontruction
y_out = un(:,1) * k(1,:) * vn(:,1)' + ...
    un(:,2) * k(2,:) * vn(:,2)' + ...
    un(:,3) * k(3,:) * vn(:,3)' ;


