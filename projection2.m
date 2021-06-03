function [y_new] = projection2(y,radius,center,N)
% this function is the projection operator used in Algorithm 2, which project any input y into a spherical set with parameters center and radius
% N is the node number (size of y)
% y_new is the projected point 

dif_y = y - ones(1,N).* center;
y_new = zeros(N,N);

ind = find(vecnorm(dif_y) <= radius);

%if vecnorm(dif_y)<=radius
    y_new(:,ind) = y(:,ind);
%else 
    ind2 = find(vecnorm(dif_y) > radius);
    dif_y_norm = dif_y(:,ind2) * radius./vecnorm(dif_y(:,ind2));
    y_new(:,ind2) = center + dif_y_norm;
%end
end
