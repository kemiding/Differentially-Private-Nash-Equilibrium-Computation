function [y_new] = projection2(y,radius,center,N)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

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
