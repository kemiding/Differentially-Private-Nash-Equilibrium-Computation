function [mean_z,mean_z2,L_a,L_b,L_g,alpha] = Algorithm2(T,s,Inf_G,Com_G,a_star,H,test_no,sigma,N,b,radius,center)
%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明

[rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H,Com_G,N, s);

mean_x = zeros(N,N,T);
mean_y = zeros(N,N,T);
mean_z = zeros(N,N,T);
mean_z2 = zeros(1,N,T);

m = 1;

L_a = radius + vecnorm(center);
L_b = max(b);
L_g = max(max(Inf_G.W));

while m<=test_no
    clear x y z py z2 gamma 
    x = zeros(N,N,T);
    y = zeros(N,N,T);
    py = zeros(N,N,T);
    z = zeros(N,N,T);
    z2 = zeros(1,N,T);
    %mean_gamma = zeros(N,1);
    for k =1:T-1
    gamma = laprnd(N,N,0,sigma);
    y(:,:,k) = x(:,:,k)+gamma; 
    z(:,:,k) = x(:,:,k)- ones(1,N).* a_star;
    z2(:,:,k) = vecnorm(z(:,:,k)).^2;
    py(:,:,k) = projection2(y(:,:,k),radius,center,N);
    for i = 1:N
    diff = ones(1,N).* py(:,i,k)- py(:,:,k); 
    x(:,i,k+1) = py(:,i,k) - sum(Com_G.W(i,:).*diff,2)- s* H(:,i)*(H(:,i)'*py(:,i,k)-(b(i))); 
    clear diff
    end
    end
    y(:,:,T) = x(:,:,T); 
    z(:,:,T) = x(:,:,T)- ones(1,N).* a_star;
    z2(:,:,T) = vecnorm(z(:,:,T)).^2;
    mean_x = ((m-1)*mean_x+x)/m;
    mean_y = ((m-1)*mean_y+y)/m;
    mean_z = ((m-1)*mean_z+z)/m;
    mean_z2 = ((m-1)*mean_z2+z2)/m;
    % mean_gamma = ((m-1)*mean_gamma+gamma)/m;  %%% test the zero mean of noise
    m=m+1
    clear k i
end














end