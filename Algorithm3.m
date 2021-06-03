function [mean_z,mean_z2,upp_bound] = Algorithm3(T,s,Com_G,a_star,H,test_no,sigma,N,b)
%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明

[rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H,Com_G,N, s);

mean_x = zeros(N,N,T);
mean_y = zeros(N,N,T);
mean_z = zeros(N,N,T);
mean_z2 = zeros(1,N,T);

%test_no = 1;
%sigma =1; % noise variance

m = 1;
while m<=test_no
    clear x y z z2 gamma 
    x = zeros(N,N,T);
    y = zeros(N,N,T);
    z = zeros(N,N,T);
    z2 = zeros(1,N,T);
    omega = laprnd(N,1,0,sigma);
    eta = laprnd(N*(N+1)/2,1,0,sigma);
    %x(:,:,1)= 0.5 * ones(N,N);
    %mean_gamma = zeros(N,1);
    H_new = zeros(N,N,N);
    for i=1:N
    H_new(:,:,i)=proj_func(N,H(:,i),eta,eye(N)-H',i);
    rank(H_new(:,:,i))
    end
    
    for k =1:T-1
    y(:,:,k) = x(:,:,k); 
    z(:,:,k) = x(:,:,k)- ones(1,N).* a_star;
    z2(:,:,k) = vecnorm(z(:,:,k)).^2;
    for i = 1:N
    diff = ones(1,N).* y(:,i,k)- y(:,:,k);
    x(:,i,k+1) = y(:,i,k) - sum(Com_G.W(i,:).*diff,2)- s*H_new(:,:,i)*y(:,i,k) + s*H(:,i)*(b(i)+omega(i)); 
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

upp_bound = N* s^2*sigma^2*h_M^2/(1-alpha)^2;

end

