function [rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H,ComG,N, s)

%   check the satisfactory of parameters' conditions 
%   ComG: the structure of the communication graph
%   N: number of node

A = (H')^(-1);
if min(A,[],'all')>=0
    fprintf('Assumption 1 is satisfied as the inverse of I-G is nonnegative.\n');
else
     error('Assumption 1 is not satisfied as the inverse of I-G is not nonnegative!');
end

ComG_L = ComG.L;

lambda = sort(eig(ComG_L));
lambda_n = max(lambda);
lambda_0 = min(lambda);

%[bool,in]=gsp_check_connectivity_undirected( ComG); %% another method to check the connectivity of the communication graph

%if bool ==1
if lambda_0>=0
fprintf('Communication graph is connected.\n');
else 
    in
    error('Communication graph is not connected!')
end

if lambda_n<2
     fprintf('Maximum eigenvalue of communication graph is less than 2.\n');
else
     error('Maximum eigenvalue of communication graph is greater than or equal to 2!');
end

%%%% compute other algorithm related parameter
M = zeros(N,N);
H_norm = zeros(N,1);
for i = 1:N
    M = M + H(:,i)*H(:,i)';
    H_norm(i) = norm(H(:,i));
end
rho_m = min(svd(M/N));  
h_M = max(H_norm);
alpha_1 = (lambda_n+2*s*h_M^2+sqrt(lambda_n^2+4*s^2*h_M^4))/2-1;
alpha_2 = 1- (sqrt((lambda(2)-s*rho_m)^2+4*s^2*h_M^4)+lambda(2)+s*rho_m)/2;
alpha = max(abs([alpha_1,alpha_2]));

s_1 = 2*(2-lambda_n)/(h_M^2*(4-lambda_n));
s_2 =  (sqrt((2-lambda(2))^2*rho_m^2+8*h_M^4*(2-lambda(2)))-(2-lambda(2))*rho_m)/(2*h_M^4);
s_m = min([s_1,s_2]);
if s>0 && s< s_m 
    fprintf('Condition of theorem 2 is satisfied.\n');
else
    fprintf('Condition of theorem 2 fails!');
end

end

