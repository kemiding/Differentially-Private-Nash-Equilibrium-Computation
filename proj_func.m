%%%%%%%%%%%%projection function %%%%%%%%%%%%%%%%%


function H_new = proj_func(N,h_i,eta,G,i)
%%%%h_i is a vector of dimension N
%%%%G is a matrix of N*N
%%%%eta is a noise vector of dimension n_dim = N*(N+1)/2

%%%%bar q_i = 1.5*|h_i|^2, \underline{q}_i = 0.5*|h_i|^2

n_dim = N*(N+1)/2; %%%dimension of variables for optimization


f = -2*(matrix_vec(N,h_i*h_i')+eta);



Q = zeros(N^2,n_dim);

I_NN = eye(n_dim);

for j=1:N
   for k=1:N
       if k>=j
           Q(k+(j-1)*N,:) = I_NN(N*(j-1)+k-j*(j-1)/2,:);
       else
           Q(k+(j-1)*N,:) = I_NN(N*(k-1)+j-k*(k-1)/2,:);
       end
   end
end

%%%%%%%%%%%%define l_j, j\neq i%%%%%%%%%
%l=zeros(N,N);
l = inv(eye(N) - G);

%%%%%equality constraint Aeq*x = b_eq
A_eq = zeros(N^2,n_dim);
b_eq = zeros(N^2,1);
for m=1:N %m-th row of \mathcal{H}
    for j=1:N %j-th column of l, l_j
        if j~= i
            for k=1:N
                A_eq((m-1)*N+j,:) =  A_eq((m-1)*N+j,:)+ l(k,j)*Q(k+(m-1)*N,:);
            end
        end
    end
   
end


%%%%%inquality constraint A*x <= b
A = zeros(2,n_dim);
b = zeros(2,1);

for j=1:N
   A(1,:) =  A(1,:)- I_NN(N*(j-1)+1-(j-1)*(j-2)/2,:);
   A(2,:) = -A(1,:);
end
b=[-0.5*h_i'*h_i;1.5*h_i'*h_i]; %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = quadprog(2*eye(n_dim),f,A,b,A_eq,b_eq);

H_new = vec_matrix(N,x);

end