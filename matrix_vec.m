function x = matrix_vec(N,H)

n_dim = N*(N+1)/2;
x = zeros(n_dim,1);


for i=1:N
   for j = i:N
       x(N*(i-1)+j-i*(i-1)/2,1) =H(i,j); %i-th row j-th column
   end
end




end