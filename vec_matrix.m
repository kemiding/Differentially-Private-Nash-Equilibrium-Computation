function H_matrix = vec_matrix(N,x)

H_matrix = zeros(N,N);

%%assign values for upper-diagonal entries
for i=1:N
   for j = i:N
       H_matrix(i,j) = x(j+N*(i-1)-i*(i-1)/2) ; %i-th row j-th column
   end
end


%%assign values for strict-lower-diagonal entries
for i=2:N
   for j=1:i-1
       H_matrix(i,j) = H_matrix(j,i);
   end
end

end