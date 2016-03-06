function U = FEM1D(x)
    M = length(x);
    h = zeros(M,1);
    for i=1:M-1
        h(i) = x(i+1)-x(i);            
    end
    
    % Initializing
    A = sparse(M,M);
    A(1,1) = 1; F(1) = 0;
    A(M,M) = 1; F(M) = 0;
    A(2,2) = 1/h(1);
    F(2) = int_f_hat2(x(1),x(2));
    
    for j=2:M-2
        A(j,j) = A(j,j) + 1/h(j);
        A(j+1,j) = A(j+1,j) - 1/h(j);
        A(j,j+1) = A(j,j+1) - 1/h(j);
        A(j+1,j+1) = A(j+1,j+1) + 1/h(j);
        
        F(j) = F(j) + int_f_hat1(x(j),x(j+1));
        F(j+1) = F(j+1) + int_f_hat2(x(j),x(j+1));
    end
    
    A(M-1,M-1) = A(M-1,M-1) + 1/h(M-1);
    F(M-1) = F(M-1) + int_f_hat2(x(M-1),x(M));
    
    U = A\F';
    error = U - u_exact(x)';
    
    %l2norm_of_error = norm(error,2)
end


