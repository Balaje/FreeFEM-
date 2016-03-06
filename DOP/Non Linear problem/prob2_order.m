%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -u'' + (u')^3 = 0     0<x<1                    %
%     (u+u')(0) = 3/sqrt(2), u'(1) = 0.5          %
%                                                 %
%        Exact u(x) = sqrt(2*(1+x))               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

format long

N = 40;

a = 0; b = 1;

for p=1:5
    h(p) = (b-a)/N;
    x = a:h(p):b;
    K = zeros(N+1,N+1);
    J = zeros(N+1,N+1);
    F = zeros(N+1,1);

    k = [1/h(p) -1/h(p); -1/h(p) 1/h(p)];
    error = 100;
    tol = 10^-10;

    U = zeros(N+1,1);
    RHS = zeros(N+1,1);

    for j=1:N
        K([j j+1],[j j+1]) = K([j j+1],[j j+1]) + k;
    end
    K(1,1) = K(1,1)-1;

    while error > tol
        J = zeros(N+1);
        F = zeros(N+1,1);
        for j=1:N
            J([j j+1],[j j+1]) = J([j j+1],[j j+1]) + ...
                [1/h(p) + 3*(U(j+1)-U(j))^2/(2*h(p)^2), -1/h(p) - 3*(U(j+1)-U(j))^2/(2*h(p)^2);...
                -1/h(p) + 3*(U(j+1)-U(j))^2/(2*h(p)^2), 1/h(p) - 3*(U(j+1)-U(j))^2/(2*h(p)^2)];
            F([j j+1],1) = F([j j+1],1) + [(U(j+1)-U(j))^3/(2*h(p)^2); (U(j+1)-U(j))^3/(2*h(p)^2)];
        end
    
        % Boundary conditons;
        RHS(1) = -3/sqrt(2);
        RHS(N+1) = 0.5;   
        J(1,1) = J(1,1) - 1;
    
        FU = K*U-F-RHS;
    
        DELTA = J\FU;
        U2 = U - DELTA;
        error = max(abs(U2-U)); % Compute the error.
    
        U = U2;
    end
    exact = zeros(N+1,1);
    for j=1:N+1
        exact(j) = sqrt(2*(1+x(j)));
    end
    
    error1(p) = max(abs(U-exact));
    N = N*2;
end

for j=1:p-1
    order(j) = log(error1(j)/error1(j+1))/log(2);
end

order'

figure(1)
plot(x,U,'*',x,exact);
legend('Approx','Exact');
figure(2)
plot(log(h),log(error1));