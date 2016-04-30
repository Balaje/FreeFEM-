%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -d/dx(du/dx)-2u*du/dx = 0                         %
%           u(0) = 1, u(1) = 1.5                     %   
%  Using FEM to solve the non-linear problem         %   
%   Exact Solution = 1/(1+x);                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 3; % 2 Elements

a = 0; b = 1;

for p=1:3
    h(p) = (b-a)/N;
    x = a:h(p):b;
    K = zeros(N+1,N+1); % Global Stiffness Matrix
    J = zeros(N+1,N+1); % Global Jacobian Matrix
    F = zeros(N+1,1);   % Global Force Vector

    % Use Newton Raphson Method to solve the non-linear system of equations

    error = 1000;
    tol = 10^-10;
    U = zeros(N+1,1);
    while error > tol
        % Assembling the local stiffness and local jacobian matrices
        J = zeros(N+1);
        K = zeros(N+1);
        for j=1:N
            K([j j+1],[j j+1]) = K([j j+1],[j j+1]) + ...
                [1/h(p)+2*U(j)/3+2*U(j+1)/6, -1/h(p)-2*U(j)/3-2*U(j+1)/6 ; ... 
                -1/h(p)+2*U(j)/6+2*U(j+1)/3, 1/h(p)-2*U(j)/6-2*U(j+1)/3]; % Using Linear Elements
         
            J([j j+1],[j j+1]) = J([j j+1],[j j+1]) + ...
                [1/h(p)+4*U(j)/3-2*U(j+1)/6, -1/h(p)-2*U(j)/3-4*U(j+1)/6 ; ...
                -1/h(p)+4*U(j)/6+2*U(j+1)/6, 1/h(p)+2*U(j)/6-4*U(j+1)/3]; % Using Linear Elements
        end
    
        % Getting the LHS of F(U) = 0 for U^n+1 = U^n - J^-1*F
        F = K*U;
        % Applying the boundary conditions
        F(1) = U(1)-1;
        F(N+1) = U(N+1)-0.5;
        J(1,1) = 1;
        J(1,2) = 0;
        J(N+1,N) = 0;
        J(N+1,N+1) = 1;
    
        % Performing the Newton-Raphson Step
        DELTA = J\F;
        U2 = U - DELTA;
        error = max(abs(U2-U)); % Compute the error.
    
        U = U2;
    end
    exact = zeros(N+1,1);
    for j=1:N+1
        exact(j) = 1/(1+x(j));
    end
    
    error1(p) = max(abs(U-exact));
    N=N*2;
end

for j=1:p-1
    order(j) = log(error1(j)/error1(j+1))/log(2);
end

order'

figure(1)
plot(x,U,'*',x,exact);
legend('Exact Solution','Approximate Solution')
grid on
xlabel('x');
ylabel('u');
str = strcat(num2str(size(U,1)-1),' elements');
title(str)

figure(2)
plot(log(h),log(error1));
grid on
xlabel('log(h)');
ylabel('log(error)');
title('Order of convergence plot');

