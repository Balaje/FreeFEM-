%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -d/dx(du/dx)-2u*du/dx = 0                         %
%           u(0) = 1, u(1) = 1.5                     %   
%  Using FEM to solve the non-linear problem         %   
%   Exact Solution = 1/(1+x);                        %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 25; % 25 Elements

a = 0; b = 1;

h = (b-a)/N;
x = a:h:b;
K = zeros(N+1,N+1); % Global Stiffness Matrix
J = zeros(N+1,N+1); % Global Jacobian Matrix
F = zeros(N+1,1);   % Global Force Vector

% Use Newton Raphson Method to solve the non-linear system of equations

error = 1000;
tol = 10^-8;
U = zeros(N+1,1);
while error > tol
    % Assembling the local stiffness and local jacobian matrices  
    for j=1:N
        K([j j+1],[j j+1]) = K([j j+1],[j j+1]) + ...
            [1/h+2*U(j)/3+2*U(j+1)/6, -1/h-2*U(j)/3-2*U(j+1)/6 ; ... 
             -1/h+2*U(j)/6+2*U(j+1)/3, 1/h-2*U(j)/6-2*U(j+1)/3]; % Using Linear Elements
         
        J([j j+1],[j j+1]) = J([j j+1],[j j+1]) + ...
            [1/h+4*U(j)/3-2*U(j+1)/6, -1/h-2*U(j)/3-4*U(j+1)/6 ; ...
             -1/h+4*U(j)/6+2*U(j+1)/6, 1/h+2*U(j)/6-4*U(j+1)/3]; % Using Linear Elements
    end
    U(1) = 1;
    U(N+1) = 0.5;
    % Getting the LHS of F(U) = 0 for U^n+1 = U^n - J^-1*F
    F = K*U;
    % Applying the boundary conditions
    F(1) = 0;
    F(N+1) = 0;
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

plot(x,U,'*',x,exact);
legend('Approx','Exact');
