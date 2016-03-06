function V = order()
% Script to determine the order of convergence
N = 10;
a = 0; b = 1; % Line interval
h = zeros(4,1);
error = h;
for i=1:10
    h(i) = (b-a)/N;
    x = a:h(i):b;
    
    U_app = FEM1D(x);
    U_exact = u_exact(x);
    
    error(i) = norm((U_app-U_exact'),4);
    N=N*2;
end

for p=1:i-1
    V(p) = log(error(p)/error(p+1))/log(2);
end
plot(log(h),log(error));

end