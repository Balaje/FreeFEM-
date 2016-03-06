function V = int_f_hat1(x1,x2)
    N = 100; % Number of points for simpson's rule.
    h = (x2-x1)/N; % Step size
    x = x1:h:x2; % Partition interval [x1,x2]
    
    sum = f(x(1))*hat1(x(1),x1,x2);
    for j=1:N/2-1
        sum = sum + f(x(2*j))*hat1(x(2*j),x1,x2)...
            + 4*f(x(2*j+1))*hat1(x(2*j+1),x1,x2) +...
            f(x(2*j+2))*hat1(x(2*j+2),x1,x2);
    end
    
    V = h/3*sum;
end