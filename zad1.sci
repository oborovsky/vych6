function x=shufle(m, d)
    a = m(1,:);
    b = m(2,:);
    c = m(3,:);
    s = length(b);
    et(1) = 0;
    ks(1) = 0;
    
    for i=1:s
        ks(i+1) = (-c(i)) / (a(i)*ks(i) + b(i));
        et(i+1) = (d(i) - a(i) * et(i)) / (a(i) * ks(i) + b(i));
    end
    
    x(s+1) = 0;
    
    for j=0:s-1
        x(s-j) = ks(s-j+1)*x(s-j+1) + et(s-j+1);
    end
    x(s+1) = [];
endfunction

function [m, d] = makeLinearSystem (u, t, h)
    a = [0];
    b = [1];
    c = [0];
    dd = [0];
    s = length(u);
    
    for i=2:s-1
        a(i) = 1/(6*t) - u(i)/(4*h);
        b(i) = 2/(3*t);
        c(i) = 1/(6*t) + u(i)/(4*h);    
        dd(i) = (u(i+1) + 4*u(i) + u(i-1)) / (6*t) + u(i)*(u(i-1) - u(i+1))/(4*h);
    end

    a(s) = 1/(2*t) - u(s)/(2*h);
    b(s) = 1/(2*t) + u(s)/(2*h);
    c(s) = 0;
    dd(s) = (u(s-1) + u(s))/(2*t) + u(s)*(u(s-1) - u(s))/(2*h);
    
    m = [a';b';c'];
    d = dd';
endfunction

function [e,u,y] = makeApp (T, X)
    t = 1/T;
    h = 1/X;
    u(1,:) = 0:h:1;
    y(1,:) = u(1,:);
    for i=1:T
        [m,d] = makeLinearSystem(u(i,:), t, h);
        x = shufle(m, d);
        u(i+1,:) =  x';
        y(i+1,:) = u(1,:)/(1+t*i);
    end
    e = y - u;
endfunction

N = 5;
[e,u,y] = makeApp(N^2, N);
ee(1) = abs(max(e));
for i = 2:2
    N = 2*N;
    X = N;
    T = X^2;
    [e,u,y] = makeApp(T, X);
    ee(i) = abs(max(e));
    x = 0:1/X:1;
   
    printf("ee(%d)=%f,ee(%d)=%f, p = %f\n",i-1,ee(i-1),i,ee(i), abs(log2(ee(i-1)/ee(i))));
end



