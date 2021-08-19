% https://www.mathworks.com/help/optim/ug/fsolve.html
function F = doublependulum(x)
%syms F_in %g m1 m2 M
g=10; M=1000; m1=100;m2=100; l1=20; l2=10;
F_in = 5;  % test value

% state:  x = [x x-dot theta1 theta1-dot theta2 theta2-dot]' 
N1 = M+m1*sin(x(3)).^2+m2*sin(x(5)).^2;
N2 = F_in - g*(m1*sin(x(3))*cos(x(3)) + m2*sin(x(5))*cos(x(5))) + ...
    m1*l1*sin(x(3))*x(4).^2 + m2*l2*sin(x(5))*x(6).^2;
F(1) = ( 1/N1 )*N2;
F(2) = (1/l1)*(-g*sin(x(3)) + (N2*cos(x(3)))/N1 );
F(3) = (1/l2)*(-g*sin(x(5)) + (N2*cos(x(5)))/N1 );

F(1) = subs(F(1));
F(2) = subs(F(2));
F(3) = subs(F(3));

%fun = @doublependulum;
%x0 = [0,0];
%x = fsolve(fun,x0)