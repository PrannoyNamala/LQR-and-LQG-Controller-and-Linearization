%% Project 2
% Two pendulums on cart
% Shelly Bagchi & Prannoy Namala
clc; clear; close all;

%% Solve nonlinear system
%syms x theta1 theta2 
syms m1 m2 M l1 l2 F %g
%g=-9.81; %M=1000; m1=100;m2=100; l1=20; l2=10;
g=10;

fun = @doublependulum;
X0 = [0,0,0,0,0,0];
%X = fsolve(fun,X0)
%[t,X] = ode45(fun, [0 60], X0)


%% Linearized system
% state:  x = [x x-dot theta1 theta1-dot theta2 theta2-dot]' 
A = [0 1 0 0 0 0; 
     0 0 -g*m1/M 0 -g*m2/M 0;
     0 0 0 1 0 0;
     0 0 -g/l1*(1+m1/M) 0 -g/l1*(m2/M) 0;
     0 0 0 0 0 1;
     0 0 -g/l2*(m1/M) 0 -g/l2*(1+m2/M) 0 ];

B = [0;1/M;0;1/(l1*M);0;1/(l2*M)];
C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];

% Check symbolic controllability matrix
co = simplify( [B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B] )
uncontrollable = length(A) - rank(co)  % num uncontrollable states: 0=CC

%% Determine value ranges for M,m1,m2,l1,l2
cc = simplify(det(co))
solve(cc==0, M,m1,m2,l1,l2, 'ReturnConditions', true);
[ans.M ans.m1 ans.m2 ans.l1 ans.l2]

% Equate 1st&3rd, 2nd&4th columns (nonzero terms ignored)
%solve(co(2,1)-co(2,3),co(4,1)-co(4,3),co(6,1)-co(6,3),...
%      co(1,2)-co(1,4),co(3,2)-co(3,4),co(5,2)-co(5,4),...
%      M,m1,m2,l1,l2, 'ReturnConditions', true)


%% Substitute given values
M=1000; m1=100;m2=100; l1=20; l2=10; g=10; 
A = double(subs(A)); B = double(subs(B));
% Check stabillity
eA = eig(A)
% Check controllability after substitution 
co = ctrb(A,B) 
cc = det(co)
uncontrollable = length(A) - rank(co) 


%% LQR control 

R=1;
%Q = eye(6);
Q = C'*C;
%G = -B*B';
[K1,P1,e1] = lqr(A,B,Q,R);
% Check eigenvalues for Lyapunov's indirect method (stability)
eAc1 = eig(A-B*K1)
%[P2,K2,e2] = icare(A,B,Q,R,[],[],[])  % should be same as above bc we didn't use the other params
%eig(A-B*K2)

R = 0.1;
Q = C'*C * 100; % Weights on thetas
Q(1,1) = 5000   % Weight on displacement x
[K2,P2,e2] = lqr(A,B,Q,R)
eAc2 = eig(A-B*K2)


%% Simulation of LQR Feedback System
Ac = (A-B*K2);
states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'F'};
outputs = {'x'; 'theta1'; 'theta2'};
sys_ss = ss(A ,B,C,0,'statename',states,'inputname',inputs,'outputname',outputs);
sys_cl = ss(Ac,B,C,0,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.1:300;  % 5 min
F_in = 5*ones(size(t));  % step input - constant force?
[y,t,x]=lsim(sys_ss,F_in,t);
[y_cl,t,x_cl]=lsim(sys_cl,F_in,t);

%% Visualization
figure
yyaxis left
hold on
plot(t,y(:,1), 'g-');
plot(t,y_cl(:,1), '-');
title('Linearized State Response after LQR Control')
xlabel('time (s)')
ylabel('cart position, x (m)')
yyaxis right
plot(t,y(:,2), 'c--');
plot(t,y(:,3), 'c-.');
plot(t,y_cl(:,2), '--');
plot(t,y_cl(:,3), '-.');
ylabel('load angle, \theta (radians)')
legend('x','x_{cl}','\theta_1','\theta_2','\theta_1_{cl}','\theta_2_{cl}')
hold off


%%  Part 2 - Observability
% Check observability of other output vectors
% y=x
%C1 = [1 0 0 0 0 0];
C1 = [1 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0];
ob = obsv(A,C1);
observability = rank(ob)
 % y=theta1,theta2
C2 = [0 0 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 0 1 0];
ob = obsv(A,C2);
observability = rank(ob)
 % y=x,theta2
C3 = [1 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 1 0];
ob = obsv(A,C3);
observability = rank(ob)

% Check observability of original system
% y=x,theta1,theta2
ob = obsv(sys_ss);
observability = rank(ob)


%% Construct Luenberger observer
% Pole placement based on eig(A)
poles = [-1 -2 -3 -4 -5 -6];
L = place(A',C',poles)'

% Construct new estimator system
Ae = [A             zeros(size(A));
    zeros(size(A)) (A-L*C)]; 
Be = [B;
   zeros(size(B))];
Ce = [C zeros(size(C))];
   
states_e = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot' ...
        'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
inputs = {'F'};
outputs = {'x'; 'theta1'; 'theta2'};

sys_est = ss(Ae,Be,Ce,0,'statename',states_e,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
F_in = 1*ones(size(t));  % Unit step input
%[y,t,x]=lsim(sys_ss,F_in,t);
X0 = [0,0,0,0,0,0,1,1,pi,pi,pi,pi];
[y_est,t,x_est]=lsim(sys_est,F_in,t,X0);
%% Visualization of Estimator System
figure
yyaxis left
hold on
%plot(t,y(:,1), 'g-');
plot(t,x_est(:,7), '-');
title('Error of Step Response with Observer-Based State-Feedback Control')
xlabel('time (s)')
ylabel('cart position, x (m)')
yyaxis right
%plot(t,y(:,2), 'c--');
%plot(t,y(:,3), 'c-.');
plot(t,x_est(:,9), '--');
plot(t,x_est(:,11), '-.');
ylabel('load angle, \theta (radians)')
%legend('x','x_{e}','\theta_1','\theta_2','\theta_1_{e}','\theta_2_{e}')
legend('x_{e}','\theta_1_{e}','\theta_2_{e}')
hold off

   
%% LQG Feedback Control
% Set up system with smallest output vector x from observability check
C1 = [1 0 0 0 0 0];
states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'F'};
outputs = {'x'};
sys_ss = ss(A ,B,C1,0,'statename',states,'inputname',inputs,'outputname',outputs);
sys_obs = ss(A,B,C1,0,'statename',states,'inputname',inputs,'outputname',outputs);

% Choose gain for cost
QXU = eye(6+1);
% Generate noise
w = wgn(6,1, 5);
v = wgn(1,1, 5);
QWV = [w;v]*[w' v'];

KLQG = lqg(sys_obs,QXU,QWV)

%% Visualization of LQG Feedback 
t = 0:0.01:10;
F_in = 1*ones(size(t));  % Unit step input
[y,t,x]=lsim(sys_ss,F_in,t);  % Get output from original sys
[y_e,t,x_e]=lsim(KLQG,y,t);  % Feed into LQG sys
[y,t,x]=lsim(sys_ss,y_e,t);  % Feed new F back into original sys
%X_dot = A*x;

figure
subplot(2,1, 1)
plot(t,y_e(:,1), '-');
title('LQG Feedback Control Signal')
xlabel('time (s)')
ylabel('new input force, F (N)')

subplot(2,1, 2)
hold on
plot(t,y(:,1), '-');
title('Output from LQG Feedback Control')
xlabel('time (s)')
ylabel('cart position, x (m)')
legend('x')
hold off

figure
%subplot(3,1, 3)
yyaxis left
hold on
plot(t,x(:,2), '-');
%plot(t,y(:,1), '-');
%plot(t,y_e(:,1), '-');
title('Step Response with LQG Feedback Control')
xlabel('time (s)')
ylabel('cart position, x (m)')
yyaxis right
plot(t,x(:,4), '--');
plot(t,x(:,6), '-.');
%plot(t,y(:,2), '--');
%plot(t,y(:,3), '-.');
%plot(t,y_e(:,2), '--');
%plot(t,y_e(:,3), '-.');
ylabel('load angle, \theta (radians)')
%legend('x','x_{e}','\theta_1','\theta_2','\theta_1_{e}','\theta_2_{e}')
legend('x_{f}-dd','\theta_1_{f}-dd','\theta_2_{f}-dd')
hold off


