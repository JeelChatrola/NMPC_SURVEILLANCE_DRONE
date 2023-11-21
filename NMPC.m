% This code is for RBE 502: Robot Control class's final project (at worcester polytechnic institute), which
% involves creating control algorithm for surveilliance drone robot which
% catches the enemy drone and brings it to the base....

% clearing all previous stuff
clear all; clc; close all;

% adding path of casadi, optimization library
addpath('/home/devsonni/Desktop/casadi/')
import casadi.*

%% Quadrotor Dynamics

T = 0.2; % time stamp
N = 20;  % prediction horizon

% constraints on controls and states

% constraints controls
u_min = 0;
u_max = 3;

% constraints on state
x_min = -5;
x_max = 5;

% creating state vector of the quadrotor

% position
x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z'); 
phi = SX.sym('phi'); theta = SX.sym('theta'); psi = SX.sym('psi'); 

% velocities
v1 = SX.sym('v1'); v2 = SX.sym('v2'); v3 = SX.sym('v3'); 

% angular velocities
omega_1 = SX.sym('omega_1'); omega_2 = SX.sym('omega_2'); 
omega_3 = SX.sym('omega_3');

% states of quadrotor
states_q = [x; y; z; phi; theta; psi; v1; v2; v3; omega_1; omega_2; omega_3];      

% length of state vector
n_states_q = length(states_q);

% control vector of the quadrotor
u1 = SX.sym('u1'); u2 = SX.sym('u2'); u3 = SX.sym('u3'); u4 = SX.sym('u4');
controls_q = [u1; u2; u3; u4]; 

% length of control vector
n_controls = length(controls_q);


%% Quadrotor Specifications
T_inv = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
        0 , cos(phi), -sin(phi);
        0 , sin(phi)/cos(theta), cos(phi)/cos(theta)];

Rot_CE =[ cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
      cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
               -sin(theta),                                 sin(phi)*cos(theta),                                 cos(phi)*cos(theta)];


l = 0.2;
m = 0.5;
I11 = 1.24;
I22 = 1.24;
I33 = 2.48;
g = 9.81;
sigma = 0.01;

I = diag([I11, I22, I33]);

r1 = SX.sym('r1'); r2 = SX.sym('r2'); r3 = SX.sym('r3'); 
n1 = SX.sym('n1'); n2 = SX.sym('n2'); n3 = SX.sym('n3');

r = [0; 0; 0];
n = [0; 0; 0];

%% Back to Quadrotor Dynamics 

% right hand side of the system (main chunk of the system)
rhs_u = [v1; v2; v3; 
    T_inv*[omega_1; omega_2; omega_3]; 
    [0; 0; -g] + (1/m)* Rot_CE * ([0; 0; u1 + u2 + u3 + u4] +  r );
    % last 3 variables
    Rot_CE*(I\([(u2 - u4)*l; (u3 - u1)*l;
    (u1 - u2 + u3 - u4)*sigma] + n - cross([omega_1; omega_2; omega_3], I * [omega_1; omega_2; omega_3])))];

disp(rhs_u)

f_u = Function('f_u', {states_q, controls_q}, {rhs_u});                    % Nonlinear Mapping Function f(x,u)
U = SX.sym('U', n_controls, N);                                            % Desition Variables 
P = SX.sym('P', n_states_q + 3);                                            

% This has the initial states of the quadrotor with the reference states of
% the quadrotor, here reference state is UAV's state that's needed to be
% captured
X = SX.sym('X',n_states_q,(N+1));

% Consists of future predicted states of UAV

% Filling the defined sysem parameters of UAV
X(:,1) = P(1:12); % initial state
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f_u(st,con);
    st_next = st + T*f_value;
    X(:,k+1) = st_next;
end

ff = Function('ff',{U,P},{X});

%% Objective function

obj = 0; % objective function
g = [];  % constrains of pitch angle theta

% this loop will run for prediction horizon times, in this case it's 15
for k=1:N
    stt = X(:, k);
    % obj = obj + ((stt(1) - P(13))^2) + ((stt(2) - P(14))^2) + ((stt(3) - P(15))^2); % squared distance
    % obj = obj + sqrt((stt(1) - P(13))^2) + sqrt((stt(2) - P(14))^2) + sqrt((stt(3) - P(15))^2); % sqrt distance
    obj = obj + sqrt((stt(1) - P(13))^2 + (stt(2) - P(14))^2 + (stt(3) - P(15))^2); % 3D distance, better than upper 
end

% compute the state constrains
 for k=1:N+1
     g = [g, X(1,k)]; % limit on x location
     g = [g, X(2,k)]; % limit on y location
     g = [g, X(3,k)]; % limit on z location
 end

 % disp(size(g))

% make the decision variables one column vector
OPT_variables = reshape(U, 4*N, 1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

% constraints on states (inequality constraints)
% args.lbg = -10;  args.ubg = 10; 
args.lbg(1:3:3*(N+1)) = -5;          args.ubg(1:3:3*(N+1)) = 5;       % bounds on state_x 
args.lbg(2:3:3*(N+1)) = -5;          args.ubg(2:3:3*(N+1)) = 5;       % bounds on state_y
args.lbg(3:3:3*(N+1)) = 0;           args.ubg(3:3:3*(N+1)) = 10;      % bounds on state_z

% constraints on controls (equality constraints)
args.lbx(1:n_controls:n_controls*N,1) = 0;               args.ubx(1:n_controls:n_controls*N,1) = u_max;  
args.lbx(2:n_controls:n_controls*N,1) = 0;               args.ubx(2:n_controls:n_controls*N,1) = u_max;  
args.lbx(3:n_controls:n_controls*N,1) = 0;               args.ubx(3:n_controls:n_controls*N,1) = u_max;
args.lbx(4:n_controls:n_controls*N,1) = 0;               args.ubx(4:n_controls:n_controls*N,1) = u_max;

%% Simulation starts from here

t0 = 0;
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % quadrotor's initial state
xs = [5; 0; 10];                            % reference state
xx(:,1) = x0;                               % storing location of 
% the quadrotor traveled
t(1) = t0;
u0 = zeros(N,4);                            % initial control vector


% NMPC starts form here
mpciter = 0;
xx1 = [];
xss = [];
ss = [];
ss(:,1) = xs;
sc = 0;
loop_run = 100;
u_cl = zeros(n_controls, loop_run);

main_loop = tic;
% while (sqrt((xs(1)-x0(1))^2 + (xs(1)-x0(1))^2 + (xs(1)-x0(1))^2) > 0.000000001)
for i=1:loop_run
    args.p = [x0; xs];
    args.x0 = reshape(u0',n_controls*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)', n_controls, N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution
    xx1(:,1:12,mpciter+1) = full(ff_value)';
    xss(:,1:3,mpciter+1) = full(xs)'; 
    u_cl(:, i) = u(1,:)';
    t(mpciter+1) = t0;
    sc = sc+1;
    [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, xs, sc); % get the initialization of the next optimization step
    
    ss(:,mpciter+2) = xs;
    xx(:,mpciter+2) = x0;
    mpciter
    mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop)

% loop_run = mpciter
%% Plotting Stuff

% just for  plotting 
ss1 = [];
ss1(1,1:loop_run) = 0;

figure
plot3(xx(1,1:loop_run), xx(2,1:loop_run), xx(3,1:loop_run),'b-');
xlim([-10 10]);
ylim([-10 10]);
zlim([0 10]);
xlabel('X[m]'); ylabel('Y[m]'); zlabel('Z[m]');
legend('Quadrotor');
grid on;

figure
plot(u_cl(1,:));

figure
plot(u_cl(2,:));

figure
plot(u_cl(3,:));

figure
plot(u_cl(4,:));

%{
Code for saving video file of simulation
figh = figure
for i=1:700
    hold on;
    plot3(ss(1,i),ss(2,i), ss1(1,i),xx(1,i),xx(2,i), xx(3,i),'go', 'LineWidth', 4, 'MarkerSize', 1);
    plot3(ss(1,i),ss(2,i), ss1(1,i), 'bo', 'LineWidth', 4, 'MarkerSize', 1);
    xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
    legend('Target','UAV');
    view(-18, 57);
    axis([0 1800 50 350 0 200])
    grid on;
    hold off;
    drawnow;
    
    movieVector(i) = getframe(figh);    
end


myWriter = VideoWriter('Track_1');
myWriter.FrameRate = 10;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}