% This code is for RBE 502: Robot Control class's final project (at worcester polytechnic institute), which
% involves creating control algorithm for surveilliance drone robot which
% catches the enemy drone and brings it to the base....

% clearing all previous stuff
clear all; clc; close all;

% adding path of casadi, optimization library
addpath('casadi/')
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

% forces & disturbances
r1 = SX.sym('r1'); r2 = SX.sym('r2'); r3 = SX.sym('r3'); 
n1 = SX.sym('n1'); n2 = SX.sym('n2'); n3 = SX.sym('n3');

% states of quadrotor
states_q = [x; y; z; phi; theta; psi; v1; v2; v3; omega_1; omega_2; omega_3; r1; r2; r3; n1; n2; n3];      

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

%% Back to Quadrotor Dynamics 

U = SX.sym('U', n_controls, N);                                            % Desition Variables 
P = SX.sym('P', n_states_q + 5);                                            


% right hand side of the system (main chunk of the system)
rhs_u = [v1; v2; v3; 
    T_inv*[omega_1; omega_2; omega_3]; 
    [0; 0; -g] + (1/m)* Rot_CE * ([0; 0; u1 + u2 + u3 + u4] +  [r1; r2; r3]);
    % last 3 variables
    Rot_CE*(I\([(u2 - u4)*l; (u3 - u1)*l;
    (u1 - u2 + u3 - u4)*sigma] + [n1; n2; n3] - cross([omega_1; omega_2; omega_3], I * [omega_1; omega_2; omega_3])));
    0;
    0;
    0; 
    0;
    0; 
    0];

f_u = Function('f_u', {states_q, controls_q}, {rhs_u});                    % Nonlinear Mapping Function f(x,u)

% This has the initial states of the quadrotor with the reference states of
% the quadrotor, here reference state is UAV's state that's needed to be
% captured
X = SX.sym('X',n_states_q,(N+1));

% Consists of future predicted states of UAV

% Filling the defined sysem parameters of UAV
X(:,1) = P(1:18); % initial state
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
    obj = obj + 100*sqrt((stt(1) - P(19))^2 + (stt(2) - P(20))^2 + (stt(3) - P(21))^2); % 3D distance, better than all
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
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % quadrotor's initial state

x_en = randi([-5 5]);
y_en = randi([-5 5]);
z_en = randi([3 9]);
plane = 3;

if plane == 1
    % disp("Enemy UAV Entering From Plane " + plane);
    x_en = -5;
    head = 0;
    % head = (randi([0 3])*randi([1 2])/4);
elseif plane == 4
    % disp("Enemy UAV Entering From Plane " + plane);
    y_en = 5;
    head = -(randi([0 3])*randi([1 2])/4);
elseif plane == 3
    % disp("Enemy UAV Entering From Plane " + plane);
    x_en = 5;
    head = -pi;
elseif plane == 2
    % disp("Enemy UAV Entering From Plane " + plane);
    y_en = -5;
    % head = 0;
    head = (randi([0 3])*randi([1 2])/4);
end

xs = [x_en; y_en; z_en; 0; head];                % reference state

xt = [x_en; y_en; z_en; 0; head];                % back to base
                              
% the quadrotor traveled
t(1) = t0;
u0 = zeros(N,4);                            % initial control vector


% NMPC starts form here
mpciter = 0;
xx1 = [];
xss = [];
ss = [];
tt = [];
dis = [];
ss(:,1) = xs;
xx(:,1) = x0; 
tt(:,1) = xt;
sc = 0;
loop_run = 70;
u_cl = zeros(n_controls, loop_run);
pri = 0;
pri2 = 0;
base = 10;
r = [5; 5; 5];
n = [5; 5; 5];

main_loop = tic;
while base > 0.0001
    args.p = [x0; xs];
    args.x0 = reshape(u0',n_controls*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)', n_controls, N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution
    xx1(:,1:18,mpciter+1) = full(ff_value)';
    xss(:,1:5,mpciter+1) = full(xs)'; 
    u_cl(:, mpciter+1) = u(1,:)';
    t(mpciter+1) = t0;

    [t0, x0, u0, xs, xt, sc, pri, pri2, r, n, distance] = shift1(T, t0, x0, u, f_u, xs, sc, xt, pri, pri2, r, n); % get the initialization of the next optimization step
    
    if sc == 1 && pri == 1
        base = sqrt((x0(1)-0)^2 + (x0(2)-0)^2 + (x0(3)-0)^2);
    end
    if sc == 0 && pri == 1
        base = sqrt((x0(1)-0)^2 + (x0(2)-0)^2 + (x0(3)-0)^2);
    end

    if sc == 1 % enemy got caught
        ss(:,mpciter+2) = xt;
    elseif sc == 0 && pri == 1 % enemy escaped
        ss(:, mpciter+2) = xt;
    else % no caught, no escape
        ss(:, mpciter+2) = xs;
    end

    xx(:,mpciter+2) = x0;
    dis(mpciter+1) = distance;

    mpciter
    mpciter = mpciter + 1;
end
main_loop_time = toc(main_loop)

% loop_run = mpciter
%% Plotting Stuff

% just for  plotting 
ss1 = [];
ss1(1,1:loop_run) = 0;
time = linspace(0, 0.2, (mpciter+1));

disp(size(time))
disp(size(u_cl(1,:)));

% %{
figure
% plotting quadrotor trajectory 
plot3(xx(1,1:mpciter+1), xx(2,1:mpciter+1), xx(3,1:mpciter+1), 'b-', 'LineWidth', 2);
hold on;
plot3(ss(1,1:mpciter+1), ss(2,1:mpciter+1), ss(3,1:mpciter+1), 'r-', 'LineWidth', 2);
hold on;
% plotting airspace
plot3([-5, 5], [-5, -5], [10, 10], 'k-','LineWidth', 2);
plot3([-5, 5], [5, 5], [0, 0], 'k-','LineWidth', 2);
plot3([5, 5], [-5, 5], [10, 10], 'k-','LineWidth', 2);
plot3([-5, -5], [-5, 5], [0, 0], 'k-','LineWidth', 2);
plot3([-5, 5], [-5, -5], [0, 0], 'k-','LineWidth', 2);
plot3([-5, 5], [5, 5], [10, 10], 'k-','LineWidth', 2);
plot3([5, 5], [-5, 5], [0, 0], 'k-','LineWidth', 2);
plot3([-5, -5], [-5, 5], [10, 10], 'k-','LineWidth', 2);
plot3([-5, -5], [5, 5], [0, 10], 'k-','LineWidth', 2);
plot3([5, 5], [-5, -5], [0, 10], 'k-','LineWidth', 2);
plot3([-5, -5], [-5, -5], [0, 10], 'k-','LineWidth', 2);
plot3([5, 5], [5, 5], [0, 10], 'k-','LineWidth', 2);

xlim([-10 10]);
ylim([-10 10]);
zlim([0 20]);
xlabel('X[m]'); ylabel('Y[m]'); zlabel('Z[m]');
legend('Quadrotor', 'Enemy UAV', 'AirSpace');
grid on;
%}

figure
plot(dis, 'LineWidth', 2);
xlabel('Iteration'); ylabel('Meter');
title('Distance B/W Quadrotor and UAV')

figure
plot(u_cl(1,:), 'LineWidth', 2);
xlabel('Iterations'); ylabel('u1');
title('Control u1');

figure
plot(u_cl(2,:), 'LineWidth', 2);
xlabel('Iterations'); ylabel('u2');
title('Control u2');

figure
plot(u_cl(3,:), 'LineWidth', 2);
xlabel('Iterations'); ylabel('u3');
title('Control u3');

figure
plot(u_cl(4,:), 'LineWidth', 2);
xlabel('Iterations'); ylabel('u4');
title('Control u4');

N = 10;
Q = linspace(0,2*pi,N)';
circle = 0.3*l*[cos(Q) sin(Q) zeros(N,1)];
loc = l*[1 0 0; 0 1 0; -1 0 0; 0 -1 0];

%{
% Code for saving video file of simulationfigh = figure
figh = figure
for i=1:(mpciter+1)

    plot3(xx(1,1:i), xx(2,1:i), xx(3,1:i),'g-', 'LineWidth', 1.5);
    Rotation_matrix = [cos(xx(5,i))*cos(xx(6,i)), sin(xx(4,i))*sin(xx(5,i))*cos(xx(6,i)) - cos(xx(4,i))*sin(xx(6,i)), sin(xx(4,i))*sin(xx(6,i)) + cos(xx(4,i))*sin(xx(5,i))*cos(xx(6,i));
 cos(xx(5,i))*sin(xx(6,i)), cos(xx(4,i))*cos(xx(6,i)) + sin(xx(4,i))*sin(xx(5,i))*sin(xx(6,i)), cos(xx(4,i))*sin(xx(5,i))*sin(xx(6,i)) - sin(xx(4,i))*cos(xx(6,i));
 -sin(xx(5,i)), sin(xx(4,i))*cos(xx(5,i)), cos(xx(4,i))*cos(xx(5,i))];
    for j=1:4
        ctr(j,:) = xx(1:3,i)' + loc(j,:)*Rotation_matrix';
        pose = ones(1,1)*xx(1:3,i)' + (ones(1,1)*loc(j,:) + circle)*Rotation_matrix';
        if i > 1
            % Clear the previous pose for the current j
            delete(h_poses{j});
        end
        % Plot the updated pose
        h_poses{j} = plot3(pose(:,1), pose(:,2), pose(:,3), 'b-', 'LineWidth', 2);
    end
    

    hold on;
    plot3(ss(1,1:i), ss(2,1:i), ss(3,1:i), 'r-', 'LineWidth', 1.5);
     % Add a sphere at the location specified by "ss"
    hold on;
    sphere_radius = 0.1    ; % Adjust the radius as needed
    [x_sphere, y_sphere, z_sphere] = sphere;
    h_sphere = surf(x_sphere * sphere_radius + ss(1,i), y_sphere * sphere_radius + ss(2,i), z_sphere * sphere_radius + ss(3,i));
    set(h_sphere, 'FaceColor', 'blue', 'FaceAlpha', 0.5);
    
    hold on;
    % plotting airspace
    plot3([-5, 5], [-5, -5], [10, 10], 'k-','LineWidth', 1);
    plot3([-5, 5], [5, 5], [0, 0], 'k-','LineWidth', 1);
    plot3([5, 5], [-5, 5], [10, 10], 'k-','LineWidth', 1);
    plot3([-5, -5], [-5, 5], [0, 0], 'k-','LineWidth', 1);
    plot3([-5, 5], [-5, -5], [0, 0], 'k-','LineWidth', 1);
    plot3([-5, 5], [5, 5], [10, 10], 'k-','LineWidth', 1);
    plot3([5, 5], [-5, 5], [0, 0], 'k-','LineWidth', 1);
    plot3([-5, -5], [-5, 5], [10, 10], 'k-','LineWidth', 1);
    plot3([-5, -5], [5, 5], [0, 10], 'k-','LineWidth', 1);
    plot3([5, 5], [-5, -5], [0, 10], 'k-','LineWidth', 1);
    plot3([-5, -5], [-5, -5], [0, 10], 'k-','LineWidth', 1);
    plot3([5, 5], [5, 5], [0, 10], 'k-','LineWidth', 1);
    % hold on;

    xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
    legend('Quadrotor', 'Enemy UAV', 'AirSpace');
    view(-18*(i/10), 50); %57
    axis([-10 10 -10 10 0 20])
    grid on;
    drawnow;
    % pause(0.5);
    x0 = 0;
    y0 = 0;
    width = 1920;
    height = 1080;
    set(gcf,'position', [x0,y0,width,height])
    movieVector(i) = getframe(gcf);    
    delete(h_sphere);
    % hold off;
end


myWriter = VideoWriter('test_e1', 'Uncompressed AVI');
myWriter.FrameRate = 5;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}