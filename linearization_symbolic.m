clc;
clear all;
clear;

% % Define symbolic variables and parameters
% syms z = [x1,x2,x 3,phi,theta,psi,v1,v2,v3,omega1,omega2,omega3]
% syms remaining parameters
% p = [g,l,m,sigma]
% I = [I11,I22,I33]
% n = [n1,n2,n3]
% r = [r1,r2,r3]
% u = [u1,u2,u3,u4];

z = sym('z', [1 12]);
p = sym('p', [1 4]);
u = sym('u', [1 4]);
I = sym('I', [1 3]);
n = sym('n', [1 3]);
r = sym('r', [1 3]);

I = diag([I(1), I(2), I(3)]);


Rot_CE =[ cos(z(5))*cos(z(6)), sin(z(4))*sin(z(5))*cos(z(6)) - cos(z(4))*sin(z(6)), sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6));
      cos(z(5))*sin(z(6)), cos(z(4))*cos(z(6)) + sin(z(4))*sin(z(5))*sin(z(6)), cos(z(4))*sin(z(5))*sin(z(6)) - sin(z(4))*cos(z(6));
               -sin(z(5)),                                 sin(z(4))*cos(z(5)),                                 cos(z(4))*cos(z(5))];

T_inv = [1, sin(z(4))*tan(z(5)), cos(z(4))*tan(z(5));
        0 , cos(z(4)), -sin(z(4));
        0 , sin(z(4))/cos(z(5)), cos(z(4))/cos(z(5))];

% Define the dynamics of the quadrotor system
f = [[z(7);z(8);z(9)];
     T_inv * [z(10);z(11);z(12)];
     -[0,0,p(1)].' + (1/p(3))* Rot_CE * ([0;0; u(1)+u(2)+u(3)+u(4)] +  [r(1); r(2);r(3)]);
     inv(I) * ( [(u(2) - u(4)) *p(2) ; (u(3) - u(1)) * p(2)  ; (u(1) - u(2) + u(3) - u(4)) * p(4) ] + [n(1); n(2);n(3)] - cross([z(10);z(11);z(12)], I * [z(10);z(11);z(12)]))];

disp(f);



%Linearize the system by calculating the Jacobian matrices
A = jacobian(f, z);
B = jacobian(f, u);

disp(A);
disp(B);

%===========================================================================
% Definition of the current equibilibrium state:

% Substitution of params
% Given
% z0 = [x1,x2,x3,0,0,0,0,0,0,0,0,0]
% u0 = [mg/4, mg/4, mg/4, mg/4]
for i = 1:1:4 
    A = subs(A, u(i), p(3)*p(1)/4);
    B = subs(B, u(i), p(3)*p(1)/4);
end

for i = 4:1:12 % As x1,x2,x3 are not given 
    A = subs(A, z(i), 0);
    B = subs(B, z(i), 0);
end

for i = 1:1:3 % As n1,n2,n3 are not given 
    A = subs(A, n(i), 0);
    B = subs(B, n(i), 0);
end

for i = 1:1:3 % As r1,r2,r3 are not given 
    A = subs(A, r(i), 0);
    B = subs(B, r(i), 0);
end


% Assuming there are no external forces and moments on the system

disp(A);
disp(B);



% Substituting some of the parameters

% Taking all in MKG 
%g = 9.8 (Normal Acceleration due to gravity)
%l = 0.450  (450 mm) (Normal Size Drone - 450mm chasis)
%m = 1 (kg) (Normal Weight of commercial drone)
%sigma = 0.95 (Generic thrust proportionality for commercial motors - can
%vary)

p = [9.8,0.450, 1, 0.95];

for i = 1:1:4 % Assumed params
    A = subs(A, p(i));
    B = subs(B, p(i));
end

I = [1,1, 1] * 1*(0.450)^2;

for i = 1:1:3 % Assumed params
    A = subs(A, I(i));
    B = subs(B, I(i));
end


%===========================================================================

disp(A);
disp(B);





% % Substitute the equilibrium values into A and B to obtain the linearized matrices
% A = subs(A, [z; u], [z_e; u_e]);
% B = subs(B, [z;u], [z_e; u_e]);
% 
% disp(size(A));
% disp(size(B));
% 
% disp(size(A));
% disp(size(B));
% 
% 
% Convert symbolic expressions to numeric values
A_numeric = double(A);
B_numeric = double(B);



disp(size(A_numeric));
disp(size(B_numeric));

% With feedback
% % % Define C and D matrices based on the output equations (modify as needed)
C = eye(12);  % Assuming direct measurement of all states
D = zeros(12, 4);  % No direct feedthrough
% 
% % Create the state-space system
sys = ss(A_numeric, B_numeric, C, D);
 
disp(sys)

% % Now, sys represents the linearized quadrotor system near the equilibrium point.

Ctrb = ctrb(A, B);
rank_Ctrb = rank(Ctrb);


disp(rank_Ctrb)

num_state_variables = 12;

if rank_Ctrb == num_state_variables
    disp('The linearized quadrotor system is controllable.');
else
    disp('The linearized quadrotor system is not fully controllable.');
end



