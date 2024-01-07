function [t0, x0, u0, xs, xt, sc, pri, pri2, r, n, distance] = shift1(T, t0, x0, u, f_u, xs, sc, xt, pri, pri2, r, n)
st = x0;
con = u(1,:)';
f_value = f_u(st,con);
st = st + (T*f_value);
x0 = full(st);

xs = full(xs);

% speeds, low, medium, high
% H = 2;
% L = 0.5;
% M = 1;

% v = randi([1 10])*randi([1 7])/100;
v = 0.5;
omega_2_u = randi([1 2])*randi([1 2])/100;
omega_3_u = randi([1 2])*randi([1 2])/100;

distance = sqrt((x0(1)-xt(1))^2 + (x0(2)-xt(2))^2 + (x0(3)-xt(3))^2);

if distance < 0.1
    sc = 1;
    if pri == 0
        disp('Enemy Got Caught Bringing It To The Base');
        pri = 1;
        while vecnorm(r)>1 && vecnorm(n)>0.5
            r = [randi([1 2])*randi([1 2])/4; randi([1 2])*randi([1 2])/4; randi([1 2])*randi([1 2])/5];
            n = [randi([1 2])*randi([1 2])/20; randi([1 2])*randi([1 2])/20; randi([1 2])*randi([1 2])/20];
        end
        disp(vecnorm(r));
        disp(vecnorm(n));
    end
     x0(13:15) = r;
     x0(16:18) = n;
end

if xt(1) < -5 || xt(1) > 5 || xt(2) < -5 || xt(2) > 5
    if pri2 == 0
        disp('Enemy Escaped Returning To The Base');
        pri2 = 1;
    end
    pri = 1;
end

if sc == 1
    xs = [0; 0; 0; 0; 0];
    xt = [x0(1); x0(2); x0(3); 0; 0];

elseif pri == 1 && sc == 0

    xs = [0; 0; 0; 0; 0];

    x = xt(1) + v*cos(xt(5))*cos(xt(4))*T;
    y = xt(2) + v*sin(xt(5))*cos(xt(4))*T;
    z = xt(3) + v*sin(xt(4))*T;
    theta_u = xt(4) + omega_2_u*T;
    psi_u = xt(5) + omega_3_u*T;

    xt = [x; y; z; theta_u; psi_u];

else
    x = xs(1) + v*cos(xs(5))*cos(xs(4))*T;
    y = xs(2) + v*sin(xs(5))*cos(xs(4))*T;
    z = xs(3) + v*sin(xs(4))*T;
    theta_u = xs(4) + omega_2_u*T;
    psi_u = xs(5) + omega_3_u*T;

    xs = [x; y; z; theta_u; psi_u];

    % x = xt(1) + v*cos(xt(5))*cos(xt(4))*T;
    % y = xt(2) + v*sin(xt(5))*cos(xt(4))*T;
    % z = xt(3) + v*sin(xt(4))*T;
    % theta_u = xt(4) + omega_2_u*T;
    % psi_u = xt(5) + omega_3_u*T;

    xt = [x; y; z; theta_u; psi_u];
end

% x = xs(1) + v*cos(xs(5))*cos(xs(4))*T;
% y = xs(2) + v*sin(xs(5))*cos(xs(4))*T;
% z = xs(3) + v*sin(xs(4))*T;
% theta_u = xs(4) + omega_2_u*T;
% psi_u = xs(5) + omega_3_u*T;
% 
% xs = [x; y; z; theta_u; psi_u];

t0 = t0 + T;
u0 = [u(2:size(u,1),:); u(size(u,1),:)];
end