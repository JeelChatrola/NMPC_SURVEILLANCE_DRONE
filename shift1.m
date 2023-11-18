function [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, xs, sc)
st = x0;
con = u(1,:)';
f_value = f_u(st,con);
st = st + (T*f_value);
x0 = full(st);

xs = full(xs);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end