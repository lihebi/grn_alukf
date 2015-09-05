function [y_res] = ode_h(x,len,n)

B = eye(len);
if len < n
    B(len,n)=0;
end
y_res = B*x;
end