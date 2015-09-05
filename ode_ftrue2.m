% function f_ture in ode model
function [x_next] = ode_ftrue2(x,dt,n_gene,n_pro,a,b,c,d,e)
% the model used to generate simulated data which is discrete here
dxdt = @(t,x) diffun_ftrue2(t,x,n_gene,n_pro,a,b,c,d,e);
sol = ode45(dxdt,[0:dt],x,[]);
x_next = deval(sol,dt);
%x_next = x + dt*diffun_ftrue2(x,n_gene,n_pro,a,b,c,d,e);
end

function [xdot] = diffun_ftrue2(t,x,n_gene,n_pro,a,b,c,d,e)
s = pro_cpx(x(n_gene + 1: n_gene + n_pro));
xdot1 = a*s./([ones(n_gene,1) b]*s) - c.*x(1:n_gene);
xdot2 = d.*x(1:n_gene) - e.*x(n_gene + 1 : n_gene + n_pro);
xdot = [xdot1; xdot2];
%    xdot(numel(x)) = 0;
end


