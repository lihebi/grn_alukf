% The function updates the states fixing the parameter
function [x_next] = fixed_f(x,dt,n_gene,n_pro,a_hat,b_hat,c,d,e)
% ode_f2 is the augmented function used in ukf, and it is a discrete model
%diffun_f2(x,n_gene,n_pro,a_hat,b_hat,c,d,e)
dxdt = @(t,x) diffun_f2(t,x,n_gene,n_pro,a_hat,b_hat,c,d,e);
sol = ode45(dxdt,[0:dt],x,[]);
x_next = deval(sol,dt);
%dt*diffun_f2(x,n_gene,n_pro,a_hat,b_hat,c,d,e);

end

function [xdot] = diffun_f2(t,x,n_gene,n_pro,a_hat,b_hat,c,d,e)
num_state = n_gene + n_pro;
%num_s = n_pro*(n_pro + 1)/2 + 1 ;
% 1~num_state is the true state, that is the gene and protein
% after num_state, are parameters, a, b
%     [a_hat,a_pad] = vec2mat(x(num_state + 1 : num_state + n_pro*num_s),num_s);
%     [b_hat,b_pad] = vec2mat(x(num_state + n_pro*num_s+  1 : end), num_s-1);
s = pro_cpx(x(n_gene + 1: n_gene + n_pro));
% degradtion rate is determined
xdot1 = a_hat*s./([ones(n_gene,1) b_hat]*s) - c.*x(1:n_gene);
xdot2 = d.*x(1:n_gene) - e.*x(n_gene+1:num_state);
xdot = [xdot1; xdot2];
%xdot(numel(x)) = 0;
end