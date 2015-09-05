% Recursive least squares
function [thea P] = myRLS(x,y,thea,P,lamda,gamma,epsilon)

N = numel(thea);
K = P*x/(lamda + x'*P*x);
% if gamma equals to 0, this is regular RLS, otherwise sparse RLS
spa = gamma*(1 - 1/lamda)*(eye(N) - K*x')*P*(sign(thea)./(abs(thea)+epsilon));
thea = thea + K*(y - thea'*x)' + spa;
P = 1/lamda*P - K*x'*(1/lamda)*P;

end


