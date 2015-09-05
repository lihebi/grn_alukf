clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model is as follows(x_dot is the first derivative of the x)
%           x_dot = f(x) + w
%           z   = h(x) + v
% Specifically, x = [x1,x2,y1,y2] is the state, z is the observation of gene level
%           x1_dot = (0.01 + 0.1*y2)/(1 + 10*y1 + 0.1*y2) - 0.1*x1
%           x2_dot = (0.1)/(1 + y1*y2) - 0.1*x2
%           y1_dot = 0.3*x1  - 0.5*y1
%           y2_dot = 0.5*x2 - 0.5*y2
%           z   = [x1;x2;y1;y2] + v
% v ~ N(0,R), Gaussian noise is added in the observation model

n_gene = 2;
n_pro = 2;
n_state = n_gene + n_pro;
n_cpx = n_pro*(n_pro + 1)/2 + 1; % protein and its complex
n = n_state;%2*n_cpx*n_gene + n_state -2; % argumented state dimension
q =0;
Q  = q^2*eye(n_state); % Covariance of noise state equation
if n_state > n
    Q(n,n) = 0;
end
r = 0;
R = r^2*eye(n_state); % !!!!!Covariance of noise measurement equation
Q_x = Q(1:n_state,1:n_state); % Covariance of state except for parameter
a = [0.01 0 0.1 0; 0.1 0 0 0];
b = [ 10 0.1 0;  0 0 1];
ab_true = [a(1,:) b(1,:) a(2,:) b(2,:)]';
d = [0.3 0.5]';
c = [0.1 0.1]';
e = [0.5 0.5]';
T = 40;% length of time
dt = 1;% sampling interval
N = T/dt + 1; %step size
%generate true data
% x_true = zeros(n_state,N);
%
% x_true(:,1) = ones(n_state,1);   % Initial state values
%
% options = odeset('RelTol',1e-02);   % the real error tolerance is 1e-02
%
% %   The numerical solutions are viewed as synthetic data
% sol = ode45(@diffun,[0:dt:T],x_true(:,1)',options);
%
% x_true = deval(sol,0:dt:T);
% P_obs_noise = 0 * diag(var(x_true'));
%
% z_true = x_true + sqrtm(P_obs_noise) * randn(n_state,N);
% dxdt = (diff(x_true'))'/dt;
% a = [0.1 0 0 0 0.05 0.025 0; 0.1 0.1 0 0 0.1 0 0;0.1 0 0 0.1 0 0 0 ];
% b = [0.1 0 10 0.05 0.025 0;0.1 0 0 0.1 10 0;0 0.1 0.1 0 0 0];
% ab_true = [a(1,:) b(1,:) a(2,:) b(2,:) a(3,:) b(3,:)]';
% c = [0.1 0.1 0.1]';
% d = [1 2 1]';
% e = [0.5 0.5 0.5 ]';
% a = [0.1 0 0 0 0.05 0.025 0;
%      0.1 0.1 0 0 0.1 0 0;
%      0.1 0 0.1 0 0 0 0];
% b = [1 0.1 0 10 0.05 0.025 0;
%      1 0.1 0 0 0.1 10 0;
%      1 0 0.1 0.1 0 0 0];
% c = [1 2 1]';
% d = 0.1*ones(3,1);
% e = 0.5*ones(3,1);
B = eye(n_state); %!!!!!!!!!!!!!!!!
if size(B,1) < n
    B(size(B,1),n)=0; % !!!!!!!!!!!!!!!!linear measurement, h(x) = B*x
end
% lamdda = zeros(n,1);% sparsity indicator
% lamdda(n_state + n_pro*n_cpx+ 1:n_state + 2*n_pro*n_cpx) = 0.01*ones(n_pro*n_cpx,1);
MC_N = 1;% Monte Carlo steps
% a_hat_mc = zeros(numel(a),MC_N);% estimation of a
% b_hat_mc = zeros(numel(b),MC_N);% estimation of b
mse_s = zeros(N ,MC_N);
mse_us = zeros(N,MC_N);
err_ab = zeros(N,1);
err_abu = zeros(N,1);
%ftrue2 = @(x)ode_ftrue2(x,dt,n_gene,n_pro,a,b,c,d,e);
% Monte Carlo run
for mc = 1: MC_N
    fprintf('This is %d monte carlo\n',mc);
    tic
    % initialize coefficients and states
    a_hat = randn(n_gene,n_cpx)
    b_hat = randn(n_gene,n_cpx - 1)
    s_t = 5*abs(randn(n_state,1));
    %z = s_t + r*randn(n_state,1); % observation of the x0
    %s_t = ftrue2(s_t) + r*randn(n_state,1); % x1
    x_hat = s_t;
    xu_hat = x_hat;
    P_hat = eye(n_state,n_state);%,eye(n-n_state));   % initial covariance
    Pu_hat = P_hat;
    % generate independent noise, which may not be necessary when the time
    % intervel is not long enough
    %     Q_Noise = mvnrnd(zeros(1,n_state),Q_x,N+1);% generate noise matrix
    %     Q_Noise = Q_Noise';
    %     R_Noise = mvnrnd(zeros(1,n_state),R,N+1);
    %     R_Noise = R_Noise';
    %[tt,xx] = ode45(@diffun,[0:dt:dt*N],ones(4,1),[]);
    %zV = xx' + R_Noise;
    xV = zeros(n,N); %estmated states
    xuV = xV;
    sV = zeros(n_state,N); %actual states
    aV = zeros(N,numel(a));
    bV = zeros(N,numel(b));
    thea1 = [a_hat(1,:) b_hat(1,:)]';
    P_para1 = 100*eye(numel(thea1));
    thea2 = [a_hat(2,:) b_hat(2,:);]';
    P_para2 = 100*eye(numel(thea2));
    %    thea3 = [a_hat(3,:) b_hat(3,:);]';
    %    P_para3 = 100*eye(numel(thea3));% 3 gene
    A1 = [];
    A2 = [];
    A3 = [];
    % this is an alternating way to update the state and parameter, that
    % means we fixed one and update the other, where the covariance is not
    % important?
    for k=1:N
        h = @(x)ode_h(x,n_state,n);%!!!!!!!!!!!!!!
        f_2 = @(x)fixed_f(x,dt,n_gene,n_pro,a_hat,b_hat,c,d,e);
        ftrue2 = @(x)ode_ftrue2(x,dt,n_gene,n_pro,a,b,c,d,e);
        z = h(s_t) + r*randn(size(R,1),1); %!!!!!! when observation is partial measurments, two places
        %z = z_true(:,k);
        [x_hat,P_hat]=ukf(f_2,x_hat,P_hat,h,z,Q,R);
        [xu_hat,Pu_hat]=ukf(f_2,xu_hat,Pu_hat,h,z,Q,R);
        xV(:,k) = x_hat;
        xuV(:,k) = xu_hat;
        if k >= 2
            % recursive least square
            x_pre = xV(:,k-1);
            x_temp = xV(:,k);
            cpx = pro_cpx(x_pre(n_gene + 1:n_state));
            y_pre = x_pre(n_gene + 1:n_state);
            y_temp = x_temp(n_gene + 1:n_state);
            y_temp - y_pre -d'*x_pre(1:n_gene) + e'*y_pre
            % to test the amplitude of the noise, which may be color since
            % there is a term in the denominator
            %         w1 = 1 + b_hat(1,:)*cpx(2:end)
            %         w2 = 1 + b_hat(2,:)*cpx(2:end)
            A1(:,k) = [-cpx*dt; (x_temp(1) - (1 - c(1)*dt)*x_pre(1))*cpx(2:end)];
            A2(:,k) = [-cpx*dt; (x_temp(2) - (1 - c(2)*dt)*x_pre(2))*cpx(2:end)];
            %A3(:,k) = [-cpx*dt; (x_temp(3) - (1 - c(3)*dt)*x_pre(3))*cpx(2:end)];%!!!!!!!!!
            b1(k) = (1 - c(1)*dt)*x_pre(1) - x_temp(1);
            b2(k) = (1 - c(2)*dt)*x_pre(2) - x_temp(2);
            %b3(k) = (1 - c(3)*dt)*x_pre(3) - x_temp(3);
            [thea1 P_para1] = myRLS(A1(:,k-1),b1(k-1),thea1,P_para1,0.1,0,0.1,k);%!!!!
            [thea2 P_para2] = myRLS(A2(:,k-1),b2(k-1),thea2,P_para2,0.1,0,0.1,k); %!!!!
            %[thea3 P_para3] = myRLS(A3(:,k-1),b3(k-1),thea3,P_para3,0.9,0,0.1,k);
            
            n1 = size(a,2);
            n2 = size(b,2);
            %         if k >= 10
            %         cvx_begin quiet
            %             variables aa(n1) bb(n2)
            %             minimize (([aa; bb] - thea1)'*inv(P_para1)*([aa; bb] - thea1) + 1e-6*norm(bb,1))
            %             aa >= 0;
            %             aa <= [1 ;bb]
            %         cvx_end
            %         thea1 = [aa;bb];
            %
            %         cvx_begin quiet
            %             variables aa2(n1) bb2(n2)
            %             minimize (([aa2; bb2] - thea2)'*inv(P_para2)*([aa2; bb2] - thea2) + 0*norm(bb2,1))
            %             subject to
            %             aa2 >= 0;
            %             aa2 <= [1 ;bb2]
            %         cvx_end
            %         thea2 = [aa2;bb2];
            %         end
            %
            %         cvx_begin quiet
            %             variables aa3(n1) bb3(n2)
            %             minimize (([aa3; bb3] - thea3)'*inv(P_para3)*([aa3; bb3] - thea3) + 1e-5*norm(bb3,1))
            %             subject to
            %             aa3 >= 0;
            %             aa3 <= [1 ;bb3]
            %         cvx_end
            %         thea3 = [aa3;bb3];
            %         end
            a_hat = [thea1(1:n_cpx) thea2(1:n_cpx) ]';%thea3(1:n_cpx)]';
            b_hat = [thea1(n_cpx + 1:end) thea2(n_cpx + 1:end)]';% thea3(n_cpx + 1:end)]';
        end
        %[x_hat, x_kal, P_hat] = ukf_sparse(f_2,x_hat,n_gene,n_protein,P_hat,h,B,lamdda,z,Q,R);
        sV(:,k)= s_t; % save actual state
        mse_s(k,mc) = norm(x_hat(1:n_state) - s_t,'fro');
        mse_us(k,mc) = norm(xu_hat(1:n_state) - s_t,'fro');
        %     aV(k,:) = x_hat(n_state + 1 : n_state + n_pro*n_cpx);
        %     bV(k,:) = x_hat(n_state + n_pro*n_cpx+  1 : end);
        %     auV(k,:) = xu_hat(n_state + 1 : n_state + n_pro*n_cpx);
        %     buV(k,:) = xu_hat(n_state + n_pro*n_cpx+  1 : end);
        %     a_hat = vec2mat(aV(k,:),n_cpx);
        %     b_hat = vec2mat(bV(k,:), n_cpx - 1);
        %     au_hat = vec2mat(auV(k,:),n_cpx);
        %     bu_hat = vec2mat(buV(k,:), n_cpx - 1);
        
        s_t = ftrue2(s_t) + q*randn(n_state,1);
        ab_hat = [a_hat b_hat]';
        %     abu_hat = [au_hat bu_hat]';
        ab_hat = ab_hat(:);
        %     abu_hat = abu_hat(:);
        err_ab(k) = norm(ab_hat - ab_true,'fro');
        %     err_abu(k) = norm(abu_hat - ab_true,'fro');
    end
    
    %     a_hat_mc(:,mc) = x_hat(n_state + 1 : n_state + n_protein*n_cpx);
    %     b_hat_mc(:,mc) = x_hat(n_state + n_protein*n_cpx +  1 : end);
    
    toc
    ab_mse = [inv(A1*A1')*A1*b1';inv(A2*A2')*A2*b2'];%inv(A3*A3')*A3*b3'];
    
end
%     mean_a_hat = mean(a_hat_mc');% mean of the estimated a
%     mean_b_hat = mean(b_hat_mc');% mean of the estimated b
%     std_a_hat = std(a_hat_mc'); % standard deviation of the estimated a
%     std_b_hat = std(b_hat_mc'); % standard deviation of the estimated b
%     a_hat_fl = vec2mat(mean_a_hat,n_cpx); % final average a_hat
%     b_hat_fl = vec2mat(mean_b_hat, n_cpx); % final average b_hat
cvx_begin quiet
variables aa1(n1) bb1(n2)
minimize (sum_square([aa1;bb1]'*A1-b1) + 0.01*norm(bb1))
subject to
aa1 >= 0
aa1 <= [1;bb1]
cvx_end

cvx_begin quiet
variables aa2(n1) bb2(n2)
minimize (sum_square([aa2;bb2]'*A2-b2) + 0.01*norm(bb2))
subject to
aa2 >= 0
aa2 <= [1;bb2]
cvx_end

S1 = svd(A1)
S2 = svd(A2)
tt = 0:1:N-1;
subplot(3,1,1)
plot(tt,sV(1,:),'r-',tt,sV(2,:),'r-.',tt,sV(3,:),'r.',tt,sV(4,:),'r*',tt,(xV(1,:))','b-',tt,(xV(2,:))','b-.',tt,(xV(3,:))','b.',tt,(xV(4,:))','b*')
legend('true x1','true x2','true y1', 'true y2','AU x1','AU x2','AU y1','AU y2')
xlabel('steps')
ylabel('states')
subplot(3,1,2)
plot(tt,sV(1,:),'r-',tt,sV(2,:),'r-.',tt,sV(3,:),'r.',tt,sV(4,:),'r*',tt,(xuV(1,:))','b-',tt,(xuV(2,:))','b-.',tt,(xuV(3,:))','b.',tt,(xuV(4,:))','b*')
legend('true x1','true x2','true y1', 'true y2','U x1','U x2','U y1','U y2')
xlabel('steps')
ylabel('states')
subplot(3,1,3)
plot(tt,err_ab,tt,err_abu)
legend('alternative ukf','ukf')
xlabel('steps')
ylabel('error of parameter')