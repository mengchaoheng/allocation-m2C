function [u] = dyn_alloc_m1(v,u,only_plim)
  
% WLS_ALLOC - Control allocation using weighted least squares.
%
%  [u,W,iter] = wls_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,W0,imax])
%
% Solves the weighted, bounded least-squares problem
%
%   min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
%
%   subj. to  umin <= u <= umax
%
% using an active set method.
%
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
% Wv    virtual control weighting matrix (k x k) [I]
% Wu    control weighting matrix (m x m) [I]
% ud    desired control (m x 1) [0]
% gamma weight (scalar) [1e6]
% u0    initial point (m x 1)
% W0    initial working set (m x 1) [empty]
% imax  max no. of iterations [100]
% 
%  Outputs:
%  -------
% u     optimal control
% W     optimal active set
% iter  no. of iterations (= no. of changes in the working set + 1)
%
%                            0 if u_i not saturated
% Working set syntax: W_i = -1 if u_i = umin_i
%                           +1 if u_i = umax_i
%
% See also: WLSC_ALLOC, IP_ALLOC, FXP_ALLOC, QP_SIM.
B=[-4.6718         0    4.6718         0   -4.6718         0    4.6718         0    3.2055;
         0  -14.9632         0   14.9632         0  -14.9632         0   14.9632         0;
    4.2635   -6.9597    4.2635   15.4868    4.2635   15.4868    4.2635   -6.9597         0];
% Number of variables
m = 9;
  if (only_plim)
    umin = [[1;1;1;1;1;1;1;1]*(-20*pi/180);-4.89];
    umax = [[1;1;1;1;1;1;1;1]*20*pi/180;4.89];
  else
    umin = max([[1;1;1;1;1;1;1;1]*(-20*pi/180);-4.89],u+[[1;1;1;1;1;1;1;1]*(-400*pi/180);-19.56]*0.01);
    umax = min( [[1;1;1;1;1;1;1;1]*20*pi/180;4.89],u+[[1;1;1;1;1;1;1;1]*400*pi/180;19.56]*0.01);
  end
  %========计算当前有效集============
W=zeros(m,1);
infeasible1 = (u <= umin);
W(infeasible1)=-1;
infeasible2 = (u >= umax);
W(infeasible2)= 1;
%===============期望舵机位置========================

W1sq = eye(9);
W2sq = [
    0.0100         0         0         0         0         0         0         0         0;
         0    0.0100         0         0         0         0         0         0         0;
         0         0    0.0100         0         0         0         0         0         0;
         0         0         0    0.0100         0         0         0         0         0;
         0         0         0         0    0.0100         0         0         0         0;
         0         0         0         0         0    0.0100         0         0         0;
         0         0         0         0         0         0    0.0100         0         0;
         0         0         0         0         0         0         0    0.0100         0;
         0         0         0         0         0         0         0         0    0.2500
    ];
S=[
    0         0         0;
         0   -0.0167   -0.0121;
         0         0         0;
         0    0.0167    0.0269;
         0         0         0;
         0   -0.0167    0.0269;
         0         0         0;
         0    0.0167   -0.0121;
    0.3120         0         0
];
us = S*v;
ud = (W1sq+W2sq)\(W1sq*us+W2sq*u);
%================================
% 加权系数
gam=1e6;
% 加权矩阵
Wv=eye(3);
Wu = sqrtm(W1sq+W2sq);
% 迭代次数上限
imax=100;
  % Set default values of optional arguments  
  gam_sq = sqrt(gam);
  A = [gam_sq*Wv*B ; Wu];
  b = [gam_sq*Wv*v ; Wu*ud];
  
  % Initial residual.
  d = b - A*u;
  % Determine indeces of free variables.
  i_free = W==0;
  
  % Iterate until optimum is found or maximum number of iterations
  % is reached.
  for iter = 1:imax
    % ----------------------------------------
    %  Compute optimal perturbation vector p.
    % ----------------------------------------
    
    % Eliminate saturated variables.
    A_free = A(:,i_free);
    % Solve the reduced optimization problem for free variables.
    p_free = A_free\d;
    % Zero all perturbations corresponding to active constraints.
    p = zeros(m,1);
    % Insert perturbations from p_free into free the variables.
    p(i_free) = p_free;
    
    % ----------------------------
    %  Is the new point feasible?
    % ----------------------------
    
    u_opt = u + p;
    infeasible = (u_opt < umin) | (u_opt > umax);

    if ~any(infeasible(i_free))

      % ----------------------------
      %  Yes, check for optimality.
      % ----------------------------
      
      % Update point and residual.
      u = u_opt;
      d = d - A_free*p_free;
      % Compute Lagrangian multipliers.
      lambda = W.*(A'*d);
      % Are all lambda non-negative?
      if lambda >= -eps
	% / ------------------------ \
	% | Optimum found, bail out. |
	% \ ------------------------ /
	return;
      end
      
      % --------------------------------------------------
      %  Optimum not found, remove one active constraint.
      % --------------------------------------------------
      
      % Remove constraint with most negative lambda from the
      % working set.
      [lambda_neg,i_neg] = min(lambda);
      W(i_neg) = 0;
      i_free(i_neg) = 1;
    
    else
      
      % ---------------------------------------
      %  No, find primary bounding constraint.
      % ---------------------------------------
      
      % Compute distances to the different boundaries. Since alpha < 1
      % is the maximum step length, initiate with ones.
      dist = ones(m,1);
      i_min = i_free & p<0;
      i_max = i_free & p>0;

      dist(i_min) = (umin(i_min) - u(i_min)) ./ p(i_min);
      dist(i_max) = (umax(i_max) - u(i_max)) ./ p(i_max);

      % Proportion of p to travel
      [alpha,i_alpha] = min(dist);
      % Update point and residual.
      u = u + alpha*p;
      d = d - A_free*alpha*p_free;
      
      % Add corresponding constraint to working set.
      W(i_alpha) = sign(p(i_alpha));
      i_free(i_alpha) = 0;
      
    end
  
  end