function [w_x, w_y,r] = compute_weights_cca(C_xx,C_yy,C_xy,f, gamma)
%@Cxx covariance matrix of data at time X
%@Cyy covariance matrix of data at time Y
%@Cxy covariance matrix between data at time X and y

% computes the weights w_x and w_y for Dimension f such that w_x*Data_x is maximally
% correlated with Data_y * w_y

if nargin>4
    I = eye(size(C_xx,1));
    C_xx = (1-gamma)* C_xx + gamma*I*trace(C_xx);
    I = eye(size(C_yy,1));
    C_yy = (1-gamma)* C_yy + gamma*I*trace(C_yy);
end

[w_x,D] = eigs(C_xy*inv(C_yy)* (C_xy)', C_xx, f);
[w_y,~] = eigs(C_xy'*inv(C_xx)* C_xy, C_yy, f);
if nargout>2
    r = sqrt(diag(D)');
end

% new fix for consistency
for idx = 1: f
    if w_x(:,idx)'*C_xy*w_y(:,idx) < 0
        w_x(:,idx) = -w_x(:,idx); % The sign of "lbd" is not constrained. Make sure we maximize correlation and avoid negative correlations
    end
end


end

