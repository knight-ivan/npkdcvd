function Z = Rodeo_GradNormal(x, xtest, mu, sigma, H)
%
% Rodeo Probabilistic Classification (training)
%
%
% Input:
%    x : N by d training sample matrix
%    xtest : n ny d matrix that want to be estimated 
%    mu, sigma : parameters of normal distribution
%    H : 1*d vector of bandwidths 
%
% Output:
%    H : bandwidth sequences
%    HH : final selected bandwidths
%    p : estimated probability

    [ntest, d]=size(xtest);

    N = size(x, 1);    
    Z = zeros(ntest, d);
    

    for m = 1 : ntest
        xx = xtest(m, :);
        tmp1 = sqrt(prod((H.^2 + sigma.^2) ./ sigma.^2));    % 1*1
        tmp2 = -0.5 * sum((xx-mu) .^2 ./ (sigma.^2 .* (sigma .^2 + H.^2) ./ H.^2));  % 1*1
        tmp3 = mean(((x - xx).^2 - H.^2 ./ repmat(H.^3, N, 1)) .* repmat(exp(-0.5 * sum((x - xx).^2 ./ repmat(H.^2, N, 1), 2)), 1, d)) / prod(H); % 1*d
        M = (H .* (2*(sigma .^2 + H .^2) + (xx - mu).^2)) ./ (2 * (sigma .^2 + H .^2) .^2);  % 1*d
        tmp4 = tmp3 + M .* mvksdensity(x, xx, 'Bandwidth', H);   % 1*d
        Z(m, :) = tmp1 * exp(tmp2 * tmp4);    % 1*d
    end

end





