function [Hfinal, Hseq, p] = Rodeo_localNormal(x, xtest, beta, c0)
%
% Rodeo Probabilistic Classification (training)
%
%
% Input:
%    x : N by d training sample matrix
%    xtest : ntest by d matrix that want to be estimated 
%    beta :  shrink parameter
%    c0 : constant
%
% Output:
%    Hseq : bandwidth sequences
%    Hfinal : final selected bandwidths
%    p : estimated probability

    [ntest, d]=size(xtest);

    N = size(x, 1);
    Hseq = cell(ntest, 1);
    Hfinal = zeros(ntest, d);
 
    p = zeros(ntest, 1);
    
    mu = mean(x); % 1*d
    sigma = std(x);
    
    for m = 1 : ntest
        xx = xtest(m, :);

        Hseq{m, 1} = repmat(c0/(log(log(N))), [1, d]); % 1*d

        A = 1 : d;
        t = 1;

        while ~isempty(A)
            t = t + 1;
            disp(t);
            Z = Rodeo_GradNormal(x, xx, mu, sigma, Hseq{m, 1}(t-1,:)); % 1*d
            
            z = zeros(20, d);
            for k = 1:20
                id = randsample(N, N, true);    % sampling with replacement
                z(k, :) = Rodeo_GradNormal(x(id,:), xx, mu, sigma, Hseq{m, 1}(t-1,:)); % 1*d
            end
            s2 = mean((z - mean(z)).^2); % 1*d estimated variance
            

           for j = 1 : d
                if ismember(j, A)
                    lambda = sqrt(s2(j)) * sqrt(2 * log(N*log(d))); %*
                    if abs(Z(j)) > lambda*10
                        Hseq{m, 1}(t, j) = beta * Hseq{m, 1}(t-1, j);
                    else
                        Hseq{m, 1}(t, j) = Hseq{m, 1}(t-1, j);
                        A = setdiff(A, j);
                    end

                else
                    Hseq{m, 1}(t, j) = Hseq{m, 1}(t-1, j);
                end
            end


        end

        Hfinal(m,:) = Hseq{m, 1}(end, :);
        
        tmp1 = mean(exp(-0.5 * sum((x - xx).^2 ./ repmat(Hfinal(m,:).^2, N, 1).^2, 2)));
        tmp2 = sqrt(prod(Hfinal(m,:).^2 + sigma.^2)/prod(sigma.^2));
        tmp3 = exp(-0.5 * (xx - mu) * (inv(diag(sigma.^2)) - inv(diag(Hfinal(m,:).^2 + sigma.^2))) * (xx - mu)');
        p(m) = tmp1 * tmp2 * tmp3;
    end
end


