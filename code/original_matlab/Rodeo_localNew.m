function [Hfinal, Hseq, p, diff2] = Rodeo_localNew(x, xtest, beta, c0, df2)
%
% Rodeo Probabilistic Classification (training)
%
%
% Input:
%    x : N by d training sample matrix
%    xtest : n ny d matrix that want to be estimated 
%    beta :  shrink parameter
%    c0 : constant
%
% Output:
%    Hseq : bandwidth sequences
%    Hfinal : final selected bandwidths n*d
%    p : estimated probability

    [ntest, d]=size(xtest);

    N = size(x, 1);
    Hseq = cell(ntest, 1);
    Hfinal = zeros(ntest, d);
 
    p = zeros(ntest, 1);
    
    diff2 = zeros(ntest, d);
    for m = 1 : ntest
        xx = xtest(m, :);

        Hseq{m, 1} = repmat(c0/(log(log(N))), [1, d]); % 1*d

        A = 1 : d;
        t = 1;

        while ~isempty(A)
            t = t + 1;
            % Uniform
            
            Z = D1(x, xx, Hseq{m, 1}(t-1,:));
            S = var(Z) / N; % 1*d
            T = mean(Z); % test statistic 1*d
                            
            for j = 1 : d
                if ismember(j, A)
                    lambda = sqrt(S(j)) * sqrt(2 * log(N*log(d)));
                    if abs(T(j)) > lambda
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
        p(m) = mean(GK(x, xx, Hfinal(m,:), 1));

        if df2 == 1
            diff2(m,:) = mean(D2(x, xx, Hfinal(m,:)));
        end
        
    end
end

function result = D1(X, mu, H)
    % X : N*d matrix
    % H : 1*d bandwidths
    % result : N*d

    [~, d] = size(X);
    tmp1 = ((X - mu).^2 - H.^2) ./ H.^3;  % N*d
    tmp2 = GK(X, mu, H, 1); % N*1  mvnpdf(X, mu, diag(H.^2));
    result = (2*pi)^(d/2) * tmp1 .* tmp2; % N*d 
end




function result = D2(X, mu, H)
    % X : N*d matrix
    % H : 1*d bandwidths

    tmp1 = -3 * (X - mu).^2 .* H.^(-4) + H.^(-2); % N*d
    tmp2 = ((X - mu).^2 - H.^2) ./ H.^3;  % N*d
    tmp3 = GK(X, mu, H, 1); % N*1
    result = tmp1 .* tmp3 + tmp2 .* tmp2 .* tmp3; % N*d
end


function result = GK(X, mu, H, D)

    % X : N*p matrix
    % H : 1*p bandwidths
    % D : divededd by prod of bandwidths or not; 
    
    % result : N*1 matrix
    %     N = size(X, 1);

    if D == 0
        result = mvnpdf((X-mu) ./H);
    else
        result = mvnpdf(X, mu, diag(H.^2));
    end    
end