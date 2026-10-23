function [Hfinal, Hseq, value, p, diff2] = RodeoPCNew(x, y ,xtest, beta, c0, df2, uni, nor)
%
% Rodeo Probabilistic Classification (training)
%
%
% Input:
%    x : N by d training sample matrix
%    y : N by 1 training label vector (taking 1,...,c)
%    xtest : n by d testing sample matrix
%    beta : shrink parameter
%    c0 : constant
%
% Output:
%    Hseq : bandwidth sequences -- cell(# of groups)
%    Hfinal : final selected bandwidths n * d * (# of groups)
%    value : estimated gourp
%    p : estimated probability -- n * (# of groups)
%

    [n, d]=size(xtest);
    [C,~,y_new] = unique(y); % Number of classes
    c = max(C); % Number of classes
    X = cell(c, 1);
    idx = cell(c,1);
    N = zeros(1, c); % # of samples in each group

    if uni +nor == 2
        Hfinal = cell(2, 1);
        Hseq = cell(c, 2);
        p = cell(2, 1);
        value = zeros(n, 2);
        diff2 = cell(c, 1);

        for yy = 1 : c
            tmp = y_new == yy;

            X{yy} = x(tmp, :);
    %         Y{yy} = y(tmp, :);
            idx{yy} = tmp;
            N(yy) = size(tmp, 1);

            [Hfinal{1,1}(:, :, yy), Hseq{yy,1}, p{1,1}(:, yy), diff2{yy,1}] = Rodeo_localNew(X{yy}, xtest, beta, c0, df2);
            [Hfinal{2,1}(:, :, yy), Hseq{yy,2}, p{2,1}(:, yy)] = Rodeo_localNormal(X{yy}, xtest, beta, c0);     
        end
        [~, value(:,1)] = max(p{1,1} ./ repmat(sum(p{1,1}, 2), 1, c),[],2); % value: n*1
        [~, value(:,2)] = max(p{2,1} ./ repmat(sum(p{2,1}, 2), 1, c),[],2); % value: n*1
        
    elseif uni + nor == 1
        
        Hfinal = zeros(n, d, c);
        Hseq = cell(c, 1);

        p = zeros(n, c);
        diff2 = cell(c, 1);
        
        for yy = 1 : c
            tmp = y_new == yy;

            X{yy} = x(tmp, :);
    %         Y{yy} = y(tmp, :);
            idx{yy} = tmp;
            N(yy) = size(tmp, 1);
        
            if uni == 1
                [Hfinal(:, :, yy), Hseq{yy, 1}, p(:, yy), diff2{yy,1}] = Rodeo_localNew(X{yy}, xtest, beta, c0, df2);
            else
                [Hfinal(:, :, yy), Hseq{yy, 1}, p(:, yy)] = Rodeo_localNormal(X{yy}, xtest, beta, c0);    
            end
        end
        [~, value] = max(p ./ repmat(sum(p, 2), 1, c),[],2); % value: n*1


    end

        
end



