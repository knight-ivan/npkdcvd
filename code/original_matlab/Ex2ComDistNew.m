function [Acc, AccStd, Time, TimeStd, ZZ0, ZZ1, Hmean00, Hmean11] = Ex2ComDistNew(t, dist, niter, ntrain, ntest, uni, nor)

    % Input:
    % dist : dist between groups
    
    d = 10;
    nc = 5;

    
    dir = [0.5, 0.7, -0.5];
    dir = dir/sqrt(sumsqr(dir));

    mu0 = [1, -3, 2];
    mu = zeros(5, 3);
    mu(1, :) = mu0;
    for l = 2 : nc
        mu(l, :) = mu0 +(l-1)*dist*dir;
    end

    data = cell(5, 1);
    for l = 1 : nc
        data{l,1} = zeros(1000, d);

        for k = 1:3
            data{l,1}(:, k) = normrnd(mu(l,k), (0.2), [1000, 1]);
        end
        for k = 4:d
%             data{l,1}(:, k) = normrnd(0, 1, [1000, 1]);
            data{l,1}(:, k) = rand(1000, 1);
        end
    end
    
%     ntrain = 200;
%     ntest = 150;
% 
    c = combnk(1:nc,t); % # of combination

    Acc = zeros(size(c,1), 2);
    AccStd = zeros(size(c,1), 2);
    Time = zeros(size(c,1), 2);
    TimeStd = zeros(size(c,1), 2);

    ZZ0 = cell(size(c,1), 1);
    ZZ1 = cell(size(c,1), 1);
    Hmean00 = cell(size(c,1), 1); 
    Hmean11 = cell(size(c,1), 1); 
    parfor k = 1:size(c,1)
%     for k = 1:2
        tmp = data;
        cc = c(k,:);
        X = cell2mat(tmp(cc,1));
        Y = ones(1000*t,1);
        for i = 1 : t
            Y(1000*(i-1)+1:1000*i) = i;
        end

        AAcc = zeros(niter, 2);
        Hmean0 = zeros(niter, d, t); % �C��iteration ���O�bt���� mean of bandwidths 
        Z0 = zeros(niter, d, t); % �C��iteration ���O�bt���� z score
        Hmean1 = zeros(niter, d, t); % �C��iteration ���O�bt���� mean of bandwidths 
        Z1 = zeros(niter, d, t); % �C��iteration ���O�bt���� z score
        
        Ttmp = zeros(niter, 2);
        for iter = 1 : niter
            Xtrain = zeros(t*ntrain, d);
            Ytrain = zeros(t*ntrain, 1);
            Xtest = zeros(t*ntest, d);
            Ytest = zeros(t*ntest, 1);

            id = cell(t,1);
            id_train = cell(t,1);
            id_test = cell(t,1);    

            Htrain = zeros(ntrain, d, t);
            for i = 1 : t

                id{i,1} = find(Y == i);
                id_train{i, 1} = id{i,1}(randsample(length(id{i,1}), ntrain));
                tmpX = X(id_train{i, 1}, :);
                Xtrain((i-1)*ntrain+1:i*ntrain, :) = tmpX;
                Ytrain((i-1)*ntrain+1:i*ntrain, 1) = Y(id_train{i, 1});

                diff = setdiff(id{i,1}, id_train{i, 1});
                id_test{i,1} = randsample(diff, ntest);
                Xtest((i-1)*ntest+1:i*ntest, :) = X(id_test{i, 1}, :);
                Ytest((i-1)*ntest+1:i*ntest, 1) = Y(id_test{i, 1});

                [Htrain(:, :, i), ~, ~, ~] = Rodeo_localNew(tmpX, tmpX, 0.9, 5, 0);
                Hmean0(iter, :, i) = mean(Htrain(:, :, i)); 
                % 1*d: the bandwidth mean of the training data in group i         
                Z0(iter,:,i) = zscore(Hmean0(iter, :, i)); % 1*d
            end 

            value = zeros(t*ntest, 2); % estimated group

            tic
            [tmp1, ~, value(:,1), ~] = RodeoPCNew(Xtrain, Ytrain ,Xtest, 0.9, 5, 0, uni, nor); % tmp1 ntest * d * t
            Ttmp(iter, 1) = toc;

            for i = 1:t
                Hmean1(iter, :, i) = mean(tmp1(value(:,1) == i,:,i)); 
                % 1*d: the bandwidth mean of the testing samples predicted as group i
                Z1(iter,:,i) = zscore(Hmean1(iter, :, i));    % 1*d
            end

        
            tic
            LSPC_model=LSPC_train(Xtrain',Ytrain');
            [ytest_LSPC,~] = LSPC_test(Xtest',LSPC_model);
            value(:,2) = ytest_LSPC';
            Ttmp(iter, 2) = toc;

            AAcc(iter,:) = sum(value == repmat(Ytest, [1, 2]))/length(Ytest);

            disp([k, iter]);
        end

        Acc(k,:) = mean(AAcc);
        AccStd(k,:) = std(AAcc);
        Time(k,:) = mean(Ttmp);
        TimeStd(k,:) = std(Ttmp);
        Hmean00{k,1} = Hmean0;
        Hmean11{k,2} = Hmean1;
        for i = 1:t
            ZZ0{k, 1}(i, :) = mean(Z0(:,:,i)); % ������Z  score 
            ZZ1{k, 1}(i, :) = mean(Z1(:,:,i)); % ������Z  score 
        end

    end




