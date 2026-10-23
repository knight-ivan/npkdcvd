function [Acc, AccStd, Time, TimeStd, Z1, HH, AA, truev] = Ex1Com(nc, t, niter, d, nvar, ntrain, ntest, uni, nor)

    % Input:
    % nc : # of groups in the data
    % t : # of groups used in comparision
    % niter # of iterations;

    % d : dim of predictors
    % nvar : # of related predictors
    
    % ntrain : # of training samples
    % ntest : # of testing samples 
    
    % uni :  
    % nor 
    
    % Acc & AccStd : mean and std of accuracy for each group
    % Time & TimeStd : mean and std of time for learning plus testing

    % Z1 : mean of h for each class and variables
    %      Z1{1:c,1} : t * d

    % AA : the bandwidth mean of the testing samples predicted in each
    %      group 1:t for each iteration
    %      AA{k,1} : t * niter * d

    % HH : z-score of bandwidth mean of the testing samples predicted in each
    %      group 1:t for each iteration
    %      HH{k, 1}(iter,:,i) = zscore(AA{k, 1}(i, iter, :))
    %      HH{k,1} : niter * d * t   
    
    
    v = zeros(nc, nvar);
    nv = zeros(nc, d-nvar);
    data = cell(nc, 1);

    truev = zeros(nc, d);

    for l = 1 : nc
        data{l,1} = zeros(1000, d);
        v(l,:) = l : l+nvar-1; % related vars in each group
        nv(l,:) = setdiff(1:d, v(l,:)); % non-related vars in each group
        for k = 1:nvar
            data{l,1}(:, v(l,k)) = normrnd(0.5, (0.02*k), [1000, 1]);
        end
        for k = 1:(d-nvar)
            data{l,1}(:, nv(l,k)) = rand(1000, 1);
%             data{l,1}(:, nv(l,k)) = normrnd(0, 1, [1000, 1]);
        end
        truev(l, l:l+nvar-1) = 1;
    end

    c = combnk(1:nc,t); % # of combination

    Acc = zeros(size(c,1), 2);
    AccStd = zeros(size(c,1), 2);
    Time = zeros(size(c,1), 2);
    TimeStd = zeros(size(c,1), 2);

    Z1 = cell(size(c,1), 1);
    HH = cell(size(c,1), 1);
    AA = cell(size(c,1), 1);
    parfor k = 1:size(c,1)
        tmp = data;
        cc = c(k,:);  % which groups in kth combination
        X = cell2mat(tmp(cc,1));
        Y = ones(1000*t,1);
        for i = 1 : t
            Y(1000*(i-1)+1:1000*i) = i;
        end

        AAcc = zeros(niter, 2);
        HH1 = zeros(niter, d, t); % �C��iteration ���O�bt����z score
        Ttmp = zeros(niter, 2);
 
        AA1 = zeros(t, niter, d);
        for iter = 1 : niter
            Xtrain = zeros(t*ntrain, d);
            Ytrain = zeros(t*ntrain, 1);
            Xtest = zeros(t*ntest, d);
            Ytest = zeros(t*ntest, 1);

            id = cell(t,1);
            id_train = cell(t,1);
            id_test = cell(t,1);    

            for i = 1 : t

                id{i,1} = find(Y == i);
                id_train{i, 1} = id{i,1}(randsample(length(id{i,1}), ntrain));
                Xtrain((i-1)*ntrain+1:i*ntrain, :) = X(id_train{i, 1}, :);
                Ytrain((i-1)*ntrain+1:i*ntrain, 1) = Y(id_train{i, 1});

                diff = setdiff(id{i,1}, id_train{i, 1});
                id_test{i,1} = randsample(diff, ntest);
                Xtest((i-1)*ntest+1:i*ntest, :) = X(id_test{i, 1}, :);
                Ytest((i-1)*ntest+1:i*ntest, 1) = Y(id_test{i, 1});
            end 

            value = zeros(t*ntest, 2); % estimated group

            tic
            [tmp1, ~, value(:,1), ~, ~] = RodeoPCNew(Xtrain, Ytrain ,Xtest, 0.9, 15, 0, uni, nor);           % tmp1 -- ntest * d * t
            
            Ttmp(iter, 1) = toc;

            tic
            LSPC_model=LSPC_train(Xtrain',Ytrain');
            [ytest_LSPC,~] = LSPC_test(Xtest',LSPC_model);
            value(:,2) = ytest_LSPC';
            Ttmp(iter, 2) = toc;

            
            for i = 1:t
                AA1(i, iter, :) = mean(tmp1(value(:,1) == i,:,i)); % the bandwidth mean of the testing samples predicted as group i
                HH1(iter,:,i) = zscore(AA1(i, iter, :));    % 1*d
            end

            
            AAcc(iter,:) = sum(value == repmat(Ytest, [1, 2]))/length(Ytest);

            disp([k, iter]);
        end

        Acc(k,:) = mean(AAcc);
        AccStd(k,:) = std(AAcc);
        Time(k,:) = mean(Ttmp);
        TimeStd(k,:) = std(Ttmp);
        for i = 1:t
            Z1{k,1}(i,:) = mean(HH1(:,:,i)); % ������Z  score 

        end
        HH{k,1} = HH1;
        AA{k,1} = AA1;
    end




