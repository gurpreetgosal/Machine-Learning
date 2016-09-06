%% Generate Training and Test data
clear all;
% define parameters for the data
d = 10;
n_train = 2000;

% generate  d standard gaussians each with 3000 points
X = normrnd(0,1,n_train+1000,d);

% select 2000 samples such that 1000 samples have label 1 and another
% 1000 have label -1

% first find labels for each of the 3000 samples
y=zeros(length(X(:,1)),1);
for i=1:length(X(:,1))
    if sumsqr(X(i,:)) > 9.34
        y(i) = 1;
    elseif sumsqr(X(i,:)) <= 9.34
        y(i) = -1;
    end
end

% then select 1000 for each class
ind_class1 = find(y==1);
ind_class2 = find(y==-1);

X_class1 = X(ind_class1,:);
X_class2 = X(ind_class2,:);
y_class1 = y(ind_class1);
y_class2 = y(ind_class2);

Xtrain = [X_class1(1:1000,:); X_class2(1:1000,:)];
ytrain = [y_class1(1:1000); y_class2(1:1000)];


% Generate 10000 test data samples
n_test = 10000;

% generate  d standard gaussians each with 3000 points
Xtest_p = normrnd(0,1,n_test+5000,d);

% select 2000 samples such that 1000 samples have label 1 and another
% 1000 have label -1

% first find labels for each of the 3000 samples
ytest_p=zeros(length(Xtest_p(:,1)),1);
for i=1:length(Xtest_p(:,1))
    if sumsqr(Xtest_p(i,:)) > 9.34
        ytest_p(i) = 1;
    elseif sumsqr(Xtest_p(i,:)) <= 9.34
        ytest_p(i) = -1;
    end
end

% then select 1000 for each class
ind_tclass1 = find(ytest_p==1);
ind_tclass2 = find(ytest_p==-1);

X_tclass1 = Xtest_p(ind_tclass1,:);
X_tclass2 = Xtest_p(ind_tclass2,:);
y_tclass1 = ytest_p(ind_tclass1);
y_tclass2 = ytest_p(ind_tclass2);

Xtest = [X_tclass1(1:5000,:); X_tclass2(1:5000,:)];
ytest = [y_tclass1(1:5000); y_tclass2(1:5000)];

%% Start Adaboost Algorithm

% first set the weights
w_initial = (1/n_train)*ones(n_train,1);
w = w_initial;
% Let J be the number of classifiers
J = 500;
hxj_wsum = zeros(n_train,1);
for j = 1:J
    % next we find the classifier which minimizes metric Lj
    [L(j,1),s,hxj,sign_h] = decision_stump2(Xtrain, ytrain,n_train,d,w);
    
    alpha(j,1) = log((1-L(j,1))/L(j,1));
    %update weights
    for i = 1:n_train
        w(i,1) = w(i,1)*exp(alpha(j,1)*(hxj(i,1)~=ytrain(i)) );
    end
    % Normalize updated weights
    %w = w/sum(w);
    hxj_wsum = hxj_wsum + alpha(j,1)* hxj;
    Hx_afterj = sign(hxj_wsum);
    % compute training misclassification error rate
    error_train(j) = sum(Hx_afterj~=ytrain)/n_train;
    display(error_train(j))
    display(alpha(j,1),'alpha =');
    display(L(j,1),'L =');
    disp('---------------------------------------------------------')
    
    %compute test error
    %error_test = sum(s)
    
end









