


%% Implementation of Hard margin and Soft margin Regularized Support Vector Machine Algorithm

%% Generate Plots for Problem 1 and compute misclassification error
clear all; clc;
% 3 kinds of data


% case 1) Linearly Separable Data

load('linear_new.mat');
[b, b0] = HardMarg(X,y);
[bs, b0s] = SoftMarg(X,y,0.5);
figure(1);
[plt_svm] = plotSVM_Hplane(b,b0,bs,b0s,X,y,'Hard Margin','Soft Margin');
xlabel('X1');
ylabel('X2');
title('Linear Data - SVM Separating Hyperplanes');
set(gca,'FontSize',16);
mclerr_test_hard  = compute_error(Xtest,ytest,b,b0);
mclerr_train_hard = compute_error(X,y,b,b0);
mclerr_test_soft  = compute_error(Xtest,ytest,bs,b0s);
mclerr_train_soft = compute_error(X,y,bs,b0s);
formatSpec = 'Linear Data Training error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_hard);
formatSpec = 'Linear Data Training error for Soft Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_soft);
formatSpec = 'Linear Data Test error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_test_hard);
formatSpec = 'Linear Data Test error for Soft Margin = %4.2f \n\n' ;
fprintf(formatSpec,mclerr_test_soft);
              

% case 2) SVM Quadratically Separable Data
clearvars b b0 bs b0s;
load('quadratic_new.mat');
[b, b0] = HardMarg(X,y);
[bs, b0s] = SoftMarg(X,y,0.5);
figure(2);
[plt_svm] = plotSVM_Hplane(b,b0,bs,b0s,X,y,'Hard Margin','Soft Margin');
xlabel('X1');
ylabel('X2');
title('Quadratic Data - SVM Separating Hyperplanes');
set(gca,'FontSize',16);
mclerr_test_hard  = compute_error(Xtest,ytest,b,b0);
mclerr_train_hard = compute_error(X,y,b,b0);
mclerr_test_soft  = compute_error(Xtest,ytest,bs,b0s);
mclerr_train_soft = compute_error(X,y,bs,b0s);
formatSpec = 'Quad Data Training error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_hard);
formatSpec = 'Quad Data Training error for Soft Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_soft);
formatSpec = 'Quad Data Test error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_test_hard);
formatSpec = 'Quad Data Test error for Soft Margin = %4.2f \n\n' ;
fprintf(formatSpec,mclerr_test_soft);


% case 3) SVM for Noisy Data
clearvars b b0 bs b0s;
load('noisy-linear_new.mat');
[b, b0] = HardMarg(X,y);
[bs, b0s] = SoftMarg(X,y,0.5);
figure(3);
[plt_svm] = plotSVM_Hplane(b,b0,bs,b0s,X,y,'Hard Margin','Soft Margin');
xlabel('X1');
ylabel('X2');
title('Noisy Linear Data - SVM Separating Hyperplanes');
set(gca,'FontSize',16);
mclerr_test_hard  = compute_error(Xtest,ytest,b,b0);
mclerr_train_hard = compute_error(X,y,b,b0);
mclerr_test_soft  = compute_error(Xtest,ytest,bs,b0s);
mclerr_train_soft = compute_error(X,y,bs,b0s);
formatSpec = 'Noisy Linear Data Training error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_hard);
formatSpec = 'Noisy Linear Data Training error for Soft Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_train_soft);
formatSpec = 'Noisy Linear Data Test error for Hard Margin = %4.2f \n' ;
fprintf(formatSpec,mclerr_test_hard);
formatSpec = 'Noisy Linear Data Test error for Soft Margin = %4.2f \n\n' ;
fprintf(formatSpec,mclerr_test_soft);