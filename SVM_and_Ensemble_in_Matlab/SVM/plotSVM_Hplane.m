function [plt_svm] = plotSVM_Hplane(b,b0,bs,b0s,X,y, str1,str2)

plot(X(1,y==1), X(2,y==1),'.','markersize',25);
hold on
plot(X(1,y==-1), X(2,y==-1),'.','markersize',25);
hold on
plt_svm = plot(X(1,:),(-1*b(1)/b(2))*X(1,:)-b0/b(2) );
legend(plt_svm,{str1});
hold on
plt_svm2 = plot(X(1,:),(-1*bs(1)/bs(2))*X(1,:)-b0s/bs(2) );
legend([plt_svm, plt_svm2],{str1,str2});
end