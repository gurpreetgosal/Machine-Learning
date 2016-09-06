function yhat = classify(Xtest,b,b0)

Xtest = Xtest';
tsamples = length(Xtest(:,1));
yhat = zeros(tsamples,1);

for i =1:tsamples
    if b*Xtest(i,:)'+b0>0 
        yhat(i) = 1;
    elseif b*Xtest(i,:)'+b0<0
         yhat(i) = -1;
    end
end
end