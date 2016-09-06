function mcl_err = compute_error(Xtest,ytest,b,b0)

yhat = classify(Xtest,b,b0);
error = 0;
for i=1:length(ytest)
    if ytest(i)~=yhat(i)
        error=error+1;
    end
end
mcl_err = error/length(ytest);

end