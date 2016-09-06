function [b, b0] = SoftMarg_gammaopt(X,y,gamma)

X= X';
K = X*X'; % linear kernel
ell = size(X,1);

f = -ones(ell,1);
A_ineq = -eye(ell);
a_ineq = zeros(ell,1);
A_eq = zeros(ell,ell);
A_eq(1,:) = y';
a_eq = zeros(ell,1);

% first find optimum gamma
index1 = 1;
for gamma = 1:2:300
    
    H = (K + eye(ell)/gamma).*(y*y');
    
    LB = zeros(ell,1);
    UB = gamma*ones(ell,1);
    alpha = quadprog(H,f,[],[],A_eq,a_eq,LB,UB);
    b = (alpha.*y)'*X;
    
    for i =1:length(alpha)
        if y(i)==1  && alpha(i)>0.1
            note_i = i;
        end
    end
    
    b0 = 1 - b*X(note_i,:)';
    
    % now determine class of each of the training sample and compute
    % misclassification error
    merror_abs = 0;
    for k = 1:ell
        if ((b*X(k,:)' + b0< 0) && (y(k)==1)) || ((b*X(k,:)' + b0 > 0) && (y(k)==-1))
            merror_abs = merror_abs+1;
        end
    end
    merror(index1,1) = merror_abs/ell;
    b_store(index1,:) = b;
    b0_store(index1) = b0;
    gamma_s(index1)=gamma;
    index1 = index1 +1;
    
end
[~,ind_min_error] = min(merror);

b=b_store(ind_min_error,:);
b0 = b0_store(ind_min_error);

end