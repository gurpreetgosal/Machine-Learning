function [b, b0] = SoftMarg(X,y,gamma)

X= X';
K = X*X'; % linear kernel
ell = size(X,1);

f = -ones(ell,1);
A_ineq = -eye(ell);
a_ineq = zeros(ell,1);
A_eq = zeros(ell,ell);
A_eq(1,:) = y';
a_eq = zeros(ell,1);

H = (K).*(y*y');

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

end