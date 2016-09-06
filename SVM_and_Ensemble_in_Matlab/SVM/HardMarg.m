function [b, b0] = HardMarg(X,y)
X= X';
ell = size(X,1);

f = -ones(ell,1);
% inequality constraints
A_ineq = -eye(ell);
a_ineq = zeros(ell,1);
% Equality constraints
A_eq = zeros(ell,ell);
A_eq(1,:) = y';
a_eq = zeros(ell,1);
% H matrix
H = (X*X').*(y*y');
% solve dual using quadprog
alpha = quadprog(H,f,A_ineq,a_ineq,A_eq,a_eq);
% determine beta
b = (alpha.*y)'*X;
% determine beta0
for i =1:length(alpha)
    if y(i)==1  && alpha(i)>0.1
        note_i = i;
    end
end

b0 = 1 - b*X(note_i,:)';

end