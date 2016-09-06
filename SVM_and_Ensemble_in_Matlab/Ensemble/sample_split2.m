function [hxi,sign_g] = sample_split2(X,y,s,n,j)

for i = 1:n
    if X(i,j) >= s
        hxi1(i,1) = 1;
    elseif X(i,j) < s
        hxi1(i,1) = -1;
    end
end

if (sum(hxi1~=y))/n > 0.5
    hxi = -1*hxi1;
    sign_g = -1;
else
    hxi = hxi1;
    sign_g = 1;
end

end

