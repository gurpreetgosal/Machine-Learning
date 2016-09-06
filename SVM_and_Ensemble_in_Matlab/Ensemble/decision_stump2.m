%% Implement ADABOOST with Decision Stump
function [L,s,hxi,sign_h,j_dir] = decision_stump2(X,y,n,d,w)

index1=1;
for i=1:n
    for j=1:d
        % define s
        s = X(i,j);
        % determine split
        [hxi,sign_g(index1)] = sample_split2(X,y,s,n,j);
        
        L_prov(index1,1)=(w'*(hxi~=y))/sum(w);
        note_s(index1,1) = s;
        note_hxi(:,index1) = hxi;
        index1=index1+1;
        dir_store(index1,1) = j;
    end
end

 [L, ind] = min(L_prov);
s = note_s(ind);
hxi = note_hxi(:,ind);
sign_h = sign_g(ind);
j_dir = dir_store(ind);


