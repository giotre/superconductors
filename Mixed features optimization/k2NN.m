%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X-> Materiali with T > T_threshold
% Y-> Materiali with T < T_threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function difff=k2NN(X,Y,r)

count=0;

Nsi=max(size(X));

for i=1:Nsi

    xx=[X(i,1),X(i,2)];

    [idx,~] = rangesearch(Y,xx,r);

    if isempty(idx{1})==0
 
        count=count+max(size(idx{1}));

    end

end

difff=count/Nsi;  % Average number of negative neighbors around a positive sample.