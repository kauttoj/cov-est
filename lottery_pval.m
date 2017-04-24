function pval = lottery_pval(vec1,vec2)

    if nargin==0
        vec1 = randi(5,1,10000)>4;
        vec2 = randi(4,1,10000)>3;
    end

    M = 10000;
    
    vec1 = vec1~=0;
    vec2 = vec2~=0;
    N = length(vec1);
    NN = sum(vec1);
        
    n = sum(vec1.*vec2);
    
    null = zeros(1,M);
    for i=1:M
        vec1 = 0*vec1;
        vec1(randsample(N,NN,'false'))=1;
        null(i) = sum(vec1.*vec2);
    end
    
    pval = nnz(null>=n)/length(null);

end