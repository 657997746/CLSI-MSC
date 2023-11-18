function [grps,max_err] = CLSI_MSC(X,gt)
choice_graph = 2;
choice_metric = 1;
m = length(X);  %denotes the number of views
clusters = length(unique(gt));
tol = 10^-8;
maxmu = 10^6;
rou = 1.5;
maxiter = 100;
iter = 1;
mu = 10^-2;
%% need adjusted of parameters
p = 0.5;
lambda1 = 0.0001;
lambda2 = 0.001;
lambda3 = 0.0001;
%% Initialization
N = size(X{1},2);
C = zeros(N,N);
K = zeros(N,N);
Q = zeros(N,N);
Y4 = zeros(N,N);
sumXj = zeros(N,N);
sum_C_right = zeros(N,N);
P = zeros(N,N);
P_hat = zeros(N,N);
for i = 1:m
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
    W{i} = zeros(N,N);
    E{i} = zeros(size(X{i}));
    Y1{i} = zeros(size(X{i}));
    Y2{i} = zeros(N,N);
    Y3{i} = zeros(N,N);
    M{i} = zeros(N,N);
    D{i} = zeros(N,N);
    x(i) = norm(X{i},'fro')^2;
end
xx = sum(x);
IsConverge = 0;

%% compute Lap Martix
weight = construct_weight(m,X,choice_graph,choice_metric);
degree = compute_degree(weight);
for i = 1:m
    Lap{i} = degree{i}-weight{i};
end
%% Update
grps = zeros(size(length(gt),1));
while (IsConverge == 0 && iter<maxiter+1)
    sum_Dv = 0;
    for j = 1:m
        %% W
         W{j} = (X{j}'*X{j}+2*eye(N))\(X{j}'*(X{j}-X{j}*C-E{j}+Y1{j}/mu)-(C-M{j}+Y2{j}/mu)-(diag(diag(D{j}))-D{j}+Y3{j}/mu));
        %% update E
         temp = X{j}-X{j}*(C+W{j})+Y1{j}/mu;
         E{j} = solve_l1l2(temp,1/mu);
        %% update D
         U{j} = W{j} + Y3{j}./mu;
        sum_Dv = sum_Dv+abs(D{j});
         D{j} = max(0,U{j}-(lambda1*sum_Dv)/mu).*sign(U{j});%+min(0,U{j}+(lambda2*sum_Dv)/mu);
         D{j} = D{j}-diag(diag(D{j}));
        %% M
        M{j} = mu*(C+W{j}+Y2{j}/mu)/(2*lambda2*Lap{j}+mu*eye(N));
        M{j} = max(M{j},0);
    end
    %% update C
    for j = 1:m
        sumXj = sumXj+X{j}'*X{j}+eye(N); 
        sum_C_right = sum_C_right +X{j}'*(X{j}-X{j}*W{j}-E{j}+Y1{j}/mu)-(W{j}-M{j}+Y2{j}/mu);
    end 
        C = (sumXj+eye(N))\(sum_C_right+K-Y4/mu);
    %%      update K
    [U1,S1,V1] = svd(C+Y4/mu);
    s = diag(S1);
    lambda = lambda3/mu;
    for mm = 1:length(s)
        s1(mm) = shcatten_p(s(mm),lambda, p);
    end
    s1(s1<0) = 0;
    T1 = (diag(s1));
    K = U1*T1*V1';

        %%         update multiplier Y1 Y2 Y3 Y4
        for j = 1:m
            Y1{j} = Y1{j} + mu *(X{j}- X{j}*(C+D{j})-E{j});
            Y2{j} = Y2{j} + mu *(C+W{j}-M{j});
            Y3{j} = Y3{j} + mu *(W{j} - D{j} + diag(D{j}));
        end
        Y4 = Y4 + mu*(C-K);
        mu = min(maxmu,mu*rou);
    %% check if converge
    for j = 1:m
        err1(j) = norm(X{j}- X{j}*(C+D{j})-E{j},"fro");
        err2(j) = norm(C+W{j}-M{j},"fro");
        err3(j) = norm(W{j} - D{j} + diag(D{j}),"fro");
        f(j) = norm(X{j}- X{j}*(C+D{j})-E{j},'fro')^2;
    end
    Err(iter) = sum(err1);
    max_err0 = max([err1(:);err2(:);err3(:)]);
    max_err(iter) = max([max_err0,norm(C-K,"fro")]);
    if  max_err(iter) < tol
        IsConverge = 1;
    end
    iter = iter+1;
    Z1 = zeros(N,N);
    Z3 = zeros(N,N);
    for j = 1:m
        Z1 = Z1 + (abs(M{j})+abs(M{j})');
        Z3 = Z3 + (abs(D{j})+abs(D{j})');
    end
end
 % Affinity matrix
 Z2 = (abs(C)+abs(C)')/2+Z3/(2*m);
 Z =  Z1/(2*m); 
 ZZ = (abs(C)+abs(C)')/2;
 grps = SpectralClustering(Z,clusters); 
end