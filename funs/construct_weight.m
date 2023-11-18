function weight = construct_weight(m,X,choice_graph,choice_metric)
%% Constructing the SIG matrices  构造相似性诱导图矩阵
pn = 15; % pn: number of adaptive neighbours  自适应的邻居数量
options = [];
options.k = 5;

weight = cell(1,m);  %生成一个1行m列的cell用来存放每个视角生成的矩阵S0
for i = 1:m
    if 1 == choice_graph % complete graph
        options.k = 0;
        if 1 == choice_metric
            options.WeightMode = 'Binary';
            weight{i} = constructS_KNG(X{i}', options);
        elseif 2 == choice_metric
            options.WeightMode = 'Cosine';
            weight{i} = constructS_KNG(X{i}', options);
        elseif 3 == choice_metric
            options.WeightMode = 'HeatKernel';
            weight{i} = constructS_KNG(X{i}', options);
        else
            error('Invalid input: check choice_metric');
        end
    elseif 2 == choice_graph % k-nearest graph
        if 1 == choice_metric
            options.WeightMode = 'Binary';
            weight{i} = constructS_KNG(X{i}', options);
        elseif 2 == choice_metric
            options.WeightMode = 'Cosine';
            weight{i} = constructS_KNG(X{i}', options);
        elseif 3 == choice_metric
            options.WeightMode = 'HeatKernel';
            weight{i} = constructS_KNG(X{i}', options);
        elseif 4 == choice_metric
            [weight{i}, distX_i] = constructS_PNG(X{i}, pn, 0);
        else
            error('Invalid input: check choice_metric');
        end
    else
        error('Invalid input: check choice_graph');
    end
end
end
