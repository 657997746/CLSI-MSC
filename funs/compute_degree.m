function degree1 = compute_degree(weight)
m = length(weight);
    for i = 1:m
        degree1{i} = sum(weight{i},1);
        degree1{i} =  diag(degree1{i});
    end
end