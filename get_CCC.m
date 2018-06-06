function CCC = get_CCC(expression_data)
    num_genes = size(expression_data, 2);
    num_pairs = nchoosek(num_genes, 2);
    
    corr_mat = corrcoef(expression_data);
    corr_mat = abs(corr_mat);
    
    d = sum(corr_mat, 2);
    D = diag(d);
    L = D - corr_mat;
    LI = pinv(L);
    V = sum(d);
    E = zeros(num_genes, num_genes);
    I = eye(num_genes, num_genes);
    for i = 1:num_genes
        for j = 1:num_genes
            E(i, j)=V*(transpose(I(:, i) - I(:, j))*LI*(I(:, i) - I(:, j)));
        end
    end
    E = sqrt(E);
    Y = zeros(1, num_pairs);
    n = 1;
    for i = 1:(num_genes - 1)
        for j = (i + 1):num_genes
            Y(1, n) = E(i, j);
            n = n + 1;
        end
    end
    Z = linkage(Y);
    CCC = cophenet(Z, Y);
end