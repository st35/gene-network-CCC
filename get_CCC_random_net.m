function CCC = get_CCC_random_net(expression_data, N)
    CCC_rand = zeros(1, N);
    num_genes = size(expression_data, 2);
    num_pairs = nchoosek(num_genes, 2);

    corr_mat = corrcoef(expression_data);
    corr_mat = abs(corr_mat);

    corr_mat_original = corr_mat;

    for T = 1:N
        for i = 1:num_genes
            for j = (i + 1):num_genes
                RI = randi([1, num_genes]);
                RJ = randi([(i + 1), num_genes]);
                corr_mat(i ,j) = corr_mat_original(RI, RJ);
                corr_mat(j, i) = corr_mat_original(RI, RJ);
            end
        end
               
        d = sum(corr_mat, 2);
        D = diag(d);
        L = D - corr_mat;
        LI = pinv(L);
        V = sum(d);
        E = zeros(num_genes, num_genes);
        I = eye(num_genes, num_genes);
        for i = 1:num_genes
            for j = 1:num_genes
                E(i, j) = V*(transpose(I(:, i) - I(:, j))*LI*(I(:, i) - I(:, j)));
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
        CCC_rand(1, T) = cophenet(Z, Y);
    end

    CCC = mean(CCC_rand);
end