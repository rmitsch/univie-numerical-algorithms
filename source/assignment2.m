##############################################
# For part 1.
##############################################

# Approximate condition number of A stochastically by generating random 
# vector and selecting the maximum ratio ||z|| / ||y|| to solve Az = y.
function condapprox = approximate_cond_stochastically(A, t = 1000)
    condapprox = 0;
    n = size(A)(1);

    max_ratio = 0;
    for i = 1:1:t
        y = rand(n, 1);
        z = A \ y;
        max_ratio = max(max_ratio, norm(z, 1) / norm(y, 1));
    end

    condapprox = norm(A, 1) * max_ratio

    printf("*****\n")
end

function executePart1()
    problem_sizes = [32 64] #[4 16 32 64 128 256 500 512];
    conds_approx = zeros(size(problem_sizes));
    conds = zeros(size(problem_sizes));
    condests_oct = zeros(size(problem_sizes));
    
    for i = 1:1:size(problem_sizes)(2)
        A = rand(problem_sizes(i));

        conds(i) = cond(A);
        condests(i) = condest(A, 1000);
        conds_approx(i) = approximate_cond_stochastically(A);
    end

    printf("------------\n")
    conds
    condests
    conds_approx
end

##############################################
# For part 2.
##############################################