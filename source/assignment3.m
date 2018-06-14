% Generate matrix M for preconditioning of matrix A used in applyCG(...).
% precondConfig is a structure having some or all of the following properties:
%     - method is the preconditioning method to use. Available: diag, block, cholesky_nofi, cholesky_drop.
%     - dropThreshold is ignored for all methods except cholesky_drop.
%     - alpha is ignored for all methods except Cholesky factorization-based ones.
%     - stepsize is ignored for all methods except block-diagonale matrices.
% and contains configuration settins for preconditioning A.
function M = generatePreconditioning(A, precondConfig)
    switch (precondConfig.method)
        case "diag"
            % Extract diagonale matrix.
            M = diag(diag(A), 0);
            
        case "block"
            M = zeros(0, 0);
            
            % Extract and concatenate stepsized blocks to block-diagonale matrix.
            for i = 1:precondConfig.stepsize:size(A)(1)
                % Calculate range of indices to use for current block.
                block_idx_range = i:min(i + precondConfig.stepsize - 1, size(A)(1));
                % Assemble M by extending block-diagonale matrix with newly extracted block.
                M = blkdiag(M, A(block_idx_range, block_idx_range));
            end

        case "cholesky_nofi"
            chol_options = struct();
            chol_options.type = "nofill";
            chol_options.diagcomp = precondConfig.alpha;
            
            M = ichol(A, chol_options);
            
        case "cholesky_drop"
            chol_options = struct();
            chol_options.type = "ict";
            chol_options.diagcomp = precondConfig.alpha;
            chol_options.dropTol = precondConfig.dropThreshold;
            
            M = ichol(A, chol_options);
            
        % Wrong method name: Output warning, continue with M = I.
        otherwise
            printf("*** Warning *** Chosen method not supported.");
            M = eye();
    endswitch
end

% CG without preconditioning.
% Uses the relative residual as convergence criterion.
function [numberOfIterations, runTime, relResHistory] = applyCG(A, convergence_threshold = 10^-3, M_inverse = eye())
    tic;

    % Define x as vector of ones.
    x_true = ones(size(A)(1), 1);
    % Determine b implicitly.
    b_true = A * x_true;
    
    % Set up evaluation variables.
    numberOfIterations = 0;
    runTime = 0;
    relResHistory = [];
    
    % Initial guess for x.
    x = rand(size(A)(1), 1);
    % Compute initial residual.
    r = b_true - A * x;
    % Compute initial search direction.
    s = M_inverse * r;
    % Compute initial convergence criterion.
    convergence_criterion = norm(r) / norm(b_true);
    relResHistory = [relResHistory convergence_criterion];
    
    % Continue as long as convergence criterion exceeds defined threshold. 
    while convergence_criterion > convergence_threshold
        % Buffer value of previous r, since we need values of previous and current/next r at the same time.
        prev_r = r;
        
        % Precomputed values that are used more than once.
        As = A * s;
        prev_rtMr = transpose(prev_r) * M_inverse * prev_r;

        % Step size.
        alpha = prev_rtMr / (transpose(s) * As);
        % Next iteration of solution vector.
        x = x + alpha * s;
        % Next residual.
        r = prev_r - alpha * As;
        % Next beta.
        M_inverse_r = M_inverse * r;
        beta = (transpose(r) * M_inverse_r) / prev_rtMr;
        % Next search direction.
        s = M_inverse_r + beta * s;
        
        % Compute convergence criterion for current iteration.
        convergence_criterion = norm(A * x - b_true) / norm(b_true);
        relResHistory = [relResHistory convergence_criterion];
        
        % Keep track of number of iterations.
        numberOfIterations += 1;
    end

    runTime = toc;
end