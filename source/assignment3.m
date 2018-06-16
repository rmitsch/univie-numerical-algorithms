% Generate matrix M for preconditioning of matrix A used in applyCG(...).
% precondConfig is a structure having some or all of the following properties:
%     - method is the preconditioning method to use. Available: diag, block, cholesky_nofi, cholesky_drop.
%     - dropThreshold is ignored for all methods except cholesky_drop.
%     - alpha is ignored for all methods except Cholesky factorization-based ones.
%     - stepsize is ignored for all methods except block-diagonale matrices.
% and contains configuration settins for preconditioning A.
function M_inv = generatePreconditioning(A, precondConfig)
    switch (precondConfig.method)
        case "diag"
            % Extract diagonale matrix.
            M = diag(diag(A), 0);
            M_inv = inv(M);
            
        case "block"
            M = zeros(0, 0);
            
            % Extract and concatenate stepsized blocks to block-diagonale matrix.
            for i = 1:precondConfig.stepsize:size(A)(1)
                % Calculate range of indices to use for current block.
                block_idx_range = i:min(i + precondConfig.stepsize - 1, size(A)(1));
                % Assemble M by extending block-diagonale matrix with newly extracted block.
                M = blkdiag(M, A(block_idx_range, block_idx_range));
            end
            M_inv = inv(M);

        case "cholesky_nofi"
            chol_options = struct();
            chol_options.type = "nofill";
            chol_options.diagcomp = precondConfig.alpha;
            
            L = ichol(A, chol_options);
            L_inv = inv(L);
            M_inv = transpose(L_inv) * L_inv; 
            
        case "cholesky_drop"
            chol_options = struct();
            chol_options.type = "ict";
            chol_options.diagcomp = precondConfig.alpha;
            chol_options.droptol = precondConfig.dropThreshold;
            
            L = ichol(A, chol_options);
            L_inv = inv(L);
            M_inv = transpose(L_inv) * L_inv; 
            
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

% Investigate effect of threshold for PC with 
% incomplete Cholesky factorization with threshold 
% dropping.
function executeParameterGridSearch()
    starttime = time();
    printf("Starting grid search.\n")
    
    % Define matrix file names.
    filenames = {"nos5.mtx"; "nos6.mtx"; "s3rmt3m3.mtx"};
    % Determined by line search.
    fileAlphas = [0.5 0 774];
    
    % Create structure for holding results.
    res = struct();
    
    % ----------------------------------------
    % 1. Grid search for block-diagonale M.
    % ----------------------------------------
    % {
    printf("\n\n* Grid seach for block-diagonale (block size).")
    
    precondConfig = struct();
    precondConfig.method = "block";
    
    
    % Iterate over test matrices.
    for i = 1:length(filenames)
        printf("\n  - %s", filenames{i});
        disp(" ");
        fflush(stdout);
        
        res.(filenames{i}) = struct();
        res.(filenames{i}).numberOfIterations = [];
        res.(filenames{i}).runtime = [];
        res.(filenames{i}).relativeResidualHistory = [];
        
        % Load file.
        [A, rows, cols, entries, rep, field, symm] = mmread(["data/" filenames{i}]);        
        
        % Compute CG with block-diag. M with different block sizes.
        stepsizes = [2 4 8 16 32 64:round(size(A)(1) / 20):size(A)(1) / 4];
        for stepsize = stepsizes
            precondConfig.stepsize = stepsize;
            M_inv = generatePreconditioning(A, precondConfig);
            [numIter, runtime, relResHistory] = applyCG(A, 10^-8, M_inv);
            res.(filenames{i}).numberOfIterations = [res.(filenames{i}).numberOfIterations numIter];
            res.(filenames{i}).runtime = [res.(filenames{i}).runtime runtime];
        end
        
        [opt_runtime opt_idx] = min(res.(filenames{i}).runtime);
        printf("    Optimal block size found: %f", stepsizes(opt_idx));
    end
    % }
    
    % ----------------------------------------
    % 2. Grid search for incomplete Cholesky 
    % factorization with no fill-in.
    % ----------------------------------------
    printf("\n\n* Grid seach for incomplete Cholesky factorization with no fill-in (alpha).")
    
    precondConfig = struct();
    precondConfig.method = "cholesky_nofi";
    
    % {
    % Iterate over test matrices.
    for i = 1:length(filenames)
        printf("\n  - %s", filenames{i})
        disp(" ");
        fflush(stdout);
    
        res.(filenames{i}) = struct();
        res.(filenames{i}).numberOfIterations = [];
        res.(filenames{i}).runtime = [];
        res.(filenames{i}).relativeResidualHistory = [];
        
        % Load file.
        [A, rows, cols, entries, rep, field, symm] = mmread(["data/" filenames{i}]);
        
        % Compute CG with block-diag. M with different block sizes.
        % For calculation of upperAlphaLimit: https://www.mathworks.com/help/matlab/ref/ichol.html#btuq275-3.
        upperAlphaLimit = abs(max(sum(abs(A), 2) ./ diag(A)) - 2);
        interval = upperAlphaLimit / 5;
        alphas = [0 0.0001 0.005 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2 3 4 5 100 upperAlphaLimit 500 1000]; #interval:interval:upperAlphaLimit;
        for alpha = alphas
            precondConfig.alpha = alpha;
            
            M_inv = generatePreconditioning(A, precondConfig);
            [numIter, runtime, relResHistory] = applyCG(A, 10^-8, M_inv);
            res.(filenames{i}).numberOfIterations = [res.(filenames{i}).numberOfIterations numIter];
            res.(filenames{i}).runtime = [res.(filenames{i}).runtime runtime];
        end

        [opt_runtime opt_idx] = min(res.(filenames{i}).runtime);
        printf("    Optimal alpha found: %f", alphas(opt_idx));
        printf("    Runtime: %f", opt_runtime);
    end    
    % }
    
    % ----------------------------------------
    % 3. Grid search for incomplete Cholesky 
    % factorization with drop threshold for 
    % fill-in.
    % ----------------------------------------
    % {
    printf("\n\n* Grid seach for incomplete Cholesky factorization with fill-in dropping threshold (alpha, dropping threshold).")
    
    precondConfig = struct();
    precondConfig.method = "cholesky_drop";
    
    % Iterate over test matrices.
    for i = 1:length(filenames)
        printf("\n  - %s", filenames{i})
        fflush(stdout);
    
        res.(filenames{i}) = struct();
        res.(filenames{i}).numberOfIterations = [];
        res.(filenames{i}).runtime = [];
        res.(filenames{i}).relativeResidualHistory = [];
        
        % Load file.
        [A, rows, cols, entries, rep, field, symm] = mmread(["data/" filenames{i}]);        
        
        % Set alpha (determined in previous step).
        precondConfig.alpha = fileAlphas(i);
        
        % Compute CG with Cholesky factorization and fixed alpha. Vary dropping threshold.
        % Note: Dropping threshold is used as: abs (L(i,j)) >= droptol * norm (A(j:end, j), 1)
        thresholds = [0.1 0.2 0.3 0.4 0.5 1 2 3 5 10 15 20];
        for threshold = thresholds
            precondConfig.dropThreshold = threshold;
            printf("\n      %f", threshold);
            disp("    ")
            fflush(stdout);
            
            % Average over 5 runs.
            iterCount = 0;
            runtime = 0;
            repNumber = 3;
            for j = 1:repNumber
                startTime = time();
                M_inv = generatePreconditioning(A, precondConfig);
                [numIter, runtime, relResHistory] = applyCG(A, 10^-8, M_inv);
                iterCount += numIter;
                runtime += time() - startTime;
            end
            res.(filenames{i}).numberOfIterations = [res.(filenames{i}).numberOfIterations iterCount / repNumber];
            res.(filenames{i}).runtime = [res.(filenames{i}).runtime runtime / repNumber];
        end
        
        [opt_runtime opt_idx] = min(res.(filenames{i}).runtime);
        printf("    Optimal drop treshold found: %f", thresholds(opt_idx));
        
        % -----------------------------------
        % Plot results.
        % -----------------------------------
        
        % Thresholds vs. iterations.
        figure('Position',[0, 0, 450, 350])
        grid on
        hold on
        plot(thresholds, res.(filenames{i}).numberOfIterations, 'markersize', 3, '3; Thresholds vs. iterations;o-');
        title (["" filenames{i} ": Comparison of thresholds and iterations"], 'fontsize', 16);
        xlabel('Threshold');
        ylabel('Number of iterations');   
        % Thresholds vs. runtimes.
        figure('Position',[0, 0, 450, 350])
        grid on
        hold on
        plot(thresholds, res.(filenames{i}).runtime, 'markersize', 3, '2; Thresholds vs. runtimes;x-');
        title (["" filenames{i} ": Comparison of thresholds and runtimes"], 'fontsize', 16);
        xlabel('Threshold');
        ylabel('Runtime');           
    end
    % }
    
    printf("\n\nGrid search done. Time: %i\n", time() - starttime);
end