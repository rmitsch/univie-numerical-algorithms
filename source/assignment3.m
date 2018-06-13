% CG without preconditioning.
% Uses the relative residual as convergence criterion.
function [numberOfIterations, runTime, relResHistory] = applyCG(A, convergence_threshold)
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
    s = r;
    % Compute initial convergence criterion.
    convergence_criterion = norm(r) / norm(b_true);
    relResHistory = [relResHistory convergence_criterion];
    
    % Continue as long as convergence criterion exceeds defined threshold. 
    while convergence_criterion > convergence_threshold
        % Buffer value of previous r, since we need values of previous and current/next r at the same time.
        prev_r = r;
        
        % Precomputed values that are used more than once.
        As = A * s;
        prev_rtr = transpose(prev_r) * prev_r;
        
        % Step size.
        alpha = prev_rtr / (transpose(s) * As);
        % Next iteration of solution vector.
        x = x + alpha * s;
        % Next residual.
        r = prev_r - alpha * As;
        % Next beta.
        beta = (transpose(r) * r) / prev_rtr;
        % Next search direction.
        s = r + beta * s;
        
        % Compute convergence criterion for current iteration.
        convergence_criterion = norm(A * x - b_true) / norm(b_true);
        relResHistory = [relResHistory convergence_criterion];
        
        % Keep track of number of iterations.
        numberOfIterations += 1;
    end
    
    runTime = toc;
end