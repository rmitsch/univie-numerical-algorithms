%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For part 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate condition number of A stochastically by generating random 
% vector and selecting the maximum ratio ||z|| / ||y|| to solve Az = y.
function condapprox = approximate_cond_stochastically(A, t = 1000)
    condapprox = 0;
    n = size(A)(1);

    max_ratio = 0;
    for i = 1:1:t
        y = rand(n, 1);
        z = A \ y;
        max_ratio = max(max_ratio, norm(z, 1) / norm(y, 1));
    end

    condapprox = norm(A, 1) * max_ratio;
end

% Execute code for part 1 of the programming assignment.
function executePart1()
	% --------------------------------
	% Part 1 of the task: Investigate
	% different n. 
	% --------------------------------

    problem_sizes = [5 10 25 50 100 200 300 400 500];
    conds_approx = zeros(size(problem_sizes));
    conds = zeros(size(problem_sizes));
    condests_oct = zeros(size(problem_sizes));
    
    for i = 1:1:size(problem_sizes)(2)
        A = rand(problem_sizes(i));

        conds(i) = cond(A, 1);
        condests(i) = condest(A);

        for j = 1:1:5
        	conds_approx(i) += approximate_cond_stochastically(A);
    	end
    	conds_approx(i) /= 5;
    end

    % Plot results.
	figure('Position',[0, 0, 750, 350])
	grid on
	hold on

	semilogy(problem_sizes, conds, 'markersize', 3, '2; cond(A);x-');
	semilogy(problem_sizes, condests, 'markersize', 3, '3; condest(A);o-');
	semilogy(problem_sizes, conds_approx, 'markersize', 3, '4; Randomized cond. number approx.;.-');

	legend ({
	    'cond(A) ',
	    'condest(A) ',
	    'Rand. approximation of cond(A) '
	}, 'location', 'eastoutside')
	title ("Fig. 2: Comparison of condition estimates \nfor random matrices", 'fontsize', 16);
	xlabel('n');
	ylabel('(Estimation of) condition number');

	% --------------------------------
	% Part 2 of the task: Investigate
	% different number of tries. 
	% --------------------------------

    problem_sizes = [5 10 100];
	random_tries = [1 5 10 20 50 100 200 300 500 1000 2000 5000];
	rel_cond_approx_errors = zeros(size(problem_sizes)(2), size(random_tries)(2));
	
	for i = 1:size(problem_sizes)(2)
		A = rand(problem_sizes(i));
		condA = cond(A, 1);

		for j = 1:size(random_tries)(2)
			% Calculate averaged randomized approximation to condition number.
			cond_approx = 0;
	        for k = 1:1:5
	        	cond_approx += approximate_cond_stochastically(A, random_tries(j));
	    	end
	    	cond_approx /= 5;

			rel_cond_approx_errors(i, j) = abs(condA - cond_approx) / condA;
		end
    end

    % Plot results.
	figure('Position',[0, 0, 750, 350])
	grid on
	hold on

	plot(random_tries, rel_cond_approx_errors(1,:), 'markersize', 3, '2; n = 5;x-');
	plot(random_tries, rel_cond_approx_errors(2,:), 'markersize', 3, '3; n = 10;o-');
	plot(random_tries, rel_cond_approx_errors(3,:), 'markersize', 3, '4; n = 100;.-');

	legend ({
	    'n = 5 ',
	    'n = 10 ',
	    'n = 100 '
	}, 'location', 'eastoutside')
	title ("Fig. 3: Behaviour of random. approximation of cond with growing number of tries", 'fontsize', 16);
	xlabel('Number of randomly selected vectors y');
	ylabel('Relative error: |cond(A) - cond_{approx}(A)| / cond(A)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For part 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execute code for part 2 of the programming assignment.
function executePart2()
	problem_sizes = [100:50:1500];
	number_of_repetitions = 10;
	lhs_results = zeros(size(problem_sizes));
	rhs_results = zeros(size(problem_sizes));

	for i = 1:size(problem_sizes)(2)
		% Define x as vector of ones.
		x = ones(problem_sizes(i), 1);

		% Generate random delta b.
		delta_b = rand(problem_sizes(i), 1);
		% Scale delta_b so that norm1(delta_b) = 10^-8. Since it's the 1-norm, we can just multiply with the 
		% corresponding ratio factor.
		delta_b = delta_b / (norm(delta_b, 1) / 10^-8);
		
		% Generate random input perturbation E.
		E = rand(problem_sizes(i));
		% Scale E so that norm1(delta_b) = 10^-8. See above (scaling for delta_b).
		E = E / (norm(E, 1) / 10^-8);

		% Average results for both sides over k iterations.
		for k = 1:1:number_of_repetitions
			% Generate random A.
			A = rand(problem_sizes(i));

			% Assume b to be correct after this calculation.
			b = A * x;
			
			% Calculate delta_x as A \ (b + delta_b).
			x_hat = A \ (b + delta_b);
			delta_x = abs(x - x_hat);

			% Compute left hand side.
			lhs_results(i) += norm(delta_x, 1) / norm(x, 1);

			% Compute right hand side.
			rhs_results(i) += cond(A) * (norm(delta_b, 1) / norm(b, 1) + norm(E, 1) / norm(A, 1));
		end
		lhs_results(i) /= number_of_repetitions;
		rhs_results(i) /= number_of_repetitions;
	end	

    % Plot results.	
	figure('Position',[0, 0, 750, 350])
	grid on
	hold on

	semilogy(problem_sizes, lhs_results, 'markersize', 3, '2; n = 5;x-');
	semilogy(problem_sizes, rhs_results, 'markersize', 3, '3; n = 10;o-');
	
	legend ({
	    'Left hand side (actual rel. error)',
	    'Right hand side (upper bound for rel. error)',
	}, 'location', 'eastoutside')
	title ("Fig. 4: Bounds of the relative error \n in a linear system in practice", 'fontsize', 16);
	xlabel('n');
	ylabel('Relative error (actual and upper bound)');
end