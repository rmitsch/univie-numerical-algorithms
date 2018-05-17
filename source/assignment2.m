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

    condapprox = norm(A, 1) * max_ratio;
end

function executePart1()
	# --------------------------------
	# Part 1 of the task: Investigate
	# different n. 
	# --------------------------------

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

    # Plot results.
	figure('Position',[0, 0, 750, 350])
	grid on
	hold on

	semilogy(problem_sizes, conds, '2; cond(A);x-');
	semilogy(problem_sizes, condests, "markersize", 3, '3; condest(A);o-');
	semilogy(problem_sizes, conds_approx, "markersize", 3, '4; Randomized cond. number approx.;.-');

	legend ({
	    "cond(A) ",
	    "condest(A) ",
	    "Rand. approximation of cond(A) "
	}, "location", "eastoutside")
	title ("Fig. 2: Comparison of condition estimates \nfor random matrices", "fontsize", 16);
	xlabel("n");
	ylabel("(Estimation of) condition number");

	# --------------------------------
	# Part 2 of the task: Investigate
	# different number of tries. 
	# --------------------------------

    problem_sizes = [5 10 100];
	random_tries = [1 5 10 20 50 100 200 300 500 1000 2000 5000];
	rel_cond_approx_errors = zeros(size(problem_sizes)(2), size(random_tries)(2));
	
	for i = 1:size(problem_sizes)(2)
		A = rand(problem_sizes(i));
		condA = cond(A, 1);

		for j = 1:size(random_tries)(2)
			# Calculate averaged
			cond_approx = 0;
	        for k = 1:1:5
	        	cond_approx += approximate_cond_stochastically(A, random_tries(j));
	    	end
	    	cond_approx /= 5;

			rel_cond_approx_errors(i, j) = abs(condA - cond_approx) / condA;
		end
    end

    # Plot results.
	figure('Position',[0, 0, 750, 350])
	grid on
	hold on

	plot(random_tries, rel_cond_approx_errors(1,:), '2; n = 5;x-');
	plot(random_tries, rel_cond_approx_errors(2,:), "markersize", 3, '3; n = 10;o-');
	plot(random_tries, rel_cond_approx_errors(3,:), "markersize", 3, '4; n = 100;.-');

	legend ({
	    "n = 5 ",
	    "n = 10 ",
	    "n = 100 "
	}, "location", "eastoutside")
	title ("Fig. 3: Behaviour of random. approximation of cond with growing number of tries", "fontsize", 16);
	xlabel("Number of randomly selected vectors y");
	ylabel("Relative error: |cond(A) - cond_{approx}(A)| / cond(A)");

end

##############################################
# For part 2.
##############################################