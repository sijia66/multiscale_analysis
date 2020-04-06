function [L_iter,Y,vechl_result, Y_result] = cal_GLSigRep(X_test, n_data, Y_0, M, alpha, num_iter, Y, beta, n_nodes, vechL_0)
tic
for iter_i = 1:num_iter
    
    %optimize L
    vechL = optimvar('vechL', n_nodes * (n_nodes + 1) / 2, 1);
    
    prob = optimproblem;
    prob.Objective = fcn2optimexpr(@optimize_vechL,vechL,Y,M,alpha,beta,'OutputSize',[1,1]);
    
    %first tr(L) = n
    expr_1 = @(vechL,M,n_nodes)trace(reshape(M * vechL, n_nodes,n_nodes));
    optimization_cond_1 = fcn2optimexpr(expr_1,vechL,M,n_nodes);
    prob.Constraints.const1 = optimization_cond_1 == n_nodes;
    
    %second Lij  = Lji <= 0, non diagnal elements of the lower triangular
    %matrix are smaller than 0
    optimization_cond_2 = fcn2optimexpr(@cal_non_diag,vechL,n_nodes) <= zeros(n_nodes*(n_nodes-1)/2,1);
    prob.Constraints.const2 = optimization_cond_2 ;
    
    %third constraint: sum of each row is 0
    expr_3 = @(vechL,n_nodes) reshape(M * vechL, n_nodes,n_nodes) * ones(n_nodes,1);
    optimization_cond_3 = fcn2optimexpr(expr_3,vechL,n_nodes); 
    prob.Constraints.const3 = optimization_cond_3 == zeros(n_nodes,1);
    
    solution = solve(prob,vechL_0); %solve the problem
    L_iter = reshape(M * solution.vechL,n_nodes,n_nodes);
    
    %second part of the problem is to optimize L_iter
    Y_opt = optimvar('Y_opt',n_nodes,n_data);
    prob2 = optimproblem;
    prob2.Objective = fcn2optimexpr(@optimize_Y,Y_opt,X_test,L_iter,alpha,'OutputSize',[1,1]);
    
    %initialize the problem
    solution2 = solve(prob2,Y_0); %solve the problem
    
    
    %update Y for next round
    vechL_0.vechL = solution.vechL;
    Y = solution2.Y_opt;
    Y_0.Y_opt = Y;
    
    %save the results
    vechl_result(iter_i) = optimize_vechL(solution.vechL,Y,M,alpha, beta);
    Y_result(iter_i) = optimize_Y(Y,X_test,L_iter,alpha);
    
    fprintf('Finished %d iteration with %0.2f  and %0.2f',iter_i,vechl_result(iter_i),Y_result(iter_i))
    
end
end
