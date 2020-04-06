% we are gonna write the algorithm
% this is from matlab

%iniitial setup
n_nodes = 11;
n_data = 10;
alpha = 0.1;
beta = 0.1;

L = optimvar('L',n_nodes,n_nodes);
Y = optimvar('Y',n_nodes,n_data);
prob = optimproblem;

f = @(L)trace(L);

%set the objective function
prob.Objective = fcn2optimexpr(f, L);

%% test the duplication matrix

test = full(DuplicationM(2))

vech = [1 2 3]';


%% write the objective function

%iniitial setup
n_nodes = 11;
n_data = 10;
alpha = 0.1;
beta = 0.1;

M = DuplicationM(n_nodes);

vechL = optimvar('vechL',n_nodes * (n_nodes + 1) / 2, 1);
Y = optimvar('Y',n_nodes,n_data);

prob = optimproblem;
prob.Objective = fcn2optimexpr(@optimize_vechL,vechL,Y,alpha,beta,'OutputSize',[1,1]);

%first tr(L) = n
expr_1 = @(vechL,M,n_nodes)trace(reshape(M * vechL, n_nodes,n_nodes));
optimization_cond_1 = fcn2optimexpr(expr_1,vechL,M,n_nodes);
prob.Constraints.const1 = optimization_cond_1 == n_nodes;

%second Lij  = Lji <= 0, non diagnal elements of the lower triangular
%matrix are smaller than 0
optimization_cond_2 = fcn2optimexpr(@cal_non_diag,vechL,n_nodes);
prob.Constraints.const2 = optimization_cond_2 <= zeros(n_nodes*(n_nodes-1)/2,1);

%third constraint: sum of each row is 0
expr_3 = @(vechL,n_nodes) reshape(M * vechL, n_nodes,n_nodes) * ones(n_nodes,1);
optimization_cond_3 = fcn2optimexpr(expr_3,vechL,n_nodes);

prob.Constraints.const3 = optimization_cond_3 == zeros(n_nodes,1);


