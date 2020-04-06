function min_value = optimize_vechL(vechL,Y,M,alpha, beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

expr1 = Y*Y';
% min_value =   alpha * expr1(:)'*M*vechL...
%             + beta * vechL' * (M' * M) * vechL; 
L_vec = M*vechL;
n_nodes =  sqrt(length(L_vec));

L = reshape(M*vechL,n_nodes,n_nodes);

min_value =   alpha * trace(Y'*L*Y)...
            + beta * norm(L); 


end

