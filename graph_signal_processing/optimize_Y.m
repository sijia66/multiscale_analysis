function min_value = optimize_Y(Y,X,L,alpha)
min_value = norm(X - Y) + alpha * trace(Y'*L * Y);
end