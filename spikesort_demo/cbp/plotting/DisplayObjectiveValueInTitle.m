function DisplayObjectiveValueInTitle(iteration_num, ...
                                      ResidualNormProgress, ...
                                      SparsityProgress)
% Display the objective function value and the last percent change.
F = ResidualNormProgress(iteration_num) + SparsityProgress(iteration_num);
if (iteration_num == 1)
    title(sprintf('F = %0.3f',F));
else
    F0 = ResidualNormProgress(iteration_num - 1) + ...
         SparsityProgress(iteration_num - 1);
    dF = - (F - F0) / abs(F0) * 100;
    title(sprintf('F = %0.3f Per. Decr = %0.3f', F, dF));          
end