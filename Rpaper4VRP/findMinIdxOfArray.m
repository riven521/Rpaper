function [idxSolution,idxMin] = findMinIdxOfArray(PopSolutionArray)
%fitness BF等计算种群每个染色体的Fitness
% 输入
%   PopSolutionArray    种群(含Fitness即目标值)
% 输出
%   idxMin     	        PopSolutionArray数组中最小Fitness索引号

    objArray = arrayfun(@(x) sum(x.instance.Vehicle.Cost_Total), PopSolutionArray);
    [~,idxMin] = min(objArray);   %可能存在多个均为最小值的idx, min函数取第一个
    idxSolution = PopSolutionArray(idxMin);
    
end % End of findMinIdxOfArray

