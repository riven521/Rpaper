function [idxSolution,idxMin] = findMinIdxOfArray(PopSolutionArray)
%fitness BF�ȼ�����Ⱥÿ��Ⱦɫ���Fitness
% ����
%   PopSolutionArray    ��Ⱥ(��Fitness��Ŀ��ֵ)
% ���
%   idxMin     	        PopSolutionArray��������СFitness������

    objArray = arrayfun(@(x) sum(x.instance.Vehicle.Cost_Total), PopSolutionArray);
    [~,idxMin] = min(objArray);   %���ܴ��ڶ����Ϊ��Сֵ��idx, min����ȡ��һ��
    idxSolution = PopSolutionArray(idxMin);
    
end % End of findMinIdxOfArray

