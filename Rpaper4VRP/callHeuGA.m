function PopSolutionMatrix = callHeuGA(thisInsData, Arg) 
%callHeuGA  调用GA算法   
% 输入
%   thisInsData     Ga输入算例
%   Arg             Ga算法参数
% 输出
%   PopSolutionMatrix  - 已计算Fitness的所有迭代种群Array 行:迭代数; 列:种群数;
%% function callHeuGA(thisInsData, Arg)
% Nested Functions:
% # getNewPopArray      GA迭代主体
% Local Functions:
% # initilizePopOfGA            555 获取GA算法初始种群
% # fitnessOfGA                 555 BF等启发式获取目标值作为每个解的fitness
% # findMinIdxOfArray(子函数)   找出种群中最小Fitness的索引号

% 0 GA参数传递 及初始化
nPop = 4*ceil(Arg.NPOP/4);   %必须是4的整数，如果是10，变为12
nCustomer = length(thisInsData.Customer.Demand);   %-个体编码长度 % nVehicle = length(thisInsData.Vehicle.Capacity);
PopSolutionMatrix = []; % 保留所有迭代的解 struct matrix

% 1 GA种群初始化 Population种群 包含所有解的信息 不只索引号
PopArray = initilizePopulation( thisInsData , nPop );

% 2  GA迭代开始
    % PopArray          - 未计算Fitness 给定顾客顺序 车辆顺序一般已某种固定顺序为单一
    % PopSolutionArray  - 已计算Fitness
    
for iter = 1: Arg.MAX_ITERATION+1
    PopSolutionArray = computeFitness(PopArray);  
    PopSolutionMatrix = [PopSolutionMatrix; PopSolutionArray];    
    
    if iter < Arg.MAX_ITERATION
        PopArray = getNewPopArray();
    else
        break;
    end    
end
   
%% % 嵌套 FUNCTIONS %
% ==========================================================================
% 嵌套 FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Nested function getNewPopArray()
% 
% 
    function newPopArray = getNewPopArray()
        % 迭代 (每次随机选取4个解；保留最好的解到下一次迭代;其余三个做三次变异; 无交叉; 无适应度选择)
        A = randperm(nPop); % A - 种群随机后,从中取4个; 取代选择操作
        for p = 4:4:nPop
            sub4PopSolutionArray = PopSolutionArray(A(p-3:p));  %sub4PopSolutionArray 随机取的4个染色体(含Fitness)
            sub4PopArray = PopArray(A(p-3:p));  %sub4PopArray 随机取的4个染色体(不含Fitness)
            
            [~,idxMin] = findMinIdxOfArray(sub4PopSolutionArray);
            bestOf4Route = sub4PopArray(idxMin);
            
            routeInsertionPoints = sort(ceil(nCustomer*rand(1,2)));
            I = routeInsertionPoints(1);
            J = routeInsertionPoints(2);
            
            bestCustomerIdx = bestOf4Route.Customer.Idx;
            for k = 1:4 % Mutate the Best to get Three New Routes
                tmpPop = bestCustomerIdx;
                switch k
                    case 1 %
                        newPopArray(p-3) = bestOf4Route; %- 新解的第一个解用于是上个迭代过程最优的 %Elite保留策略
                    case 2 % Flip
                        tmpPop(I:J) = bestCustomerIdx(J:-1:I);
                        newPopArray(p-2) = sortThisInsData...
                            (bestOf4Route,tmpPop,bestOf4Route.Vehicle.Idx);
                    case 3 % Swap
                        tmpPop([I J]) = bestCustomerIdx([J I]);
                        newPopArray(p-1) = sortThisInsData...
                            (bestOf4Route,tmpPop,bestOf4Route.Vehicle.Idx);
                    case 4 % Slide
                        tmpPop(I:J) = bestCustomerIdx([I+1:J I]);
                        newPopArray(p-0) = sortThisInsData...
                            (bestOf4Route,tmpPop,bestOf4Route.Vehicle.Idx);
                    otherwise % Do Nothing
                end
            end
        end
    end

end
           
%%  % 局部 FUNCTIONS %
% ==========================================================================
% 局部 FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

function PopArray = initilizePopulation(thisInsData, nPop)
%initilizePopOfGA 获取GA算法初始种群(策略选择:1车辆是否随机排序/2种群是否完全随机)
% 输入
%   thisInsData      instance算例
%   nPop             Ga种群数
% 输出
%   PopArray         instance不同顾客的种群算例
% 注意 : NOTE 两种初始化 1 全部随机 2 部分初始解由算法生成

% 初始化
idxCus = {};   idxVeh = {};
nCustomer = length(thisInsData.Customer.Demand); 
nVehicle = length(thisInsData.Vehicle.Capacity);

% 不同顾客排序
idx = (1: nCustomer )';              idxCus{end+1} = idx;     %获取正常的顺序
[~, idx] = sort(thisInsData.Customer.Demand,'descend');      idxCus{end+1} = idx;  %获取递减的V_size  最重要
[~, idx] = sort(thisInsData.Customer.Demand,'ascend');       idxCus{end+1} = idx; %获取递增的V_size 最重要

% 不同车辆排序（固定）
idx = (1 : nVehicle)';               idxVeh{end+1} = idx;     %获取正常的顺序

% 去除重复customer和vehicle索引顺序
idxCus = unique(cell2mat(idxCus)','rows')';  idxVeh = unique(cell2mat(idxVeh)','rows')';

% 重复一定次数 方便arrayfun使用
nIdxCus = size(idxCus,2); nIdxVeh = size(idxVeh,2); nTotal = nIdxCus * nIdxVeh;
idxCusArray = repelem(idxCus,1, nIdxVeh); idxVehArray = repelem(idxVeh,1, nIdxCus);
InsDataArray = repmat(thisInsData, [nTotal,1]);

for iSort = 1: nPop
    if iSort <= nTotal && 1 % TODO 增加初始化种群的参数 
        thisSortedInsData = sortThisInsData(InsDataArray(iSort), idxCusArray(:,iSort), idxVehArray(:,iSort));
    elseif iSort <= nPop
%         thisSortedInsData = sortThisInsData(thisInsData, randperm(nCustomer)', idxVeh );%车辆固定顺序 TODO 二选一
        thisSortedInsData = sortThisInsData(thisInsData,randperm(nCustomer)', randperm(nVehicle)' );%车辆随机顺序
    else
        error('EEEEE');
    end
    PopArray(iSort) =  thisSortedInsData;
end

end

function PopSolutionArray = computeFitness(PopArray)
%fitness BF等计算种群每个染色体的Fitness
% 输入
%   PopArray             种群(不含Fitness即目标值)
% 输出
%   PopSolutionArray     种群(包含Fitness即目标值)

PopSolutionArray = [];
for pIter = 1 : length(PopArray)
    PopSolutionArray = [PopSolutionArray callHeuBF(PopArray(pIter))];  %BF算法计算Fitness
    PopSolutionArray(pIter).algName = 'Genetic Algorithm';
    PopSolutionArray(pIter).algArg = 'TBA';
end
end % End of fitness