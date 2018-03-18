function PopSolutionMatrix = callHeuGA(thisInsData, Arg) 
%callHeuGA  ����GA�㷨   
% ����
%   thisInsData     Ga��������
%   Arg             Ga�㷨����
% ���
%   PopSolutionMatrix  - �Ѽ���Fitness�����е�����ȺArray ��:������; ��:��Ⱥ��;
%% function callHeuGA(thisInsData, Arg)
% Nested Functions:
% # getNewPopArray      GA��������
% Local Functions:
% # initilizePopOfGA            555 ��ȡGA�㷨��ʼ��Ⱥ
% # fitnessOfGA                 555 BF������ʽ��ȡĿ��ֵ��Ϊÿ�����fitness
% # findMinIdxOfArray(�Ӻ���)   �ҳ���Ⱥ����СFitness��������

% 0 GA�������� ����ʼ��
nPop = 4*ceil(Arg.NPOP/4);   %������4�������������10����Ϊ12
nCustomer = length(thisInsData.Customer.Demand);   %-������볤�� % nVehicle = length(thisInsData.Vehicle.Capacity);
PopSolutionMatrix = []; % �������е����Ľ� struct matrix

% 1 GA��Ⱥ��ʼ�� Population��Ⱥ �������н����Ϣ ��ֻ������
PopArray = initilizePopulation( thisInsData , nPop );

% 2  GA������ʼ
    % PopArray          - δ����Fitness �����˿�˳�� ����˳��һ����ĳ�̶ֹ�˳��Ϊ��һ
    % PopSolutionArray  - �Ѽ���Fitness
    
for iter = 1: Arg.MAX_ITERATION+1
    PopSolutionArray = computeFitness(PopArray);  
    PopSolutionMatrix = [PopSolutionMatrix; PopSolutionArray];    
    
    if iter < Arg.MAX_ITERATION
        PopArray = getNewPopArray();
    else
        break;
    end    
end
   
%% % Ƕ�� FUNCTIONS %
% ==========================================================================
% Ƕ�� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Nested function getNewPopArray()
% 
% 
    function newPopArray = getNewPopArray()
        % ���� (ÿ�����ѡȡ4���⣻������õĽ⵽��һ�ε���;�������������α���; �޽���; ����Ӧ��ѡ��)
        A = randperm(nPop); % A - ��Ⱥ�����,����ȡ4��; ȡ��ѡ�����
        for p = 4:4:nPop
            sub4PopSolutionArray = PopSolutionArray(A(p-3:p));  %sub4PopSolutionArray ���ȡ��4��Ⱦɫ��(��Fitness)
            sub4PopArray = PopArray(A(p-3:p));  %sub4PopArray ���ȡ��4��Ⱦɫ��(����Fitness)
            
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
                        newPopArray(p-3) = bestOf4Route; %- �½�ĵ�һ�����������ϸ������������ŵ� %Elite��������
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
           
%%  % �ֲ� FUNCTIONS %
% ==========================================================================
% �ֲ� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

function PopArray = initilizePopulation(thisInsData, nPop)
%initilizePopOfGA ��ȡGA�㷨��ʼ��Ⱥ(����ѡ��:1�����Ƿ��������/2��Ⱥ�Ƿ���ȫ���)
% ����
%   thisInsData      instance����
%   nPop             Ga��Ⱥ��
% ���
%   PopArray         instance��ͬ�˿͵���Ⱥ����
% ע�� : NOTE ���ֳ�ʼ�� 1 ȫ����� 2 ���ֳ�ʼ�����㷨����

% ��ʼ��
idxCus = {};   idxVeh = {};
nCustomer = length(thisInsData.Customer.Demand); 
nVehicle = length(thisInsData.Vehicle.Capacity);

% ��ͬ�˿�����
idx = (1: nCustomer )';              idxCus{end+1} = idx;     %��ȡ������˳��
[~, idx] = sort(thisInsData.Customer.Demand,'descend');      idxCus{end+1} = idx;  %��ȡ�ݼ���V_size  ����Ҫ
[~, idx] = sort(thisInsData.Customer.Demand,'ascend');       idxCus{end+1} = idx; %��ȡ������V_size ����Ҫ

% ��ͬ�������򣨹̶���
idx = (1 : nVehicle)';               idxVeh{end+1} = idx;     %��ȡ������˳��

% ȥ���ظ�customer��vehicle����˳��
idxCus = unique(cell2mat(idxCus)','rows')';  idxVeh = unique(cell2mat(idxVeh)','rows')';

% �ظ�һ������ ����arrayfunʹ��
nIdxCus = size(idxCus,2); nIdxVeh = size(idxVeh,2); nTotal = nIdxCus * nIdxVeh;
idxCusArray = repelem(idxCus,1, nIdxVeh); idxVehArray = repelem(idxVeh,1, nIdxCus);
InsDataArray = repmat(thisInsData, [nTotal,1]);

for iSort = 1: nPop
    if iSort <= nTotal && 1 % TODO ���ӳ�ʼ����Ⱥ�Ĳ��� 
        thisSortedInsData = sortThisInsData(InsDataArray(iSort), idxCusArray(:,iSort), idxVehArray(:,iSort));
    elseif iSort <= nPop
%         thisSortedInsData = sortThisInsData(thisInsData, randperm(nCustomer)', idxVeh );%�����̶�˳�� TODO ��ѡһ
        thisSortedInsData = sortThisInsData(thisInsData,randperm(nCustomer)', randperm(nVehicle)' );%�������˳��
    else
        error('EEEEE');
    end
    PopArray(iSort) =  thisSortedInsData;
end

end

function PopSolutionArray = computeFitness(PopArray)
%fitness BF�ȼ�����Ⱥÿ��Ⱦɫ���Fitness
% ����
%   PopArray             ��Ⱥ(����Fitness��Ŀ��ֵ)
% ���
%   PopSolutionArray     ��Ⱥ(����Fitness��Ŀ��ֵ)

PopSolutionArray = [];
for pIter = 1 : length(PopArray)
    PopSolutionArray = [PopSolutionArray callHeuBF(PopArray(pIter))];  %BF�㷨����Fitness
    PopSolutionArray(pIter).algName = 'Genetic Algorithm';
    PopSolutionArray(pIter).algArg = 'TBA';
end
end % End of fitness