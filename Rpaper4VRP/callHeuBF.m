function ThisSolution = callHeuBF(ThisInsData) 
%callHeuBF  调用BF算法   TODO  增加预处理后的影响
% 输入
%   ThisSortedInsData     算法计算数据集
% 输出
%   ThisSolution              包含*解*的数据集结构体
%% function callHeuBF(ThisInsData) 
% Local Functions:
% # getAvailableVehicles
% # findBestFitVehicle    TODO 待写 ThisSortedInsData

ThisSolution.algName = 'Best Fit Algorithm';
ThisSolution.algArg = 'TBA';

da = ThisInsData;
nCustomer = length(da.Customer.Demand);
nVehicle = length(da.Vehicle.Capacity);
iwff = zeros(nVehicle, nCustomer);  %wff 代表某列(供应商)放入某行(可用车数)

% 循环为每个iCus找出最合适vehicle索引并放入iCus
for iCus = 1 : nCustomer
    % 1 找出可以放入customer i的车辆集合
    [idxVehNewArray,idxVehUsedArray,idxVehArray] = getVehicleIdx(iCus, iwff, da);    
    % 2 找出可以放入customer i的最合适的那个车辆
    thisBsetVehicle = findBestFitVehicle(iCus, iwff, da, ...
        idxVehNewArray,idxVehUsedArray,idxVehArray);                    
    % 3 将iCus放入车辆序号为thisBsetVehicle的车,更新对应iwff数据
    iwff(thisBsetVehicle,iCus) = 1;    
    % 4 输出过程
%     fprintf('Node %d (%d) Veh %d (%d; %d)   \n ', iCus, da.Customer.Demand(iCus), thisBsetVehicle, da.Vehicle.Capacity(thisBsetVehicle) , da.Vehicle.MaxPoint(thisBsetVehicle))
end %结束BF算法


if ~isCompletedWff(iwff),  error('不是所有顾客都放入车辆中'); end

% 更新Cost并返回
%   WFF获取成功：计算Cost：需要wff数据
da.Veh_Cus.wff = iwff;
ThisSolution.instance = computeCost(da); %NOTE 当算法继续运行需要Cost时必须运算;BF可不运算 FIXME 参数可选比较好
% printstruct(ThisSolution.instance.Vehicle)

end


%%  % 局部 FUNCTIONS %
% ==========================================================================
% 局部 FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Local function 555 getAvailableVehicles
%
%
function [idxAvailableNewArray,idxAvailableUsedArray,idxAvailableArray] = getVehicleIdx(iCus, iwff, ThisSortedInsData)
%getAvailableVehicle 找出满足所有约束的车辆集合 ThisSortedInsData
%

da = ThisSortedInsData;

% 约束1: 满足容量约束的Logical array
AssignedCapa = iwff * da.Customer.Demand;
ResidueCapa = da.Vehicle.Capacity - AssignedCapa;

idxCapacity = ( ResidueCapa >= da.Customer.Demand(iCus) );
idxNewCapacity = (ResidueCapa >= da.Customer.Demand(iCus)) ... 
                            & (ResidueCapa == da.Vehicle.Capacity );
idxUsedCapacity = (ResidueCapa >= da.Customer.Demand(iCus)) ...
                            & ( ResidueCapa < da.Vehicle.Capacity );
if isempty(idxCapacity), error('没有车辆满足(容量)要求，请增加可用车数量'); end

% 约束2: 满足最大点数的Logical array
AssignedPoint = sum(iwff,2);

idxPoint = ( (da.Vehicle.MaxPoint - AssignedPoint ) > 0 );
if isempty(idxPoint), error('没有车辆满足(最大点数)要求，请增加最大可访问点数'); end

% 约束3: 满足顾客i与车辆内其他顾客连通性的Logical array  比较难理解
%   获取带你iCus与其它customer的连通性行, 并复制若干车辆个 (方便下面矩阵相乘)
nVeh = length(da.Vehicle.Capacity);
customerIMat = repmat( da.Customer.Compatible(iCus,:), [nVeh 1]); 

%   iwff和customerIMat两个0-1矩阵相乘; 行相加后仍为0, 即满足连通性要求
idxCompatible = ( sum(iwff .* customerIMat, 2) == 0 );
if isempty(idxCompatible), error('没有车辆满足(顾客连通性)要求，请增加顾客连通性'); end

% 返回: 满足所有以上约束的车辆Logical array返回
idxAvailableArray = all( [idxCapacity idxPoint idxCompatible], 2);
idxAvailableNewArray = all( [idxNewCapacity idxPoint idxCompatible], 2);
idxAvailableUsedArray = all( [idxUsedCapacity idxPoint idxCompatible], 2);

if (idxAvailableArray ~= any([idxAvailableNewArray idxAvailableUsedArray],2)),  error('新旧车之和不等于总可用车辆数'); end
if (sum(idxAvailableArray) == 0),  error('无车可用,请增加车辆数或降低连通性或聚类数'); end

end

%% 2 Local function 555 findBestFitVehicle
%
%
function thisIdxValue = findBestFitVehicle(iCus, iwff, ThisSortedInsData, idxVehNewArray,idxVehUsedArray,idxVehArray)
%从UsedVehicle找BestFit的Idx；否则从NewVehicl找 ThisSortedInsData
da = ThisSortedInsData;

% 初始化计算 $ vehicleAbsCost = |Vehicle_{basecost} - iCus_{basecost} |  $
vehicleBaseCost = max( iwff .* da.Veh_Cus.Cost_Num, [], 2) ; %当前Vehicle中基价(当前车内所有点中最大值)
vehicleAbsCost = abs( da.Veh_Cus.Cost_Num(:,iCus) - vehicleBaseCost ); %当前车辆与基价的绝对值

%如果Used车辆存在,必须找一个放
if any(idxVehUsedArray)
    thisMinAbsCost = min( vehicleAbsCost ( idxVehUsedArray )); %找到used Vehicle中最小值
    % 找到 1 等于最小值 2 蕴含在idxAvailableUsed 中的IDX值
    idxArray = (vehicleAbsCost == thisMinAbsCost);
    idxArray = logical(idxArray .* idxVehUsedArray);
elseif any(idxVehNewArray)
    thisMinAbsCost = min(vehicleAbsCost(idxVehNewArray));
    % 找到 1 等于最小值 2 蕴含在idxAvailableNew 中的IDX值
    idxArray = (vehicleAbsCost == thisMinAbsCost);
    idxArray = logical(idxArray .* idxVehNewArray);
else
    error('没有任何车辆可用，请增加可用车数量');
end
%从新车或旧车中找都可以
% %         minAbsCost = min(AbsCost(idxAvailable));

idxValueArray = find(idxArray);  %找到logical对应的idx位置
if isempty(idxValueArray)
    error('没有任何车辆可用，代码有bug'); 
else
    thisIdxValue = idxValueArray(1);  % NOTE 当由多个均为最小值的车辆时 选取第一个
%     thisIdx = thisIdx(randsample(numel(thisIdxValue),1)); 当由多个均为最小值的车辆时 选取随机一个
end

if numel(thisIdxValue) ~=1, error('ERROR'); end

end
