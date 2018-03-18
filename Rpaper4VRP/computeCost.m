function ThisInsDataWithCost = computeCost(ThisInsData)
%computeCost  计算ThisInsData中每条路径的费用(必须有更新的wff)
% 输入
%   ThisInsData                  完整含wff的数据结构体 
% 输出
%   ThisInsDataWithCost    增加含Vehilce.Cost_***的数据结构体
%% function ThisInsDataWithCost = computeCost(ThisInsData)

da = ThisInsData;
wff = da.Veh_Cus.wff;

if ~isCompletedWff(wff),  warning('Wff不是完整的解，正在计算非完整解的成本'); end

% 分别计算 1 基价 2 增点费 3 总费用 4 使用容量（剩余容量） 5 总访问数 6 访问点  （均和wff相关）
da.Vehicle.Cost_Base = max(wff .* da.Veh_Cus.Cost_Num,[],2);

da.Vehicle.Cost_Point = da.Vehicle.PointCost .* ( sum(wff,2) - 1);
da.Vehicle.Cost_Point(da.Vehicle.Cost_Point < 0)=0;

da.Vehicle.Cost_Total = da.Vehicle.Cost_Base + da.Vehicle.Cost_Point;

da.Vehicle.Cost_TotalDemand = wff * da.Customer.Demand;
da.Vehicle.Cost_ResidueCapacity = da.Vehicle.Capacity - da.Vehicle.Cost_TotalDemand;

da.Vehicle.Cost_nCustomer = sum(wff,2);

da.Vehicle.Cost_Route = ...
    arrayfun(@(x) find(wff(x,:)==1),(1:length(da.Vehicle.Capacity)),...
    'UniformOutput',false);

if any(da.Vehicle.Cost_Base<0) || any(da.Vehicle.Cost_Point <0) || any(da.Vehicle.Cost_Total)<0 ...
        || any(da.Vehicle.Cost_TotalDemand)<0 || any(da.Vehicle.Cost_ResidueCapacity)<0 || any(da.Vehicle.Cost_nCustomer)<0
    fprintf("存在路径计算错误");
end

ThisInsDataWithCost = da;
end