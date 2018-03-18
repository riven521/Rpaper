function printSolu(thisSolution)
%printSolu 输出解 包含每条路径信息和总路径信息
% 输入
%   thisSolution    打印完整的解(在Vehicle.Cost_***中)

nVeh = length(thisSolution.instance.Vehicle.Capacity);
wff = thisSolution.instance.Veh_Cus.wff;
if ~isCompletedWff(wff), warning('警告: 解不完整'); end

usedVehIdx = find(sum(wff,2)>0);

r = {}; %每条路径信息
r{end+1} = usedVehIdx;
r{end+1} = thisSolution.instance.Vehicle.Capacity(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_nCustomer(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_Total(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_Route(usedVehIdx);

t = {}; %所有路径信息
t{end+1} = nVeh;
t{end+1} = numel(usedVehIdx);
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx))/...
    sum(thisSolution.instance.Vehicle.Capacity(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Capacity(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_nCustomer(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_Total(usedVehIdx));

% fprintf('采用算法 %s \n', thisSolution.algName); fprintf('算法参数 %s \n', thisSolution.algArg);

for x=1: numel(usedVehIdx)
    fprintf('Vehicle %d (Capacity =%d; Load = %d; Point = %d; Cost = %d )  visit Customer %s \n',...
        r{1}(x),r{2}(x),r{3}(x),r{4}(x),r{5}(x),mat2str(r{6}{x}));
end

fprintf('共%d俩车,使用%d辆;总装载率=%d; 总容量 = %d; 总装载 = %d; 总数量 = %d; 总费用 = %d  \n', t{1},t{2},t{3},t{4},t{5},t{6},t{7});

end

