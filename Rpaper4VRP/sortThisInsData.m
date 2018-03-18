function ThisInsDataSorted = sortThisInsData(ThisInsData, thisidxCus, thisidxVeh)
%sortThisInsData 对数据文件按customer和vehicle顺序排序   thisidxCus
% 输入
%   ThisInsData          待排序的instance数据
%   thisidxCus           给定customer顺序
%   thisidxVeh           给定Vehicle顺序
% 输出
%   ThisInsDataSorted    排序后的instance数据

% NOTE 轻易不能调整insData数据结构; 调整后必须修改该函数 

if size(thisidxCus,2)>1 || size(thisidxVeh,2)>1
%     thisidxCus
%     thisidxVeh
    error('EEEEE');  
end

cus = (thisidxCus);
veh = (thisidxVeh);
da = ThisInsData;

%算例d依据顾客顺序idxCus和车辆顺序idxVeh变化 TODO 后续可能要增加Cluster顺序变化
da.Customer.Idx = thisidxCus;   %给定顾客顺序
da.Customer.Demand = da.Customer.Demand(cus,:);
da.Customer.IdxCluster = da.Customer.IdxCluster(cus,:);
da.Customer.Compatible = da.Customer.Compatible(cus,cus);
da.Customer.Coord = da.Customer.Coord(cus,:);
da.Customer.Distance = da.Customer.Distance(cus,cus);

da.Depot.CustomerDistance = da.Depot.CustomerDistance(cus,:);

da.Node.Demand = [da.Depot.Demand; da.Customer.Demand];
da.Node.Coord = [da.Depot.Coord; da.Customer.Coord];
da.Node.Distance(1,:) = [ 0; da.Depot.CustomerDistance]';
da.Node.Distance(:,1) = [ 0; da.Depot.CustomerDistance];
da.Node.Distance(2:end,2:end) = da.Customer.Distance;

da.Vehicle.Idx = thisidxVeh;    %给定车辆顺序
da.Vehicle.Capacity = da.Vehicle.Capacity(veh,:);
da.Vehicle.FixCost = da.Vehicle.FixCost(veh,:);
da.Vehicle.MaxPoint = da.Vehicle.MaxPoint(veh,:);
da.Vehicle.PointCost = da.Vehicle.PointCost(veh,:);
da.Vehicle.VariableCost = da.Vehicle.VariableCost(veh,:);
da.Vehicle.UnitCost = da.Vehicle.UnitCost(veh,:);

da.Veh_Cus.Cost_Num = da.Veh_Cus.Cost_Num(veh,cus);
%     insDataNew.Veh_Cus.Cost_Type= unique(insDataNew.Veh_Cus.Cost_Num,
%     'rows'); %如需要按类型的车辆数据结构 再启用该语句

ThisInsDataSorted = da;
end
