function ThisInsDataSorted = sortThisInsData(ThisInsData, thisidxCus, thisidxVeh)
%sortThisInsData �������ļ���customer��vehicle˳������   thisidxCus
% ����
%   ThisInsData          �������instance����
%   thisidxCus           ����customer˳��
%   thisidxVeh           ����Vehicle˳��
% ���
%   ThisInsDataSorted    ������instance����

% NOTE ���ײ��ܵ���insData���ݽṹ; ����������޸ĸú��� 

if size(thisidxCus,2)>1 || size(thisidxVeh,2)>1
%     thisidxCus
%     thisidxVeh
    error('EEEEE');  
end

cus = (thisidxCus);
veh = (thisidxVeh);
da = ThisInsData;

%����d���ݹ˿�˳��idxCus�ͳ���˳��idxVeh�仯 TODO ��������Ҫ����Cluster˳��仯
da.Customer.Idx = thisidxCus;   %�����˿�˳��
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

da.Vehicle.Idx = thisidxVeh;    %��������˳��
da.Vehicle.Capacity = da.Vehicle.Capacity(veh,:);
da.Vehicle.FixCost = da.Vehicle.FixCost(veh,:);
da.Vehicle.MaxPoint = da.Vehicle.MaxPoint(veh,:);
da.Vehicle.PointCost = da.Vehicle.PointCost(veh,:);
da.Vehicle.VariableCost = da.Vehicle.VariableCost(veh,:);
da.Vehicle.UnitCost = da.Vehicle.UnitCost(veh,:);

da.Veh_Cus.Cost_Num = da.Veh_Cus.Cost_Num(veh,cus);
%     insDataNew.Veh_Cus.Cost_Type= unique(insDataNew.Veh_Cus.Cost_Num,
%     'rows'); %����Ҫ�����͵ĳ������ݽṹ �����ø����

ThisInsDataSorted = da;
end
