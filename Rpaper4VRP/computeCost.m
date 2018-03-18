function ThisInsDataWithCost = computeCost(ThisInsData)
%computeCost  ����ThisInsData��ÿ��·���ķ���(�����и��µ�wff)
% ����
%   ThisInsData                  ������wff�����ݽṹ�� 
% ���
%   ThisInsDataWithCost    ���Ӻ�Vehilce.Cost_***�����ݽṹ��
%% function ThisInsDataWithCost = computeCost(ThisInsData)

da = ThisInsData;
wff = da.Veh_Cus.wff;

if ~isCompletedWff(wff),  warning('Wff���������Ľ⣬���ڼ����������ĳɱ�'); end

% �ֱ���� 1 ���� 2 ����� 3 �ܷ��� 4 ʹ��������ʣ�������� 5 �ܷ����� 6 ���ʵ�  ������wff��أ�
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
    fprintf("����·���������");
end

ThisInsDataWithCost = da;
end