function printSolu(thisSolution)
%printSolu ����� ����ÿ��·����Ϣ����·����Ϣ
% ����
%   thisSolution    ��ӡ�����Ľ�(��Vehicle.Cost_***��)

nVeh = length(thisSolution.instance.Vehicle.Capacity);
wff = thisSolution.instance.Veh_Cus.wff;
if ~isCompletedWff(wff), warning('����: �ⲻ����'); end

usedVehIdx = find(sum(wff,2)>0);

r = {}; %ÿ��·����Ϣ
r{end+1} = usedVehIdx;
r{end+1} = thisSolution.instance.Vehicle.Capacity(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_nCustomer(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_Total(usedVehIdx);
r{end+1} = thisSolution.instance.Vehicle.Cost_Route(usedVehIdx);

t = {}; %����·����Ϣ
t{end+1} = nVeh;
t{end+1} = numel(usedVehIdx);
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx))/...
    sum(thisSolution.instance.Vehicle.Capacity(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Capacity(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_TotalDemand(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_nCustomer(usedVehIdx));
t{end+1} = sum(thisSolution.instance.Vehicle.Cost_Total(usedVehIdx));

% fprintf('�����㷨 %s \n', thisSolution.algName); fprintf('�㷨���� %s \n', thisSolution.algArg);

for x=1: numel(usedVehIdx)
    fprintf('Vehicle %d (Capacity =%d; Load = %d; Point = %d; Cost = %d )  visit Customer %s \n',...
        r{1}(x),r{2}(x),r{3}(x),r{4}(x),r{5}(x),mat2str(r{6}{x}));
end

fprintf('��%d����,ʹ��%d��;��װ����=%d; ������ = %d; ��װ�� = %d; ������ = %d; �ܷ��� = %d  \n', t{1},t{2},t{3},t{4},t{5},t{6},t{7});

end

