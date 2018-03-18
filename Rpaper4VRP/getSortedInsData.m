function SortedInsDataArray = getSortedInsData(thisInsData)
%getSortedInsData ��ȡ�����ļ���ͬ�����ĸ�������
% ����
%   thisInsData                     �������instance����
% ���
%   SortedInsDataArray     ������instance����
%% function getSortedInsData(thisInsData)

idxCus = {};   idxVeh = {};

% NOTE ���� �޸�customer��vehicle˳��
% idx : Customer������˳��,�� 1:n ��һ��˳��
idx = (1:length(thisInsData.Customer.Demand))';              idxCus{end+1} = idx;     %��ȡ������˳��
[~, idx] = sort(thisInsData.Customer.Demand,'descend');      idxCus{end+1} = idx;  %��ȡ�ݼ���V_size  ����Ҫ
[~, idx] = sort(thisInsData.Customer.Demand,'ascend');         idxCus{end+1} = idx; %��ȡ������V_size ����Ҫ

idx = (1:length(thisInsData.Vehicle.Capacity))';               idxVeh{end+1} = idx;     %��ȡ������˳��
% [~,idx] = sort(thisInsData.Vehicle.Capacity,'ascend');         idxVeh{end+1} = idx;
% [~,idx] = sort(thisInsData.Vehicle.Capacity,'ascend');            idxVeh{end+1} = idx;
% [~,idx] = sort(thisInsData.Vehicle.UnitCost,'ascend');            idxVeh{end+1} = idx;

% ȥ���ظ���customer��vehicle����˳��
idxCus = unique(cell2mat(idxCus)','rows')';
idxVeh = unique(cell2mat(idxVeh)','rows')';

% �ظ�һ������ ����arrayfunʹ��
nIdxCus = size(idxCus,2);
nIdxVeh = size(idxVeh,2);
nTotal = nIdxCus * nIdxVeh;
idxCusArray = repelem(idxCus,1, nIdxVeh);   %     celldisp(idxCus1);
idxVehArray = repelem(idxVeh,1, nIdxCus);
InsDataArray = repmat(thisInsData, [nTotal,1]);

for iSort = 1:nTotal
    thisSortedInsData = sortThisInsData(InsDataArray(iSort), idxCusArray(:,iSort), idxVehArray(:,iSort));
    SortedInsDataArray(iSort) =  thisSortedInsData;
end
    % ��ͬ����forѭ��
    % SortedInsDataArray = arrayfun(@sortThisInsData, ...
    %     InsDataArray, idxCusArray, idxVehArray);
    %     printstruct(SortedInsDataArray);
end
