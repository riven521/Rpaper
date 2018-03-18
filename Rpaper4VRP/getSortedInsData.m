function SortedInsDataArray = getSortedInsData(thisInsData)
%getSortedInsData 获取数据文件不同排序后的更新数据
% 输入
%   thisInsData                     待排序的instance数据
% 输出
%   SortedInsDataArray     排序后的instance数据
%% function getSortedInsData(thisInsData)

idxCus = {};   idxVeh = {};

% NOTE 核心 修改customer和vehicle顺序
% idx : Customer的排序顺序,即 1:n 的一个顺序
idx = (1:length(thisInsData.Customer.Demand))';              idxCus{end+1} = idx;     %获取正常的顺序
[~, idx] = sort(thisInsData.Customer.Demand,'descend');      idxCus{end+1} = idx;  %获取递减的V_size  最重要
[~, idx] = sort(thisInsData.Customer.Demand,'ascend');         idxCus{end+1} = idx; %获取递增的V_size 最重要

idx = (1:length(thisInsData.Vehicle.Capacity))';               idxVeh{end+1} = idx;     %获取正常的顺序
% [~,idx] = sort(thisInsData.Vehicle.Capacity,'ascend');         idxVeh{end+1} = idx;
% [~,idx] = sort(thisInsData.Vehicle.Capacity,'ascend');            idxVeh{end+1} = idx;
% [~,idx] = sort(thisInsData.Vehicle.UnitCost,'ascend');            idxVeh{end+1} = idx;

% 去除重复的customer和vehicle索引顺序
idxCus = unique(cell2mat(idxCus)','rows')';
idxVeh = unique(cell2mat(idxVeh)','rows')';

% 重复一定次数 方便arrayfun使用
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
    % 等同上面for循环
    % SortedInsDataArray = arrayfun(@sortThisInsData, ...
    %     InsDataArray, idxCusArray, idxVehArray);
    %     printstruct(SortedInsDataArray);
end
