function ThisSolution = callHeuBF(ThisInsData) 
%callHeuBF  ����BF�㷨   TODO  ����Ԥ������Ӱ��
% ����
%   ThisSortedInsData     �㷨�������ݼ�
% ���
%   ThisSolution              ����*��*�����ݼ��ṹ��
%% function callHeuBF(ThisInsData) 
% Local Functions:
% # getAvailableVehicles
% # findBestFitVehicle    TODO ��д ThisSortedInsData

ThisSolution.algName = 'Best Fit Algorithm';
ThisSolution.algArg = 'TBA';

da = ThisInsData;
nCustomer = length(da.Customer.Demand);
nVehicle = length(da.Vehicle.Capacity);
iwff = zeros(nVehicle, nCustomer);  %wff ����ĳ��(��Ӧ��)����ĳ��(���ó���)

% ѭ��Ϊÿ��iCus�ҳ������vehicle����������iCus
for iCus = 1 : nCustomer
    % 1 �ҳ����Է���customer i�ĳ�������
    [idxVehNewArray,idxVehUsedArray,idxVehArray] = getVehicleIdx(iCus, iwff, da);    
    % 2 �ҳ����Է���customer i������ʵ��Ǹ�����
    thisBsetVehicle = findBestFitVehicle(iCus, iwff, da, ...
        idxVehNewArray,idxVehUsedArray,idxVehArray);                    
    % 3 ��iCus���복�����ΪthisBsetVehicle�ĳ�,���¶�Ӧiwff����
    iwff(thisBsetVehicle,iCus) = 1;    
    % 4 �������
%     fprintf('Node %d (%d) Veh %d (%d; %d)   \n ', iCus, da.Customer.Demand(iCus), thisBsetVehicle, da.Vehicle.Capacity(thisBsetVehicle) , da.Vehicle.MaxPoint(thisBsetVehicle))
end %����BF�㷨


if ~isCompletedWff(iwff),  error('�������й˿Ͷ����복����'); end

% ����Cost������
%   WFF��ȡ�ɹ�������Cost����Ҫwff����
da.Veh_Cus.wff = iwff;
ThisSolution.instance = computeCost(da); %NOTE ���㷨����������ҪCostʱ��������;BF�ɲ����� FIXME ������ѡ�ȽϺ�
% printstruct(ThisSolution.instance.Vehicle)

end


%%  % �ֲ� FUNCTIONS %
% ==========================================================================
% �ֲ� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Local function 555 getAvailableVehicles
%
%
function [idxAvailableNewArray,idxAvailableUsedArray,idxAvailableArray] = getVehicleIdx(iCus, iwff, ThisSortedInsData)
%getAvailableVehicle �ҳ���������Լ���ĳ������� ThisSortedInsData
%

da = ThisSortedInsData;

% Լ��1: ��������Լ����Logical array
AssignedCapa = iwff * da.Customer.Demand;
ResidueCapa = da.Vehicle.Capacity - AssignedCapa;

idxCapacity = ( ResidueCapa >= da.Customer.Demand(iCus) );
idxNewCapacity = (ResidueCapa >= da.Customer.Demand(iCus)) ... 
                            & (ResidueCapa == da.Vehicle.Capacity );
idxUsedCapacity = (ResidueCapa >= da.Customer.Demand(iCus)) ...
                            & ( ResidueCapa < da.Vehicle.Capacity );
if isempty(idxCapacity), error('û�г�������(����)Ҫ�������ӿ��ó�����'); end

% Լ��2: ������������Logical array
AssignedPoint = sum(iwff,2);

idxPoint = ( (da.Vehicle.MaxPoint - AssignedPoint ) > 0 );
if isempty(idxPoint), error('û�г�������(������)Ҫ�����������ɷ��ʵ���'); end

% Լ��3: ����˿�i�복���������˿���ͨ�Ե�Logical array  �Ƚ������
%   ��ȡ����iCus������customer����ͨ����, ���������ɳ����� (��������������)
nVeh = length(da.Vehicle.Capacity);
customerIMat = repmat( da.Customer.Compatible(iCus,:), [nVeh 1]); 

%   iwff��customerIMat����0-1�������; ����Ӻ���Ϊ0, ��������ͨ��Ҫ��
idxCompatible = ( sum(iwff .* customerIMat, 2) == 0 );
if isempty(idxCompatible), error('û�г�������(�˿���ͨ��)Ҫ�������ӹ˿���ͨ��'); end

% ����: ������������Լ���ĳ���Logical array����
idxAvailableArray = all( [idxCapacity idxPoint idxCompatible], 2);
idxAvailableNewArray = all( [idxNewCapacity idxPoint idxCompatible], 2);
idxAvailableUsedArray = all( [idxUsedCapacity idxPoint idxCompatible], 2);

if (idxAvailableArray ~= any([idxAvailableNewArray idxAvailableUsedArray],2)),  error('�¾ɳ�֮�Ͳ������ܿ��ó�����'); end
if (sum(idxAvailableArray) == 0),  error('�޳�����,�����ӳ������򽵵���ͨ�Ի������'); end

end

%% 2 Local function 555 findBestFitVehicle
%
%
function thisIdxValue = findBestFitVehicle(iCus, iwff, ThisSortedInsData, idxVehNewArray,idxVehUsedArray,idxVehArray)
%��UsedVehicle��BestFit��Idx�������NewVehicl�� ThisSortedInsData
da = ThisSortedInsData;

% ��ʼ������ $ vehicleAbsCost = |Vehicle_{basecost} - iCus_{basecost} |  $
vehicleBaseCost = max( iwff .* da.Veh_Cus.Cost_Num, [], 2) ; %��ǰVehicle�л���(��ǰ�������е������ֵ)
vehicleAbsCost = abs( da.Veh_Cus.Cost_Num(:,iCus) - vehicleBaseCost ); %��ǰ��������۵ľ���ֵ

%���Used��������,������һ����
if any(idxVehUsedArray)
    thisMinAbsCost = min( vehicleAbsCost ( idxVehUsedArray )); %�ҵ�used Vehicle����Сֵ
    % �ҵ� 1 ������Сֵ 2 �̺���idxAvailableUsed �е�IDXֵ
    idxArray = (vehicleAbsCost == thisMinAbsCost);
    idxArray = logical(idxArray .* idxVehUsedArray);
elseif any(idxVehNewArray)
    thisMinAbsCost = min(vehicleAbsCost(idxVehNewArray));
    % �ҵ� 1 ������Сֵ 2 �̺���idxAvailableNew �е�IDXֵ
    idxArray = (vehicleAbsCost == thisMinAbsCost);
    idxArray = logical(idxArray .* idxVehNewArray);
else
    error('û���κγ������ã������ӿ��ó�����');
end
%���³���ɳ����Ҷ�����
% %         minAbsCost = min(AbsCost(idxAvailable));

idxValueArray = find(idxArray);  %�ҵ�logical��Ӧ��idxλ��
if isempty(idxValueArray)
    error('û���κγ������ã�������bug'); 
else
    thisIdxValue = idxValueArray(1);  % NOTE ���ɶ����Ϊ��Сֵ�ĳ���ʱ ѡȡ��һ��
%     thisIdx = thisIdx(randsample(numel(thisIdxValue),1)); ���ɶ����Ϊ��Сֵ�ĳ���ʱ ѡȡ���һ��
end

if numel(thisIdxValue) ~=1, error('ERROR'); end

end
