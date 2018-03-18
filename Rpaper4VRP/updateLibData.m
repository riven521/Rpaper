function InsDataArray = updateLibData(LibDataArray, ArgIns, fileType)
%updateLibDatas  ����LibDataArrayΪInsDataArray
% ����
%   LibDataArray     lib�ļ���ȡ�����ݽṹ�� ��Ҫ����
%   ArgIns                instance��ز����ṹ�� PHO GAMMA
%   fileType             fileType�ļ�����
% ���
%   InsDataArray    ���ݼ��ṹ��
%% function updateLibData(LibDataArray, ArgIns, fileType)
% Local Functions:
% # VRP_Update
% # CCVRP_Update    TODO ��д
% # checkInsDataArray

%% 
%{ InsDataArray.Node = struct('Coord',[],'Demand',[],'Distance',[]);
% InsDataArray.Customer = struct('Coord',[],'Demand',[],'IdxCluster',[],'Distance',[],'Compatible',[]);
% InsDataArray.Depot = struct('Coord',[],'Demand',[],'CustomerDistance',[],'ClusterDistance',[]);
% InsDataArray.Cluster = struct('Coord',[],'Demand',[],'Distance',[],'Compatible',[]);
% InsDataArray.Vehicle = struct('Capacity',[],'FixCost',[],'VariableCost',[],...
%     'PointCost',[],'MaxPoint',[],'Route',[],'Cost',[]); %Route����wff
% InsDataArray.Veh_Cus = struct('Cost_Num',[],'Cost_Type',[],'wff',[]);
%} InsDataArray = repmat(InsDataArray,[length(LibDataArray) 1]);
%%
% ����LibDataArrayΪInsDataArray
ArgInsArray  = repmat(ArgIns,[length(LibDataArray) 1]);
if strcmp('vrp', fileType)
    InsDataArray = arrayfun(@VRP_Update,LibDataArray,ArgInsArray);   %     InsDataArray = arrayfun(@VRP_Update,LibDataArray,InsDataArray,ArgInsArray); 
elseif strcmp('ccvrp', fileType)
    InsDataArray = arrayfun(@CCVRP_Update,LibDataArray,ArgInsArray);
end

% �ж�InsDataArray�Ƿ����
if all(arrayfun(@checkInsDataArray,InsDataArray)),  fprintf("All data is corrected ! \n ");  end

fprintf('Updating file with updateLibDatas() ... done \n');

end

%%  % �ֲ� FUNCTIONS %
% ==========================================================================
% �ֲ� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Local function 555 VRP_Update
%
%
function InsDataArray = VRP_Update(ThisLibData,ThisArgIns)
% NOTE ��Ҫ����: ����VRP's Lib����

%%
% Get Node(Demand;Coord); Cluster(Demand;Coord); Customer(Coor; Demand);
%% Number
NodeNumber = length(ThisLibData.Node_Coord);
DepotNumber = 1; % FIXME:��Depot��ֻһ������λ��0/1ʱҪ�޸�;��DEPOT_SECTION����
CustomerNumber = NodeNumber - DepotNumber;
% VehicleNumber = sum(libData.VEHICLE_Number); %libData.VEHICLE_Number��ÿ�����ͳ���������
%% Coord
InsDataArray.Node.Coord = ThisLibData.Node_Coord;
InsDataArray.Depot.Coord = ThisLibData.Node_Coord(DepotNumber,:);
InsDataArray.Customer.Coord = ThisLibData.Node_Coord(DepotNumber+1:end,:);
%% Demand
InsDataArray.Node.Demand = ThisLibData.Node_Demand;
InsDataArray.Depot.Demand = InsDataArray.Node.Demand(DepotNumber,:);
InsDataArray.Customer.Demand = InsDataArray.Node.Demand(DepotNumber+1:end,:);

%% Cluster
% ʹ�ò���PHO
[InsDataArray.Cluster.Coord, InsDataArray.Cluster.Demand, InsDataArray.Customer.IdxCluster ] = getCluster();
InsDataArray.Depot.IdxCluster = max(InsDataArray.Customer.IdxCluster)+1;  %FIXME ����can����depot�ľ���Ϊ��0��/������Ϊ��MAX+1��

    function [cluCoord,cluDemand,idx] = getCluster()
        n = getClusterN();
        [idx, cluCoord,~,~] = kmeans(InsDataArray.Customer.Coord, n); % 2 NOTE K-means�㷨��ȡ��������Coord
        cluDemand = arrayfun(@(x) sum(InsDataArray.Customer.Demand (idx == x)), (1:n)' ); % 5 ����ÿ������������� Demand;
        
        % ����ÿ�����������������ʱ���ã�
        % NumberOneCluster = arrayfun(@(x) sum( D.Customer.IdxCluster == x ),(1:n)');
        
        %% Ƕ�� function getClusterN()
        %
        %
        function ClusterNumber = getClusterN()
            % NOTE ˵�� ���������ɼ���{
            % ���ȷ�������� - ����2016CIE����PHO*Q;
            % ���ȷ�������������;
            %} case1: ��������<=���������; case2: ��������ȫ��>��������� case3: �����������ڻ�С�����������
            
            minClusterNumber = 1;
            maxClusterNumber = NodeNumber - 1;
            
            PHO = ThisArgIns.PHO;
            if (PHO > 0 && PHO < 1)
                %   ClusterNumber = ceil(PHO*sum(TSP_Demand(2:end)')/TSP_Capacity);
                ClusterNumber = ceil(PHO * maxClusterNumber);   % NOTE ClusterNumber ������ ���Բ��� para.PHO [0,1]) ���������С1��,��಻������ĸ���(��������)
            elseif PHO == 0 %���PHO=0, ֻ��һ������; PHO��Ϊ0Ҳ����ֻ��һ������
                ClusterNumber = minClusterNumber;
            elseif PHO == 1 %isinf(PHO)
                ClusterNumber = maxClusterNumber;
            else
                error('PHO = %d, ����[0,1]֮��', PHO);
            end
            
            if (ClusterNumber < minClusterNumber || ClusterNumber > maxClusterNumber)
                warning('ClusterNumber = %d, ������Χ', ClusterNumber);
            else
                %                 fprintf('ClusterNumber = %d is in [1, %d] \n',ClusterNumber,maxClusterNumber);
            end
        end
        
    end

%% Compatible
% ʹ�ò���GAMMA
[InsDataArray.Cluster.Compatible,InsDataArray.Customer.Compatible] = getCompatible();

%% Ƕ�� function getCompatible()
%
%
    function [M_Cludepot,M_depot] = getCompatible()
        % NOTE ˵�� ��ͨ�����ɼ���{ % GAMMA Խ��, 0�ĸ���Խ��, ��ͨ��Խ��
        % GET M_Cludepot  ������ͨ�� (0:����ͨ; 1:������ͨ)
        n = length(InsDataArray.Cluster.Demand);
        if n==1
            M_Cludepot = 0;
            M_depot = zeros(CustomerNumber,CustomerNumber);
        else
            GAMMA = ThisArgIns.GAMMA;  % ClusterNumber=5 % GAMMA = 0.5;
            % ����totalZeros :һ����Ҫ��0�ĸ���,���ڷ�Χ:[minZeros,maxZeros]
            % GAMMA = 0.5; Ӧǡ��ֻ��һ����ͨ��
            minZeros = n; %�Խ�����0������
            maxZeros = n * n;
            totalZeros = minZeros + GAMMA * (maxZeros - minZeros);
            
            numZeros = minZeros;
            M_Cludepot = ones(n,n);
            M_Cludepot = M_Cludepot-diag(diag(M_Cludepot));
            
            % NOTE ����:������ͨ�Ծ��� M_Cludepot
            % ��ʼ���ͳ��Խ������Ϊ1 ѭ������0�ĸ���until����totalZeros
            while 1
                x = randi([1 n],1);
                y = randi([1 n],1);
                if numZeros >= totalZeros,  break;  end
                if x ~= y && M_Cludepot(x,y)~=0
                    M_Cludepot(x,y) = 0;
                    M_Cludepot(y,x) = 0;
                    numZeros = numZeros + 2;
                else
                    %   fprintf('x = %d = y = %d or M_Cludepot(x,y) = %d \n',x,y,M_Cludepot(x,y))
                end
            end
            
            
            % GET M_depot �����ͨ�� ���ھ�����ͨ�Ե���չ
            M_depot = zeros(CustomerNumber,CustomerNumber);
            for tempi = 1 : CustomerNumber - 1
                for tempj = tempi + 1 : CustomerNumber
                    clusteri = InsDataArray.Customer.IdxCluster(tempi);
                    clusterj = InsDataArray.Customer.IdxCluster(tempj);
                    if M_Cludepot(clusteri,clusterj) == 1
                        M_depot(tempi,tempj) = 1;
                        M_depot(tempj,tempi) = 1;
                    end
                end
            end
            
            % �������������GAMMA,��֤���ݺ����Էǳ���Ҫ
            ifsymmetric1 = (M_Cludepot == M_Cludepot');
            ifsymmetric2 = (M_depot == M_depot');
            if (totalZeros < minZeros || totalZeros > maxZeros || numZeros < totalZeros)
                error('totalZeros = %d, ������Χ', totalZeros);
            elseif (GAMMA <0 || GAMMA >1)
                error('GAMMA = %d, ������Χ [0,1]', GAMMA);
            elseif (~all(ifsymmetric1(:)) && ~all(ifsymmetric2(:)))
                error('M_Cludepot��M_depot���ǶԳƾ���');
            else
                %             fprintf('CORRECT : totalZeros = %d in [%d, %d] and GAMMA = %d, \n',totalZeros,minZeros,maxZeros,GAMMA);
            end
        end
        
    end

%% Distance
% ����ŷʽ����
[InsDataArray.Node.Distance,InsDataArray.Customer.Distance,...
    InsDataArray.Depot.CustomerDistance,InsDataArray.Cluster.Distance,...
    InsDataArray.Depot.ClusterDistance] = getDistance();

%% Ƕ�� function getDistance(libData,insData,arg)
%
%
    function [NodeNodeDist,CusCusdist,DepCusdist,CluCludist,DepCludist] = getDistance()
        
        NodeNodeDist = calEuclidean(InsDataArray.Node.Coord); %����Node����������ŷ�Ͼ���
        CusCusdist = NodeNodeDist(DepotNumber+1:end,DepotNumber+1:end); % Customer֮��ľ���(����Depot)
        DepCusdist = NodeNodeDist(DepotNumber+1:end,1); % Depot����Custoemr��ľ���(����depot����)
        
        ClusterNode_Coord = [InsDataArray.Depot.Coord; InsDataArray.Cluster.Coord];
        DepCluDepCluDist = calEuclidean(ClusterNode_Coord);
        CluCludist = DepCluDepCluDist(DepotNumber+1:end,DepotNumber+1:end);  % Clusters֮��ľ���(����Depot)
        DepCludist = DepCluDepCluDist(DepotNumber+1:end,1);  % Depot��Clusters�ľ���(����Depot)
        % 3 ���Depot��ÿ�����ĵľ���-��Ϊ�̶��ɱ� TSP_ClusterEdgeWeight
        %% function calEuclidean(libData,insData,arg)
        %
        %
        function euclDist = calEuclidean(coord)
            % coord - n*2���� ����ֵ
            % euclDist - n*n���� ��������
            num = length(coord);
            euclDist = zeros([num, num]);
            for iedge = 1:num
                % create full matrix with edge weights (distances), every distance twice,
                Diff1 = repmat(coord(iedge,:), [num 1]);
                Diff2 = coord;
                euclDist(iedge,:) = round(sqrt(sum(((Diff1 - Diff2).^2),2))); % NOTE ����������ŷʽ���벢ȡ�� ����ȡ���Ƿ����?
            end
        end
        
    end

%% Vehicle
% n*1���� �����복���
% (Capacity;FixCost;VariableCost;PointCost;MaxPoint)
InsDataArray.Vehicle.Capacity = ...
    repelem(ThisLibData.VEHICLE_Capacity,ThisLibData.VEHICLE_Number);
InsDataArray.Vehicle.FixCost =  ...
    repelem(ThisLibData.VEHICLE_FixCost,ThisLibData.VEHICLE_Number);
InsDataArray.Vehicle.VariableCost =  ...
    repelem(ThisLibData.VEHICLE_VariableCost,ThisLibData.VEHICLE_Number);
InsDataArray.Vehicle.PointCost =  ...
    repelem(ThisLibData.VEHICLE_FeePoint,ThisLibData.VEHICLE_Number);
InsDataArray.Vehicle.MaxPoint=  ...
    repelem(ThisLibData.VEHICLE_MaxPoint,ThisLibData.VEHICLE_Number);

%% Veh_Cus
% Cost���㳵�����þ���
% D.Depot.CustomerDistance'; % ��depot����������滻depot��customer���� FIXME
M_FixCost = repmat(InsDataArray.Vehicle.FixCost, 1, CustomerNumber); % ÿ����ͬ

%WAY1 NOTE ��depot����������
%{         M_VarCost = insData.Vehicle.VariableCost * insData.Depot.CustomerDistance'; }

%WAY2 NOTE ��depot������������
tmp = arrayfun(@(x) InsDataArray.Depot.ClusterDistance...
    (InsDataArray.Customer.IdxCluster(x)),...
    (1:CustomerNumber));
M_VarCost = InsDataArray.Vehicle.VariableCost * tmp;  

InsDataArray.Veh_Cus.Cost_Num = M_FixCost + M_VarCost; %�����

%��ͬM_cost �����������ĳɱ��� NOTE:���������ȫ�����û���}
%{ insData.Veh_Cus.Cost_Type= unique(insData.Veh_Cus.Cost_Num, 'rows'); 

%% Vehicle new (����unitCost)
InsDataArray.Vehicle.UnitCost =  ... 
    nanmean(InsDataArray.Veh_Cus.Cost_Num,2) ./ InsDataArray.Vehicle.Capacity;

fprintf('Updating file (%s) ... done \n',ThisLibData.NAME );

% TO DO 8 GET volume and weight ������������� Ŀǰ���� ���������???
%{ P_n = CustomerNumber;
% P_NbType = length(V_MaxCapacity);
% V_volVehicle = (ones(1,P_NbType)*P_n)';
% V_weightVehicle = (ones(1,P_NbType)*P_n)';
% V_volCustomer = (ones(1,P_n)*1);
%} V_weightCustomer =  (ones(1,P_n)*1);

end

%% 2 Local function 333 checkInsDataArray
%
%
% TODO �������º˲����ݵĺ���
function flag = checkInsDataArray(da)
flag = 1;
% �жϾ�����������Ƿ��д�
Flag1 = any( da.Cluster.Demand <= 0 );
Flag2 = sum(da.Cluster.Demand) ~= sum(da.Customer.Demand);
if Flag1 || Flag2, error('�������������ڴ���'); end

% �ж��ܷ����Ƿ���ڸ���
Flag1 = any( da.Veh_Cus.Cost_Num(:) <= 0 );
if Flag1,  error('�ܷ��ô��ڸ���'); end

Flag1 = any( da.Cluster.Demand > max(da.Vehicle.Capacity) );
if Flag1,  warning('���ھ���������>���������'); end

% �ж�M_cost��������Ƿ����� (δ����Cost_Type)
%         Flag2 = size(insData.Veh_Cus.Cost_Type,1) ~= length(unique(insData.Vehicle.Capacity));
%         Flag3 = size(insData.Veh_Cus.Cost_Type,2) ~= length(insData.Customer.Demand);
%         if Flag1 || Flag2 || Flag3, warning('�����������������ڴ���'); end

% Check Depot �� cluster ֮���Ƿ���С��0�ľ���
if any(da.Depot.ClusterDistance(:) <= 0) || any(da.Cluster.Distance(:) < 0)
    error('Depot��Cluster�ľ��� �� cluster֮����� ���� С��0');
end

%         % Check insData�ṹ���Ƿ��ɿ�ȱ�ֶ�
%         if any(structfun (@(x) any(structfun (@isempty,x)),insData))
%             error('Updating file has something wrong ! \n' );
%         end
end

