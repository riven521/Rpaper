function InsDataArray = updateLibData(LibDataArray, ArgIns, fileType)
%updateLibDatas  更新LibDataArray为InsDataArray
% 输入
%   LibDataArray     lib文件读取的数据结构体 需要更新
%   ArgIns                instance相关参数结构体 PHO GAMMA
%   fileType             fileType文件类型
% 输出
%   InsDataArray    数据集结构体
%% function updateLibData(LibDataArray, ArgIns, fileType)
% Local Functions:
% # VRP_Update
% # CCVRP_Update    TODO 待写
% # checkInsDataArray

%% 
%{ InsDataArray.Node = struct('Coord',[],'Demand',[],'Distance',[]);
% InsDataArray.Customer = struct('Coord',[],'Demand',[],'IdxCluster',[],'Distance',[],'Compatible',[]);
% InsDataArray.Depot = struct('Coord',[],'Demand',[],'CustomerDistance',[],'ClusterDistance',[]);
% InsDataArray.Cluster = struct('Coord',[],'Demand',[],'Distance',[],'Compatible',[]);
% InsDataArray.Vehicle = struct('Capacity',[],'FixCost',[],'VariableCost',[],...
%     'PointCost',[],'MaxPoint',[],'Route',[],'Cost',[]); %Route来自wff
% InsDataArray.Veh_Cus = struct('Cost_Num',[],'Cost_Type',[],'wff',[]);
%} InsDataArray = repmat(InsDataArray,[length(LibDataArray) 1]);
%%
% 更新LibDataArray为InsDataArray
ArgInsArray  = repmat(ArgIns,[length(LibDataArray) 1]);
if strcmp('vrp', fileType)
    InsDataArray = arrayfun(@VRP_Update,LibDataArray,ArgInsArray);   %     InsDataArray = arrayfun(@VRP_Update,LibDataArray,InsDataArray,ArgInsArray); 
elseif strcmp('ccvrp', fileType)
    InsDataArray = arrayfun(@CCVRP_Update,LibDataArray,ArgInsArray);
end

% 判断InsDataArray是否合理
if all(arrayfun(@checkInsDataArray,InsDataArray)),  fprintf("All data is corrected ! \n ");  end

fprintf('Updating file with updateLibDatas() ... done \n');

end

%%  % 局部 FUNCTIONS %
% ==========================================================================
% 局部 FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Local function 555 VRP_Update
%
%
function InsDataArray = VRP_Update(ThisLibData,ThisArgIns)
% NOTE 重要函数: 更新VRP's Lib数据

%%
% Get Node(Demand;Coord); Cluster(Demand;Coord); Customer(Coor; Demand);
%% Number
NodeNumber = length(ThisLibData.Node_Coord);
DepotNumber = 1; % FIXME:当Depot不只一个或不在位置0/1时要修改;按DEPOT_SECTION计数
CustomerNumber = NodeNumber - DepotNumber;
% VehicleNumber = sum(libData.VEHICLE_Number); %libData.VEHICLE_Number是每种类型车辆的数量
%% Coord
InsDataArray.Node.Coord = ThisLibData.Node_Coord;
InsDataArray.Depot.Coord = ThisLibData.Node_Coord(DepotNumber,:);
InsDataArray.Customer.Coord = ThisLibData.Node_Coord(DepotNumber+1:end,:);
%% Demand
InsDataArray.Node.Demand = ThisLibData.Node_Demand;
InsDataArray.Depot.Demand = InsDataArray.Node.Demand(DepotNumber,:);
InsDataArray.Customer.Demand = InsDataArray.Node.Demand(DepotNumber+1:end,:);

%% Cluster
% 使用参数PHO
[InsDataArray.Cluster.Coord, InsDataArray.Cluster.Demand, InsDataArray.Customer.IdxCluster ] = getCluster();
InsDataArray.Depot.IdxCluster = max(InsDataArray.Customer.IdxCluster)+1;  %FIXME 或许can设置depot的聚类为第0个/可设置为第MAX+1个

    function [cluCoord,cluDemand,idx] = getCluster()
        n = getClusterN();
        [idx, cluCoord,~,~] = kmeans(InsDataArray.Customer.Coord, n); % 2 NOTE K-means算法获取聚类质心Coord
        cluDemand = arrayfun(@(x) sum(InsDataArray.Customer.Demand (idx == x)), (1:n)' ); % 5 计算每个聚类的总需求 Demand;
        
        % 计算每个聚类的总数量（暂时不用）
        % NumberOneCluster = arrayfun(@(x) sum( D.Customer.IdxCluster == x ),(1:n)');
        
        %% 嵌套 function getClusterN()
        %
        %
        function ClusterNumber = getClusterN()
            % NOTE 说明 聚类数生成技巧{
            % 如何确定聚类数 - 文献2016CIE采用PHO*Q;
            % 如何确定聚类最大容量;
            %} case1: 聚类容量<=最大车型容量; case2: 聚类容量全部>最大车型容量 case3: 聚类容量大于或小于最大车辆容量
            
            minClusterNumber = 1;
            maxClusterNumber = NodeNumber - 1;
            
            PHO = ThisArgIns.PHO;
            if (PHO > 0 && PHO < 1)
                %   ClusterNumber = ceil(PHO*sum(TSP_Demand(2:end)')/TSP_Capacity);
                ClusterNumber = ceil(PHO * maxClusterNumber);   % NOTE ClusterNumber 聚类数 来自参数 para.PHO [0,1]) 聚类个数最小1个,最多不超过点的个数(单独成类)
            elseif PHO == 0 %如果PHO=0, 只有一个聚类; PHO不为0也可能只有一个聚类
                ClusterNumber = minClusterNumber;
            elseif PHO == 1 %isinf(PHO)
                ClusterNumber = maxClusterNumber;
            else
                error('PHO = %d, 不在[0,1]之间', PHO);
            end
            
            if (ClusterNumber < minClusterNumber || ClusterNumber > maxClusterNumber)
                warning('ClusterNumber = %d, 超出范围', ClusterNumber);
            else
                %                 fprintf('ClusterNumber = %d is in [1, %d] \n',ClusterNumber,maxClusterNumber);
            end
        end
        
    end

%% Compatible
% 使用参数GAMMA
[InsDataArray.Cluster.Compatible,InsDataArray.Customer.Compatible] = getCompatible();

%% 嵌套 function getCompatible()
%
%
    function [M_Cludepot,M_depot] = getCompatible()
        % NOTE 说明 连通性生成技巧{ % GAMMA 越大, 0的个数越多, 连通性越好
        % GET M_Cludepot  聚类连通性 (0:可连通; 1:不可连通)
        n = length(InsDataArray.Cluster.Demand);
        if n==1
            M_Cludepot = 0;
            M_depot = zeros(CustomerNumber,CustomerNumber);
        else
            GAMMA = ThisArgIns.GAMMA;  % ClusterNumber=5 % GAMMA = 0.5;
            % 计算totalZeros :一共需要的0的个数,属于范围:[minZeros,maxZeros]
            % GAMMA = 0.5; 应恰好只有一半连通性
            minZeros = n; %对角线上0的数量
            maxZeros = n * n;
            totalZeros = minZeros + GAMMA * (maxZeros - minZeros);
            
            numZeros = minZeros;
            M_Cludepot = ones(n,n);
            M_Cludepot = M_Cludepot-diag(diag(M_Cludepot));
            
            % NOTE 技巧:构造连通性矩阵 M_Cludepot
            % 初始化和除对角线外均为1 循环增加0的个数until大于totalZeros
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
            
            
            % GET M_depot 点点连通性 基于聚类连通性的扩展
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
            
            % 由于有输入参数GAMMA,验证数据合理性非常重要
            ifsymmetric1 = (M_Cludepot == M_Cludepot');
            ifsymmetric2 = (M_depot == M_depot');
            if (totalZeros < minZeros || totalZeros > maxZeros || numZeros < totalZeros)
                error('totalZeros = %d, 超出范围', totalZeros);
            elseif (GAMMA <0 || GAMMA >1)
                error('GAMMA = %d, 超出范围 [0,1]', GAMMA);
            elseif (~all(ifsymmetric1(:)) && ~all(ifsymmetric2(:)))
                error('M_Cludepot或M_depot不是对称矩阵');
            else
                %             fprintf('CORRECT : totalZeros = %d in [%d, %d] and GAMMA = %d, \n',totalZeros,minZeros,maxZeros,GAMMA);
            end
        end
        
    end

%% Distance
% 计算欧式距离
[InsDataArray.Node.Distance,InsDataArray.Customer.Distance,...
    InsDataArray.Depot.CustomerDistance,InsDataArray.Cluster.Distance,...
    InsDataArray.Depot.ClusterDistance] = getDistance();

%% 嵌套 function getDistance(libData,insData,arg)
%
%
    function [NodeNodeDist,CusCusdist,DepCusdist,CluCludist,DepCludist] = getDistance()
        
        NodeNodeDist = calEuclidean(InsDataArray.Node.Coord); %计算Node坐标两两的欧氏距离
        CusCusdist = NodeNodeDist(DepotNumber+1:end,DepotNumber+1:end); % Customer之间的距离(不含Depot)
        DepCusdist = NodeNodeDist(DepotNumber+1:end,1); % Depot到各Custoemr间的距离(不含depot本身)
        
        ClusterNode_Coord = [InsDataArray.Depot.Coord; InsDataArray.Cluster.Coord];
        DepCluDepCluDist = calEuclidean(ClusterNode_Coord);
        CluCludist = DepCluDepCluDist(DepotNumber+1:end,DepotNumber+1:end);  % Clusters之间的距离(不含Depot)
        DepCludist = DepCluDepCluDist(DepotNumber+1:end,1);  % Depot到Clusters的距离(不含Depot)
        % 3 获得Depot到每个质心的距离-作为固定成本 TSP_ClusterEdgeWeight
        %% function calEuclidean(libData,insData,arg)
        %
        %
        function euclDist = calEuclidean(coord)
            % coord - n*2矩阵 坐标值
            % euclDist - n*n矩阵 两两距离
            num = length(coord);
            euclDist = zeros([num, num]);
            for iedge = 1:num
                % create full matrix with edge weights (distances), every distance twice,
                Diff1 = repmat(coord(iedge,:), [num 1]);
                Diff2 = coord;
                euclDist(iedge,:) = round(sqrt(sum(((Diff1 - Diff2).^2),2))); % NOTE 计算向量间欧式距离并取整 向上取整是否合适?
            end
        end
        
    end

%% Vehicle
% n*1向量 数据与车相关
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
% Cost计算车辆费用矩阵
% D.Depot.CustomerDistance'; % 以depot到聚类距离替换depot到customer距离 FIXME
M_FixCost = repmat(InsDataArray.Vehicle.FixCost, 1, CustomerNumber); % 每列相同

%WAY1 NOTE 与depot点距离正相关
%{         M_VarCost = insData.Vehicle.VariableCost * insData.Depot.CustomerDistance'; }

%WAY2 NOTE 与depot聚类距离正相关
tmp = arrayfun(@(x) InsDataArray.Depot.ClusterDistance...
    (InsDataArray.Customer.IdxCluster(x)),...
    (1:CustomerNumber));
M_VarCost = InsDataArray.Vehicle.VariableCost * tmp;  

InsDataArray.Veh_Cus.Cost_Num = M_FixCost + M_VarCost; %简单相加

%等同M_cost （按车辆类别的成本） NOTE:车辆类别完全按费用划分}
%{ insData.Veh_Cus.Cost_Type= unique(insData.Veh_Cus.Cost_Num, 'rows'); 

%% Vehicle new (增加unitCost)
InsDataArray.Vehicle.UnitCost =  ... 
    nanmean(InsDataArray.Veh_Cus.Cost_Num,2) ./ InsDataArray.Vehicle.Capacity;

fprintf('Updating file (%s) ... done \n',ThisLibData.NAME );

% TO DO 8 GET volume and weight 构建体积和重量 目前无用 但必须放入???
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
% TODO 完善如下核查数据的函数
function flag = checkInsDataArray(da)
flag = 1;
% 判断聚类需求计算是否有错
Flag1 = any( da.Cluster.Demand <= 0 );
Flag2 = sum(da.Cluster.Demand) ~= sum(da.Customer.Demand);
if Flag1 || Flag2, error('聚类需求计算存在错误'); end

% 判断总费用是否存在负数
Flag1 = any( da.Veh_Cus.Cost_Num(:) <= 0 );
if Flag1,  error('总费用存在负数'); end

Flag1 = any( da.Cluster.Demand > max(da.Vehicle.Capacity) );
if Flag1,  warning('存在聚类总需求>最大车辆容量'); end

% 判断M_cost矩阵计算是否有误 (未计算Cost_Type)
%         Flag2 = size(insData.Veh_Cus.Cost_Type,1) ~= length(unique(insData.Vehicle.Capacity));
%         Flag3 = size(insData.Veh_Cus.Cost_Type,2) ~= length(insData.Customer.Demand);
%         if Flag1 || Flag2 || Flag3, warning('车辆费用需求计算存在错误'); end

% Check Depot 与 cluster 之间是否由小于0的距离
if any(da.Depot.ClusterDistance(:) <= 0) || any(da.Cluster.Distance(:) < 0)
    error('Depot到Cluster的距离 或 cluster之间距离 存在 小于0');
end

%         % Check insData结构体是否由空缺字段
%         if any(structfun (@(x) any(structfun (@isempty,x)),insData))
%             error('Updating file has something wrong ! \n' );
%         end
end

