function insDatas = updateLibDatas(libDatas,insArg,fileType)
%updateLibDatas：更新Lib类型文本数据
%   libDatas
%   insArg
%   fileType
%%
% *NOTE - 技巧: 初始化嵌套类结构体* 
%
%% function readLibData(libFolder,fileType)
% 
% # VRP_Update
% # CCVRP_Update
% # checkInsDatas
% 
%% 
insDatas.Node = struct('Coord',[],'Demand',[],'Distance',[]);
insDatas.Customer = struct('Coord',[],'Demand',[],'IdxCluster',[],'Distance',[],'Compatible',[]);
insDatas.Depot = struct('Coord',[],'Demand',[],'CustomerDistance',[],'ClusterDistance',[]);
insDatas.Cluster = struct('Coord',[],'Demand',[],'Distance',[],'Compatible',[]);
insDatas.Vehicle = struct('Capacity',[],'FixCost',[],'VariableCost',[],...
    'PointCost',[],'MaxPoint',[]);
insDatas.Veh_Cus = struct('FixCost',[],'VarCost',[],'Cost_Num',[],'Cost_Type',[]);

insDatas = repmat(insDatas,[length(libDatas) 1]);
insArgs = repmat(insArg,[length(libDatas) 1]);

%%
if strcmp('vrp', fileType)
    insDatas = arrayfun(@VRP_Update,libDatas,insDatas,insArgs); 
elseif strcmp('ccvrp', fileType)
    insDatas = arrayfun(@CCVRP_Update,libDatas,insDatas,insArgs);
end
fprintf('Updating file with updateLibDatas() ... done \n');

%%
% 判断lib数据: 是否为合理数据
if all(arrayfun(@checkInsDatas,insDatas)), fprintf("All data is corrected ! "); end

%% function VRP_Update(libData,insData,arg)
% 
% 
    function insData = VRP_Update(libData,insData,arg)
    % NOTE 重要函数: 更新VRP's Lib数据
        
       %% 
        % Get Node(Demand;Coord); Cluster(Demand;Coord); Customer(Coor; Demand);
       %% Number
        NodeNumber = length(libData.Node_Coord);
        DepotNumber = 1; % FIXME:当Depot不只一个或不在位置0/1时要修改;按DEPOT_SECTION计数
        CustomerNumber = NodeNumber - DepotNumber;
        % VehicleNumber = sum(libData.VEHICLE_Number); %libData.VEHICLE_Number是每种类型车辆的数量
       %% Coord
        insData.Node.Coord = libData.Node_Coord;
        insData.Depot.Coord = libData.Node_Coord(DepotNumber,:);
        insData.Customer.Coord = libData.Node_Coord(DepotNumber+1:end,:);
       %% Demand
        insData.Node.Demand = libData.Node_Demand;
        insData.Depot.Demand = insData.Node.Demand(DepotNumber,:);
        insData.Customer.Demand = insData.Node.Demand(DepotNumber+1:end,:);
        
       %% Cluster 
        % 使用参数pho
        [insData.Cluster.Coord, insData.Cluster.Demand, insData.Customer.IdxCluster ] = getCluster();
        
        function [cluCoord,cluDemand,idx] = getCluster()
            n = getClusterN(); 
            [idx, cluCoord,~,~] = kmeans(insData.Customer.Coord, n); % 2 NOTE K-means算法获取聚类质心Coord
            cluDemand = arrayfun(@(x) sum(insData.Customer.Demand (idx == x)), (1:n)' ); % 5 计算每个聚类的总需求 Demand;
            
            % 计算每个聚类的总数量（暂时不用）
            % NumberOneCluster = arrayfun(@(x) sum( D.Customer.IdxCluster == x ),(1:n)');
            
           %% function getClusterN()
            % 
            %              
            function ClusterNumber = getClusterN()
                % NOTE 说明 聚类数生成技巧{
                % 如何确定聚类数 - 文献2016CIE采用pho*Q;
                % 如何确定聚类最大容量;
                %} case1: 聚类容量<=最大车型容量; case2: 聚类容量全部>最大车型容量 case3: 聚类容量大于或小于最大车辆容量
                
                minClusterNumber = 1;
                maxClusterNumber = NodeNumber - 1;
                
                pho = arg.pho;
                if (pho > 0 && pho < 1)
                    %   ClusterNumber = ceil(pho*sum(TSP_Demand(2:end)')/TSP_Capacity);
                    ClusterNumber = ceil(pho * maxClusterNumber);   % NOTE ClusterNumber 聚类数 来自参数 para.pho [0,1]) 聚类个数最小1个,最多不超过点的个数(单独成类)
                elseif pho == 0 %如果pho=0, 只有一个聚类; pho不为0也可能只有一个聚类
                    ClusterNumber = minClusterNumber;
                elseif pho == 1 %isinf(pho)
                    ClusterNumber = maxClusterNumber;
                else
                    error('pho = %d, 不在[0,1]之间', pho);
                end
                
                if (ClusterNumber < minClusterNumber || ClusterNumber > maxClusterNumber)
                    warning('ClusterNumber = %d, 超出范围', ClusterNumber);
                else
                    %                 fprintf('ClusterNumber = %d is in [1, %d] \n',ClusterNumber,maxClusterNumber);
                end
            end
            
        end
        
       %% Compatible 
        % 使用参数gamma
        [insData.Cluster.Compatible,insData.Customer.Compatible] = getCompatible();
       %% function getCompatible()
        %
        %
        function [M_Cludepot,M_depot] = getCompatible()
            % NOTE 说明 连通性生成技巧{ % gamma 越大, 0的个数越多, 连通性越好
            % GET M_Cludepot  聚类连通性 (0:可连通; 1:不可连通)
            n = length(insData.Cluster.Demand);
            gamma = arg.gamma;  % ClusterNumber=5 % gamma = 0.5;
            % 计算totalZeros :一共需要的0的个数,属于范围:[minZeros,maxZeros]
            % gamma = 0.5; 应恰好只有一半连通性
            minZeros = n; %对角线上0的数量
            maxZeros = n * n;
            totalZeros = minZeros + gamma * (maxZeros - minZeros);
            
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
            % FIXME 去除for
            M_depot = zeros(CustomerNumber,CustomerNumber);
            for tempi = 1 : CustomerNumber - 1
                for tempj = tempi + 1 : CustomerNumber
                    clusteri = insData.Customer.IdxCluster(tempi);
                    clusterj = insData.Customer.IdxCluster(tempj);
                    if M_Cludepot(clusteri,clusterj) == 1
                        M_depot(tempi,tempj) = 1;
                        M_depot(tempj,tempi) = 1;
                    end
                end
            end
            
            % 由于有输入参数gamma,验证数据合理性非常重要
            ifsymmetric1 = (M_Cludepot == M_Cludepot');
            ifsymmetric2 = (M_depot == M_depot');
            if (totalZeros < minZeros || totalZeros > maxZeros || numZeros < totalZeros)
                error('totalZeros = %d, 超出范围', totalZeros);
            elseif (gamma <0 || gamma >1)
                error('gamma = %d, 超出范围 [0,1]', gamma);
            elseif (~all(ifsymmetric1(:)) && ~all(ifsymmetric2(:)))
                error('M_Cludepot或M_depot不是对称矩阵');
            else
                %             fprintf('CORRECT : totalZeros = %d in [%d, %d] and gamma = %d, \n',totalZeros,minZeros,maxZeros,gamma);
            end
        end
        
       %% Distance 
        % 计算欧式距离
        [insData.Node.Distance,insData.Customer.Distance,...
            insData.Depot.CustomerDistance,insData.Cluster.Distance,...
            insData.Depot.ClusterDistance] = getDistance();
        
       %% function getDistance(libData,insData,arg)
        % 
        %         
        function [NodeNodeDist,CusCusdist,DepCusdist,CluCludist,DepCludist] = getDistance()
            
            NodeNodeDist = calEuclidean(insData.Node.Coord); %计算Node坐标两两的欧氏距离
            CusCusdist = NodeNodeDist(DepotNumber+1:end,DepotNumber+1:end); % Customer之间的距离(不含Depot)
            DepCusdist = NodeNodeDist(DepotNumber+1:end,1); % Depot到各Custoemr间的距离(不含depot本身)
            
            ClusterNode_Coord = [insData.Depot.Coord; insData.Cluster.Coord];
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
        % n*1向量 数据与车型相关，与Customer无关
        % (Capacity;FixCost;VariableCost;PointCost;MaxPoint)
        insData.Vehicle.Capacity = ...
            repelem(libData.VEHICLE_Capacity,libData.VEHICLE_Number);
        insData.Vehicle.FixCost =  ...
            repelem(libData.VEHICLE_FixCost,libData.VEHICLE_Number);
        insData.Vehicle.VariableCost =  ...
            repelem(libData.VEHICLE_VariableCost,libData.VEHICLE_Number);
        insData.Vehicle.PointCost =  ...
            repelem(libData.VEHICLE_FeePoint,libData.VEHICLE_Number);
        insData.Vehicle.MaxPoint=  ...
            repelem(libData.VEHICLE_MaxPoint,libData.VEHICLE_Number);
       %% Veh_Cus
        % Cost计算车辆费用矩阵
        % D.Depot.CustomerDistance'; % 以depot到聚类距离替换depot到customer距离 FIXME
        insData.Veh_Cus.FixCost = repmat(insData.Vehicle.FixCost, 1, CustomerNumber); % 每列相同
        insData.Veh_Cus.VarCost = insData.Vehicle.VariableCost * insData.Depot.CustomerDistance'; %NOTE 与距离正相关
        insData.Veh_Cus.Cost_Num = insData.Veh_Cus.FixCost + insData.Veh_Cus.VarCost; %简单相加
        insData.Veh_Cus.Cost_Type= unique(insData.Veh_Cus.Cost_Num, 'rows'); %等同M_cost （按车辆类别的成本）
        
       %% fprintf
        fprintf('Updating file (%s) ... done \n',libData.NAME );
        
        %{ TO DO 8 GET volume and weight 构建体积和重量 目前无用 但必须放入???
        % P_n = CustomerNumber;
        % P_NbType = length(V_MaxCapacity);
        % V_volVehicle = (ones(1,P_NbType)*P_n)';
        % V_weightVehicle = (ones(1,P_NbType)*P_n)';
        % V_volCustomer = (ones(1,P_n)*1);
        %} V_weightCustomer =  (ones(1,P_n)*1);
    end

%% function checkInsDatas(libData,insData,arg)
% 
% 
    % TODO 完善如下核查数据的函数
    function flag = checkInsDatas(insData)
        flag = 1;
        % 判断聚类需求计算是否有错
        Flag1 = any( insData.Cluster.Demand <= 0 );
        Flag2 = any( insData.Cluster.Demand > max(insData.Vehicle.Capacity) );
        Flag3 = sum(insData.Cluster.Demand) ~= sum(insData.Customer.Demand);
        if Flag1 || Flag2 || Flag3, warning('聚类需求计算存在错误'); end
        
        % 判断总费用是否存在负数
        Flag1 = any( insData.Veh_Cus.Cost_Num(:) <= 0 ); %总费用存在负数
        
        % 判断M_cost矩阵计算是否有误
        Flag2 = size(insData.Veh_Cus.Cost_Type,1) ~= length(unique(insData.Vehicle.Capacity));
        Flag3 = size(insData.Veh_Cus.Cost_Type,2) ~= length(insData.Customer.Demand);
        if Flag1 || Flag2 || Flag3, warning('车辆费用需求计算存在错误'); end
        
        % Check Depot 与 cluster 之间是否由小于0的距离
        if any(insData.Depot.ClusterDistance(:) <= 0) || any(insData.Cluster.Distance(:) < 0)
            warning('Depot到Cluster的距离 或 cluster之间距离 存在 小于0');
            flag = 0;
        end
        
        % Check insData结构体是否由空缺字段
        if any(structfun (@(x) any(structfun (@isempty,x)),insData))
            error('Updating file has something wrong ! \n' );
        end
        
        
    end
end

