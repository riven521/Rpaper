function insDatas = updateLibDatas(libDatas,insArg,fileType)
%updateLibDatas������Lib�����ı�����
%   libDatas
%   insArg
%   fileType
%%
% *NOTE - ����: ��ʼ��Ƕ����ṹ��* 
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
% �ж�lib����: �Ƿ�Ϊ��������
if all(arrayfun(@checkInsDatas,insDatas)), fprintf("All data is corrected ! "); end

%% function VRP_Update(libData,insData,arg)
% 
% 
    function insData = VRP_Update(libData,insData,arg)
    % NOTE ��Ҫ����: ����VRP's Lib����
        
       %% 
        % Get Node(Demand;Coord); Cluster(Demand;Coord); Customer(Coor; Demand);
       %% Number
        NodeNumber = length(libData.Node_Coord);
        DepotNumber = 1; % FIXME:��Depot��ֻһ������λ��0/1ʱҪ�޸�;��DEPOT_SECTION����
        CustomerNumber = NodeNumber - DepotNumber;
        % VehicleNumber = sum(libData.VEHICLE_Number); %libData.VEHICLE_Number��ÿ�����ͳ���������
       %% Coord
        insData.Node.Coord = libData.Node_Coord;
        insData.Depot.Coord = libData.Node_Coord(DepotNumber,:);
        insData.Customer.Coord = libData.Node_Coord(DepotNumber+1:end,:);
       %% Demand
        insData.Node.Demand = libData.Node_Demand;
        insData.Depot.Demand = insData.Node.Demand(DepotNumber,:);
        insData.Customer.Demand = insData.Node.Demand(DepotNumber+1:end,:);
        
       %% Cluster 
        % ʹ�ò���pho
        [insData.Cluster.Coord, insData.Cluster.Demand, insData.Customer.IdxCluster ] = getCluster();
        
        function [cluCoord,cluDemand,idx] = getCluster()
            n = getClusterN(); 
            [idx, cluCoord,~,~] = kmeans(insData.Customer.Coord, n); % 2 NOTE K-means�㷨��ȡ��������Coord
            cluDemand = arrayfun(@(x) sum(insData.Customer.Demand (idx == x)), (1:n)' ); % 5 ����ÿ������������� Demand;
            
            % ����ÿ�����������������ʱ���ã�
            % NumberOneCluster = arrayfun(@(x) sum( D.Customer.IdxCluster == x ),(1:n)');
            
           %% function getClusterN()
            % 
            %              
            function ClusterNumber = getClusterN()
                % NOTE ˵�� ���������ɼ���{
                % ���ȷ�������� - ����2016CIE����pho*Q;
                % ���ȷ�������������;
                %} case1: ��������<=���������; case2: ��������ȫ��>��������� case3: �����������ڻ�С�����������
                
                minClusterNumber = 1;
                maxClusterNumber = NodeNumber - 1;
                
                pho = arg.pho;
                if (pho > 0 && pho < 1)
                    %   ClusterNumber = ceil(pho*sum(TSP_Demand(2:end)')/TSP_Capacity);
                    ClusterNumber = ceil(pho * maxClusterNumber);   % NOTE ClusterNumber ������ ���Բ��� para.pho [0,1]) ���������С1��,��಻������ĸ���(��������)
                elseif pho == 0 %���pho=0, ֻ��һ������; pho��Ϊ0Ҳ����ֻ��һ������
                    ClusterNumber = minClusterNumber;
                elseif pho == 1 %isinf(pho)
                    ClusterNumber = maxClusterNumber;
                else
                    error('pho = %d, ����[0,1]֮��', pho);
                end
                
                if (ClusterNumber < minClusterNumber || ClusterNumber > maxClusterNumber)
                    warning('ClusterNumber = %d, ������Χ', ClusterNumber);
                else
                    %                 fprintf('ClusterNumber = %d is in [1, %d] \n',ClusterNumber,maxClusterNumber);
                end
            end
            
        end
        
       %% Compatible 
        % ʹ�ò���gamma
        [insData.Cluster.Compatible,insData.Customer.Compatible] = getCompatible();
       %% function getCompatible()
        %
        %
        function [M_Cludepot,M_depot] = getCompatible()
            % NOTE ˵�� ��ͨ�����ɼ���{ % gamma Խ��, 0�ĸ���Խ��, ��ͨ��Խ��
            % GET M_Cludepot  ������ͨ�� (0:����ͨ; 1:������ͨ)
            n = length(insData.Cluster.Demand);
            gamma = arg.gamma;  % ClusterNumber=5 % gamma = 0.5;
            % ����totalZeros :һ����Ҫ��0�ĸ���,���ڷ�Χ:[minZeros,maxZeros]
            % gamma = 0.5; Ӧǡ��ֻ��һ����ͨ��
            minZeros = n; %�Խ�����0������
            maxZeros = n * n;
            totalZeros = minZeros + gamma * (maxZeros - minZeros);
            
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
            % FIXME ȥ��for
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
            
            % �������������gamma,��֤���ݺ����Էǳ���Ҫ
            ifsymmetric1 = (M_Cludepot == M_Cludepot');
            ifsymmetric2 = (M_depot == M_depot');
            if (totalZeros < minZeros || totalZeros > maxZeros || numZeros < totalZeros)
                error('totalZeros = %d, ������Χ', totalZeros);
            elseif (gamma <0 || gamma >1)
                error('gamma = %d, ������Χ [0,1]', gamma);
            elseif (~all(ifsymmetric1(:)) && ~all(ifsymmetric2(:)))
                error('M_Cludepot��M_depot���ǶԳƾ���');
            else
                %             fprintf('CORRECT : totalZeros = %d in [%d, %d] and gamma = %d, \n',totalZeros,minZeros,maxZeros,gamma);
            end
        end
        
       %% Distance 
        % ����ŷʽ����
        [insData.Node.Distance,insData.Customer.Distance,...
            insData.Depot.CustomerDistance,insData.Cluster.Distance,...
            insData.Depot.ClusterDistance] = getDistance();
        
       %% function getDistance(libData,insData,arg)
        % 
        %         
        function [NodeNodeDist,CusCusdist,DepCusdist,CluCludist,DepCludist] = getDistance()
            
            NodeNodeDist = calEuclidean(insData.Node.Coord); %����Node����������ŷ�Ͼ���
            CusCusdist = NodeNodeDist(DepotNumber+1:end,DepotNumber+1:end); % Customer֮��ľ���(����Depot)
            DepCusdist = NodeNodeDist(DepotNumber+1:end,1); % Depot����Custoemr��ľ���(����depot����)
            
            ClusterNode_Coord = [insData.Depot.Coord; insData.Cluster.Coord];
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
        % n*1���� �����복����أ���Customer�޹�
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
        % Cost���㳵�����þ���
        % D.Depot.CustomerDistance'; % ��depot����������滻depot��customer���� FIXME
        insData.Veh_Cus.FixCost = repmat(insData.Vehicle.FixCost, 1, CustomerNumber); % ÿ����ͬ
        insData.Veh_Cus.VarCost = insData.Vehicle.VariableCost * insData.Depot.CustomerDistance'; %NOTE ����������
        insData.Veh_Cus.Cost_Num = insData.Veh_Cus.FixCost + insData.Veh_Cus.VarCost; %�����
        insData.Veh_Cus.Cost_Type= unique(insData.Veh_Cus.Cost_Num, 'rows'); %��ͬM_cost �����������ĳɱ���
        
       %% fprintf
        fprintf('Updating file (%s) ... done \n',libData.NAME );
        
        %{ TO DO 8 GET volume and weight ������������� Ŀǰ���� ���������???
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
    % TODO �������º˲����ݵĺ���
    function flag = checkInsDatas(insData)
        flag = 1;
        % �жϾ�����������Ƿ��д�
        Flag1 = any( insData.Cluster.Demand <= 0 );
        Flag2 = any( insData.Cluster.Demand > max(insData.Vehicle.Capacity) );
        Flag3 = sum(insData.Cluster.Demand) ~= sum(insData.Customer.Demand);
        if Flag1 || Flag2 || Flag3, warning('�������������ڴ���'); end
        
        % �ж��ܷ����Ƿ���ڸ���
        Flag1 = any( insData.Veh_Cus.Cost_Num(:) <= 0 ); %�ܷ��ô��ڸ���
        
        % �ж�M_cost��������Ƿ�����
        Flag2 = size(insData.Veh_Cus.Cost_Type,1) ~= length(unique(insData.Vehicle.Capacity));
        Flag3 = size(insData.Veh_Cus.Cost_Type,2) ~= length(insData.Customer.Demand);
        if Flag1 || Flag2 || Flag3, warning('�����������������ڴ���'); end
        
        % Check Depot �� cluster ֮���Ƿ���С��0�ľ���
        if any(insData.Depot.ClusterDistance(:) <= 0) || any(insData.Cluster.Distance(:) < 0)
            warning('Depot��Cluster�ľ��� �� cluster֮����� ���� С��0');
            flag = 0;
        end
        
        % Check insData�ṹ���Ƿ��ɿ�ȱ�ֶ�
        if any(structfun (@(x) any(structfun (@isempty,x)),insData))
            error('Updating file has something wrong ! \n' );
        end
        
        
    end
end

