%% 主函数文件
% 今后的模板文件，用于美观Matlab发布。
%% function S = VRP_Main()
% 
% # par参数结构体初始化
% # readLibData
% # updateLibData
% 
function S = VRP_Main()
%   初始化
clc;  clear; close all;  format long g; format bank; %NOTE 不被MATLAB CODE 支持
rng(1); % NOTE 是否随机的标志

fileType = 'vrp';
pars = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
pars.ins = struct('pho',0.2,'gamma',0.2); %pho正比聚类数；gamma正比连通性
pars.alg = struct('name','BF');
% 获取lib数据: 从TempInstance中读取后缀为vrp的文件到nested结构体libDatas
libDatas = readLibData('TempInstance',fileType); %TODO 后期修改参数二未以提出author命名

% 更新lib数据: 从libDates获取更新后可用的nested结构体insDatas
insDatas = updateLibDatas(libDatas,pars.ins,fileType);
clear fileType libDatas;

% 获取排序后的ins数据（555 函数内指定车辆和顾客顺序）
insIdxDatas = arrayfun(@sortInsDatas,insDatas); 
% printstruct(insIdxDatas);

% 最终ins数据调用算法
soltions.instance = [];
soltions.algName = [];
soltions.algArg = [];
soltions.solWff = [];
soltions.solCost = [];
soltions.solFeasible = [];
soltions = repmat(soltions,[length(insIdxDatas) 1]);

soltions = arrayfun(@callHeuBF1,insIdxDatas,soltions);
% printstruct(soltions);

end

function solu = callHeuBF1(inst,solu) 

fprintf("Running algorithm BestFit ... \n ");
 
%迭代安排customer到vehicle
customerNumber = length(inst.Customer.Demand);
vehicleNumber = length(inst.Vehicle.Capacity);
wff = zeros(vehicleNumber, customerNumber); %wff 代表某列(供应商)放入某行(可用车数) 供应商数*V_NbvehicleV_Nbvehicle -------------------

for curCustomer = 1:customerNumber   
    % 找出可以放入customer i的车辆集合
    [idxAvailableNew,idxAvailableUsed,idxAvailable] = getAvailableVehicles();
    
    % 找出可以放入customer i的最合适的那个车辆
    bsetVehicle = findBestFitVehicle();                    

    % 将customer i放入车辆序号为idx的车,更新对应wff数据
    putIdxVehicle();
end

if ~all(sum(wff))
    error('不是所有顾客都放入车辆中');
else
    solu.algName = 'Best Fit Algorithm';
    solu.algArg = 'TBA';
    solu.instance = inst;
    solu.solFeasible = 1;
    solu.solWff = wff;
    solu.solCost = calCost(inst,wff);
    fprintf("Running algorithm BestFit ... done \n ");
end


%% 
%
% 嵌套函数
%
    function [idxAvailableNew,idxAvailableUsed,idxAvailable] = getAvailableVehicles()
        %找出满足容量约束的所有车辆
        %
        %计算车辆可用容量      
        Vehicle.AssignedCapa = wff * inst.Customer.Demand;
        Vehicle.ResidueCapa = inst.Vehicle.Capacity - Vehicle.AssignedCapa;
        
        % 约束1: 满足容量约束的Logical array
        idxCapacity = ( Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer) );
        idxNewVehicles = (Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer)) & ( Vehicle.ResidueCapa == inst.Vehicle.Capacity );
        idxUsedVehicles = (Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer)) & ( Vehicle.ResidueCapa < inst.Vehicle.Capacity );
        if isempty(idxCapacity), error('没有车辆满足(容量)要求，请增加可用车数量'); end
        
        % 约束2: 满足最大点数的Logical array
        AssignedPointsEachVeh = sum(wff,2);
        idxMaxPoint = ( (inst.Vehicle.MaxPoint - AssignedPointsEachVeh ) > 0 );
        if isempty(idxMaxPoint), error('没有车辆满足(最大点数)要求，请增加最大可访问点数'); end
        
%         inst.Customer.Compatible
        % 约束3: 满足顾客i与车辆内其他顾客连通性的Logical array
        customerIMat = repmat( inst.Customer.Compatible(curCustomer,:), [vehicleNumber 1]);
        idxCompatible = ( sum(wff.*customerIMat,2) == 0 );
        if isempty(idxCompatible), error('没有车辆满足(顾客连通性)要求，请增加顾客连通性'); end
        
        % 满足所有约束的车辆Logical array返回
        idxAvailable = [idxCapacity idxMaxPoint idxCompatible];
        idxAvailableNew = [idxNewVehicles idxMaxPoint idxCompatible];
        idxAvailableUsed = [idxUsedVehicles idxMaxPoint idxCompatible];
        idxAvailable = all(idxAvailable,2);
        idxAvailableNew = all(idxAvailableNew,2);
        idxAvailableUsed = all(idxAvailableUsed,2);
        
    end

    function idx = findBestFitVehicle()
        %计算车辆基价 放入当前Vehicle中所有Customer中最大的
        Vehicle.BaseCost = max( wff .* inst.Veh_Cus.Cost_Num, [], 2) ;
        AbsCost = abs( inst.Veh_Cus.Cost_Num(:,curCustomer) - Vehicle.BaseCost );
        
        %如果Used车辆存在,必须找一个放
        if any(idxAvailableUsed)
            minAbsCost = min(AbsCost(idxAvailableUsed));
        elseif any(idxAvailableNew)
            minAbsCost = min(AbsCost(idxAvailableNew));
        else
            error('没有任何车辆可用，请增加可用车数量');
        end
        %从新车或旧车中找都可以
        % %         minAbsCost = min(AbsCost(idxAvailable));
        
        idx = find(AbsCost == minAbsCost); 
        if isempty(idx), error('EEEEEE'); end
        if numel(idx)>1, idx = idx(1); end
        fprintf('Node %d (demand %d) in Vehicle %d (capacity = %d;maxpoint = %d) compatible = 未知  \n ',...
            curCustomer, inst.Customer.Demand(curCustomer), idx, inst.Vehicle.Capacity(idx) , inst.Vehicle.MaxPoint(idx))
    end

    function putIdxVehicle()
        wff(bsetVehicle,curCustomer) = 1;
        calCost(inst,wff); %NOTE 当算法继续运行需要Cost时必须运算;BF可不运算
    end

    function Cost = calCost(d,wff)
        Veh_Cus_BaseCost = wff .* d.Veh_Cus.Cost_Num;
        Veh_BaseCost = max(Veh_Cus_BaseCost,[],2);
        Veh_PointCost = d.Vehicle.PointCost .* ( sum(wff,2) - 1);
        Veh_PointCost(Veh_PointCost < 0)=0;
        Veh_TotalCost = Veh_BaseCost + Veh_PointCost;
        TotalCost = sum(Veh_TotalCost);
        if any(Veh_TotalCost<0) || any(Veh_TotalCost <0) || any(Veh_PointCost)<0
            fprintf("error cost");
        else
            Cost = struct('TotalCost',TotalCost,'Veh_BaseCost',Veh_BaseCost,...
                'Veh_PointCost',Veh_PointCost,'Veh_TotalCost',Veh_TotalCost);
            fprintf("TotalCost of thisInstance = %d \n ", TotalCost);
        end
    end

end

function insDataNew = sortInsDatas(insData)
    [~,idxCus] = sort(insData.Customer.Demand);
    [~,idxVeh] = sort(insData.Vehicle.Capacity);
    insDataNew = getIdxInstance(idxCus,idxVeh,insData);
    
        function d = getIdxInstance(idxCus,idxVeh,d)
        %算例d依据顾客顺序idxCus和车辆顺序idxVeh变化 TODO 后续可能要增加Cluster顺序变化
        d.Customer.Demand = d.Customer.Demand(idxCus,:);
        d.Customer.IdxCluster = d.Customer.IdxCluster(idxCus,:);
        d.Customer.Compatible = d.Customer.Compatible(idxCus,idxCus);
        d.Customer.Coord = d.Customer.Coord(idxCus,:);
        d.Customer.Distance = d.Customer.Distance(idxCus,idxCus);
        
        d.Depot.CustomerDistance = d.Depot.CustomerDistance(idxCus,:);
        
        d.Node.Demand = [d.Depot.Demand; d.Customer.Demand];
        d.Node.Coord = [d.Depot.Coord; d.Customer.Coord];   
        d.Node.Distance(1,:) = [ 0; d.Depot.CustomerDistance]';
        d.Node.Distance(:,1) = [ 0; d.Depot.CustomerDistance];
        d.Node.Distance(2:end,2:end) = d.Customer.Distance;

        d.Vehicle.Capacity = d.Vehicle.Capacity(idxVeh,:);
        d.Vehicle.FixCost = d.Vehicle.FixCost(idxVeh,:);
        d.Vehicle.MaxPoint = d.Vehicle.MaxPoint(idxVeh,:);
        d.Vehicle.PointCost = d.Vehicle.PointCost(idxVeh,:);
        d.Vehicle.VariableCost = d.Vehicle.VariableCost(idxVeh,:);

        d.Veh_Cus.Cost_Num = d.Veh_Cus.Cost_Num(idxVeh,idxCus);
        d.Veh_Cus.FixCost = d.Veh_Cus.FixCost(idxVeh,idxCus);
        d.Veh_Cus.VarCost = d.Veh_Cus.VarCost(idxVeh,idxCus);        
        d.Veh_Cus.Cost_Type= unique(d.Veh_Cus.Cost_Num, 'rows');
    end

end









%% 555 重要参数说明
% idec : 改变V_size M_cost M_depot wff的顺序 (修改与工件相关的变量顺序)    
    %- idxV_UsedVehicle - 已usedbin :vector
    %- idxV_NewVehicle - 新bin :vector
    %- k - 从众多vehicle找出最佳安排车辆
    %- i - 当前待安排的item
    %- M_BCapacity(2,k) - 找出最佳车辆j对应的最佳车型

% 输入: x: 节点顺序; y: 车辆顺序; inst: 原始数据;
% 输出: *sol*
% sol.wff = wff;  行:
% sol.M_BCapacity = M_BCapacity;
% sol.isFeasible = 1;

%%
% % % [idxRow,idxCol,~] = find( wff(:,sum(wff,1)>0)>0 );  tmp = [idxRow  idxCol];     tmp % length(tmp)
% % %%
% % function [sol] = calHeuBF(sol) 
% % % isFlagFirstRandom = 1; %default is 0 means 空车时选择第一辆车 1 means 空车时随机选择一辆车
% % %% 0 数据初始化 Initialization
% % % rand('seed',1);  
% % 
% % %%
% % idec = x;  %idec: 真正的item索引:非原索引 按照pop中个体(供应商)顺序
% % idec = 1:length(x); % 55555555 测试
% % jdec = y;
% % [~,irec]=sort(idec);  %  irec 要复原的时候用的顺序 is the recovery ordering irec只在启发式结束后的wff中使用
% % [~,jrec]=sort(jdec);  
% % 
% %     M_cost2 = inst.mcost;
% %     idecM_cost2=M_cost2(jdec,idec); %jdec
% %     
% %     M_depot = inst.ndepot;
% %     idecM_depot = M_depot(idec,idec); %**** 对兼容性矩阵也要做相应调整 *** 下面只改ud即可  % note that ud=V_size(idec); and V_size = ud(irec);
% % 
% %     V_size = inst.node(2,:);
% %     idecV_size=V_size(:,idec);   %ud 按GA给定调整的供应商顺序 是复原后的items顺序
% % 
% %     V_MaxItemNb_2 = inst.vehicle(4,:); %修改顺序
% %     jdecV_MaxItemNb_2 = V_MaxItemNb_2(:,jdec);
% %     
% % % V_size_splitFlag = inst.nsplitflag;
% % % idecV_size_splitFlag=V_size_splitFlag(:,idec); %160831 update 增加是否拆分工件判断
% % % M_BCapacity = inst.vehicle(1:3,:);
% % % M_BCapacity = y;      % 与车型顺序相关 V_rff,V_cff顺序已为调整后顺序; 车型顺序不变
% % %     M_cost = inst.vncost;
% % %     idecM_cost=M_cost(:,idec); % idecM_cost 按idec调整后的M_cost------------
% % %  V_MaxItemNb = inst.vitem' %170723 改V为P % P_MaxItemNb = inst.vitem(1); %数组和单值的关系
% % %  V_MaxItemNb_2 = V_MaxItemNb(M_BCapacity(2,:)) %V_MaxItemNb_2: 对应车数量的最大可访问点数
% % 
% % n = size(inst.node,2);   %n: item数量
% % m = size(inst.vehicle,2);  %m: bin数量
% % wff=zeros(n, m); %wff 代表某行(供应商)放入某列(可用车数) 供应商数*V_NbvehicleV_Nbvehicle -------------------
% % % 长度为1*V_Nbvehicle 代表车辆剩余容量和车辆基价费用
% % V_rff = inst.vehicle(1, :);  %V_rff向量 - 当前所有V_Nbvehicle的随时变化的剩余容量 初始化为最大容量 - 1*V_Nbvehicle ----------------------------
% % V_rff = V_rff(:, jdec);
% % V_cff = zeros(1,m);  % V_cff 当前所有V_Nbvehicle的随时变化的cost 初始化为0的费用 1*V_Nbvehicle ----------------------------
% % %% NEW - 1 增加部分解partsol参数后的代码
% % tmpN = [];
% % if ~isempty(partsol) %更新部分解的V_cff 和 V_rff; 还需更新工件数量,不能是n个了
% %     wff = partsol.wff(idec, jdec);  %部分解wff顺序需更新
% % %     M_BCapacity = partsol.M_BCapacity; %应该把M_BCapacity变成y,rff和cff,V_maxNb也要跟着变,太困难, 暂不考虑把.
% % %      V_rff=M_BCapacity(1,:);
% %     idxV_UsedVehicle = find(sum(wff,1)>0); %idxV_UsedVehicle:部分解中使用的车idx
% %     [idxRow,idxCol,~] = find( wff(:,idxV_UsedVehicle)>0 ); %idxRow: 子矩阵中item的索引, idxCol: 所在列
% %     for tmp = 1: max(idxCol)  %每辆used车更新容量和成本
% %         tmpN= [tmpN idxRow(idxCol==tmp)']; %已经安排的item索引号
% %         tmpCost = max(idecM_cost2( idxV_UsedVehicle(tmp), idxRow(idxCol==tmp) )); %555 车辆 与 工件 的最大基价 = 车辆基价 %     tmpCost = max(idecM_cost( M_BCapacity(2, idxV_UsedVehicle(tmp)), idxRow(idxCol==tmp) )); %555 车型(单值) 与 工件(向量) 的最大基价 = 车辆基价
% %         V_cff(idxV_UsedVehicle(tmp)) = max(V_cff(idxV_UsedVehicle(tmp)), tmpCost); %更新车辆基价费用
% %                         %         idxRow(idxCol==tmp)
% %                         %         idecV_size(idxRow(idxCol==tmp))
% %                         %         sum(idecV_size(idxRow(idxCol==tmp)))
% %                         %         V_rff(idxV_UsedVehicle(tmp))
% %                         %         sum(wff)
% %         V_rff(idxV_UsedVehicle(tmp)) = V_rff(idxV_UsedVehicle(tmp)) - sum(idecV_size(idxRow(idxCol==tmp))); % 555 剩余容量 = 原容量-已安排item容量和
% %         if length(V_rff(idxV_UsedVehicle(tmp))) ~= length(sum(idecV_size(idxRow(idxCol==tmp)))) error('EEE'); end
% %     end
% % %     V_rff
% % %     V_cff
% % end
% % % V_rff(V_rff<0)
% % NoVehicleFlag = 0;       %- 1代表无车可用 如出现1 代表非可行解
% % 
% % %%
% % for i=1:n, %从p的新顺序逐个安排供应商
% %     if ismember(i, tmpN)
% %         continue;
% %     end
% %     %% Step1 逐步筛选出i可放入的Used过的车辆 : idxV_UsedVehicle
% %     idxV_UsedVehicle = find(sum(wff,1)>0); %1.1 所有used车辆集合
% %     if any(idxV_UsedVehicle) %如果非空
% %         idxV_UsedVehicle = idxV_UsedVehicle( V_rff(idxV_UsedVehicle) >= idecV_size(i) );  %1.2 约束1: 车辆剩余容量>=i尺寸
% %         if any(idxV_UsedVehicle)
% %             idxV_UsedVehicle =  idxV_UsedVehicle( sum(wff(:,idxV_UsedVehicle)) <= jdecV_MaxItemNb_2(idxV_UsedVehicle) );    %1.3 约束2: 每辆车已有items数<= 对应车型最大点数
% %                  if length(jdecV_MaxItemNb_2(idxV_UsedVehicle)) ~= length(sum(wff(:,idxV_UsedVehicle))) error('1.3'); end
% %             if any(idxV_UsedVehicle)
% %                 [idxRow,idxCol,~] = find( wff(:,idxV_UsedVehicle)>0 ); %idxRow: 子矩阵中item的索引, idxCol: 所在列
% %                 tmpV = [];
% %                 for tmp = 1: max(idxCol)
% %                     if ~all( idecM_depot(i, idxRow(idxCol==tmp)) ==0 )  %如果不是全部连通
% %                         tmpV = [tmpV tmp];
% %                     end
% %                 end
% %                 idxV_UsedVehicle(tmpV) = [];  %1.4 约束3: 去除不能连通的车辆
% %             end
% %         end
% %     end
% %     %% Step 2: 是否UsedVehicle内存在车辆
% %     if any(idxV_UsedVehicle) %如果存在: 那就选择最好的best bin; 否则随机选择一个bin        
% %         [~,k]=min(abs( V_cff(idxV_UsedVehicle) - idecM_cost2(idxV_UsedVehicle,i)' ) ); %% 555 最重要的条件(Best Fit) %         [~,k]=min(abs( V_cff(idxV_UsedVehicle) - idecM_cost(M_BCapacity(2,idxV_UsedVehicle),i)' ) ); %% 555 最重要的条件(Best Fit)
% %         k = idxV_UsedVehicle(k);                  %第k辆车被选中
% % %         disp('From Used Bin');
% %     else
% %         %找出new bin, 即新车
% %         idxV_NewVehicle = find(sum(wff,1)==0); %向量: 空车且满足容量约束
% %         idxV_NewVehicle = idxV_NewVehicle( V_rff(idxV_NewVehicle) >= idecV_size(i) );  %防止新车容量 < item容量情况
% %         if ~any(idxV_NewVehicle), 
% %             NoVehicleFlag = 1;
% %             warning('Not Enough Vehicles');
% %             break;
% %         end  %车型不够在此出错55555
% %         k = randsample(length(idxV_NewVehicle),1);
% %         k = idxV_NewVehicle(k);                                  %第k辆新车被选中
% % %         disp('From New Bin');
% %     end
% %     
% %     %% Step 3:i放入k车中并更新    
% %     wff(i, k) = 1;    %- 将该item(i)放入该bin(idxCol)
% %     V_rff(k)= V_rff(k) - idecV_size(i); %- 更新容量
% %     V_cff(k) = max( idecM_cost2( k, i),  V_cff(k));  %- 更新车辆基价 %      V_cff(k) = max( idecM_cost( M_BCapacity(2,k), i ),  V_cff(k))  %- 更新车辆基价
% %                if find(V_rff<0),        error('Capacity constraint is broken!');     end %-破坏容量约束不应该遇到
% % %                if M_BCapacity(2,k) > length(V_MaxItemNb), error('M_BCapacity constraint is broken!'); end    
% % end
% % 
% % % wff复原
% % wff = wff(irec,jrec); % sort it back to the original assignment order
% % 
% % sol.wff = wff;
% % % sol.M_BCapacity = inst.vehicle(1:2, :);
% % sol.isFeasible = 1;
% % 
% % if  NoVehicleFlag,
% %     sol.isFeasible =  0;
% % end
% % 
% % if ~(sum(sum(wff,2))==size(M_cost2,2)),
% %     sol.isFeasible =  0;
% % end
% % 
% % %  sol.isFeasible
% % end
