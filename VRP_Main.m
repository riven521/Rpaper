%% �������ļ�
% ����ģ���ļ�����������Matlab������
%% function S = VRP_Main()
% 
% # par�����ṹ���ʼ��
% # readLibData
% # updateLibData
% 
function S = VRP_Main()
%   ��ʼ��
clc;  clear; close all;  format long g; format bank; %NOTE ����MATLAB CODE ֧��
rng(1); % NOTE �Ƿ�����ı�־

fileType = 'vrp';
pars = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
pars.ins = struct('pho',0.2,'gamma',0.2); %pho���Ⱦ�������gamma������ͨ��
pars.alg = struct('name','BF');
% ��ȡlib����: ��TempInstance�ж�ȡ��׺Ϊvrp���ļ���nested�ṹ��libDatas
libDatas = readLibData('TempInstance',fileType); %TODO �����޸Ĳ�����δ�����author����

% ����lib����: ��libDates��ȡ���º���õ�nested�ṹ��insDatas
insDatas = updateLibDatas(libDatas,pars.ins,fileType);
clear fileType libDatas;

% ��ȡ������ins���ݣ�555 ������ָ�������͹˿�˳��
insIdxDatas = arrayfun(@sortInsDatas,insDatas); 
% printstruct(insIdxDatas);

% ����ins���ݵ����㷨
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
 
%��������customer��vehicle
customerNumber = length(inst.Customer.Demand);
vehicleNumber = length(inst.Vehicle.Capacity);
wff = zeros(vehicleNumber, customerNumber); %wff ����ĳ��(��Ӧ��)����ĳ��(���ó���) ��Ӧ����*V_NbvehicleV_Nbvehicle -------------------

for curCustomer = 1:customerNumber   
    % �ҳ����Է���customer i�ĳ�������
    [idxAvailableNew,idxAvailableUsed,idxAvailable] = getAvailableVehicles();
    
    % �ҳ����Է���customer i������ʵ��Ǹ�����
    bsetVehicle = findBestFitVehicle();                    

    % ��customer i���복�����Ϊidx�ĳ�,���¶�Ӧwff����
    putIdxVehicle();
end

if ~all(sum(wff))
    error('�������й˿Ͷ����복����');
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
% Ƕ�׺���
%
    function [idxAvailableNew,idxAvailableUsed,idxAvailable] = getAvailableVehicles()
        %�ҳ���������Լ�������г���
        %
        %���㳵����������      
        Vehicle.AssignedCapa = wff * inst.Customer.Demand;
        Vehicle.ResidueCapa = inst.Vehicle.Capacity - Vehicle.AssignedCapa;
        
        % Լ��1: ��������Լ����Logical array
        idxCapacity = ( Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer) );
        idxNewVehicles = (Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer)) & ( Vehicle.ResidueCapa == inst.Vehicle.Capacity );
        idxUsedVehicles = (Vehicle.ResidueCapa >= inst.Customer.Demand(curCustomer)) & ( Vehicle.ResidueCapa < inst.Vehicle.Capacity );
        if isempty(idxCapacity), error('û�г�������(����)Ҫ�������ӿ��ó�����'); end
        
        % Լ��2: ������������Logical array
        AssignedPointsEachVeh = sum(wff,2);
        idxMaxPoint = ( (inst.Vehicle.MaxPoint - AssignedPointsEachVeh ) > 0 );
        if isempty(idxMaxPoint), error('û�г�������(������)Ҫ�����������ɷ��ʵ���'); end
        
%         inst.Customer.Compatible
        % Լ��3: ����˿�i�복���������˿���ͨ�Ե�Logical array
        customerIMat = repmat( inst.Customer.Compatible(curCustomer,:), [vehicleNumber 1]);
        idxCompatible = ( sum(wff.*customerIMat,2) == 0 );
        if isempty(idxCompatible), error('û�г�������(�˿���ͨ��)Ҫ�������ӹ˿���ͨ��'); end
        
        % ��������Լ���ĳ���Logical array����
        idxAvailable = [idxCapacity idxMaxPoint idxCompatible];
        idxAvailableNew = [idxNewVehicles idxMaxPoint idxCompatible];
        idxAvailableUsed = [idxUsedVehicles idxMaxPoint idxCompatible];
        idxAvailable = all(idxAvailable,2);
        idxAvailableNew = all(idxAvailableNew,2);
        idxAvailableUsed = all(idxAvailableUsed,2);
        
    end

    function idx = findBestFitVehicle()
        %���㳵������ ���뵱ǰVehicle������Customer������
        Vehicle.BaseCost = max( wff .* inst.Veh_Cus.Cost_Num, [], 2) ;
        AbsCost = abs( inst.Veh_Cus.Cost_Num(:,curCustomer) - Vehicle.BaseCost );
        
        %���Used��������,������һ����
        if any(idxAvailableUsed)
            minAbsCost = min(AbsCost(idxAvailableUsed));
        elseif any(idxAvailableNew)
            minAbsCost = min(AbsCost(idxAvailableNew));
        else
            error('û���κγ������ã������ӿ��ó�����');
        end
        %���³���ɳ����Ҷ�����
        % %         minAbsCost = min(AbsCost(idxAvailable));
        
        idx = find(AbsCost == minAbsCost); 
        if isempty(idx), error('EEEEEE'); end
        if numel(idx)>1, idx = idx(1); end
        fprintf('Node %d (demand %d) in Vehicle %d (capacity = %d;maxpoint = %d) compatible = δ֪  \n ',...
            curCustomer, inst.Customer.Demand(curCustomer), idx, inst.Vehicle.Capacity(idx) , inst.Vehicle.MaxPoint(idx))
    end

    function putIdxVehicle()
        wff(bsetVehicle,curCustomer) = 1;
        calCost(inst,wff); %NOTE ���㷨����������ҪCostʱ��������;BF�ɲ�����
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
        %����d���ݹ˿�˳��idxCus�ͳ���˳��idxVeh�仯 TODO ��������Ҫ����Cluster˳��仯
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









%% 555 ��Ҫ����˵��
% idec : �ı�V_size M_cost M_depot wff��˳�� (�޸��빤����صı���˳��)    
    %- idxV_UsedVehicle - ��usedbin :vector
    %- idxV_NewVehicle - ��bin :vector
    %- k - ���ڶ�vehicle�ҳ���Ѱ��ų���
    %- i - ��ǰ�����ŵ�item
    %- M_BCapacity(2,k) - �ҳ���ѳ���j��Ӧ����ѳ���

% ����: x: �ڵ�˳��; y: ����˳��; inst: ԭʼ����;
% ���: *sol*
% sol.wff = wff;  ��:
% sol.M_BCapacity = M_BCapacity;
% sol.isFeasible = 1;

%%
% % % [idxRow,idxCol,~] = find( wff(:,sum(wff,1)>0)>0 );  tmp = [idxRow  idxCol];     tmp % length(tmp)
% % %%
% % function [sol] = calHeuBF(sol) 
% % % isFlagFirstRandom = 1; %default is 0 means �ճ�ʱѡ���һ���� 1 means �ճ�ʱ���ѡ��һ����
% % %% 0 ���ݳ�ʼ�� Initialization
% % % rand('seed',1);  
% % 
% % %%
% % idec = x;  %idec: ������item����:��ԭ���� ����pop�и���(��Ӧ��)˳��
% % idec = 1:length(x); % 55555555 ����
% % jdec = y;
% % [~,irec]=sort(idec);  %  irec Ҫ��ԭ��ʱ���õ�˳�� is the recovery ordering irecֻ������ʽ�������wff��ʹ��
% % [~,jrec]=sort(jdec);  
% % 
% %     M_cost2 = inst.mcost;
% %     idecM_cost2=M_cost2(jdec,idec); %jdec
% %     
% %     M_depot = inst.ndepot;
% %     idecM_depot = M_depot(idec,idec); %**** �Լ����Ծ���ҲҪ����Ӧ���� *** ����ֻ��ud����  % note that ud=V_size(idec); and V_size = ud(irec);
% % 
% %     V_size = inst.node(2,:);
% %     idecV_size=V_size(:,idec);   %ud ��GA���������Ĺ�Ӧ��˳�� �Ǹ�ԭ���items˳��
% % 
% %     V_MaxItemNb_2 = inst.vehicle(4,:); %�޸�˳��
% %     jdecV_MaxItemNb_2 = V_MaxItemNb_2(:,jdec);
% %     
% % % V_size_splitFlag = inst.nsplitflag;
% % % idecV_size_splitFlag=V_size_splitFlag(:,idec); %160831 update �����Ƿ��ֹ����ж�
% % % M_BCapacity = inst.vehicle(1:3,:);
% % % M_BCapacity = y;      % �복��˳����� V_rff,V_cff˳����Ϊ������˳��; ����˳�򲻱�
% % %     M_cost = inst.vncost;
% % %     idecM_cost=M_cost(:,idec); % idecM_cost ��idec�������M_cost------------
% % %  V_MaxItemNb = inst.vitem' %170723 ��VΪP % P_MaxItemNb = inst.vitem(1); %����͵�ֵ�Ĺ�ϵ
% % %  V_MaxItemNb_2 = V_MaxItemNb(M_BCapacity(2,:)) %V_MaxItemNb_2: ��Ӧ�����������ɷ��ʵ���
% % 
% % n = size(inst.node,2);   %n: item����
% % m = size(inst.vehicle,2);  %m: bin����
% % wff=zeros(n, m); %wff ����ĳ��(��Ӧ��)����ĳ��(���ó���) ��Ӧ����*V_NbvehicleV_Nbvehicle -------------------
% % % ����Ϊ1*V_Nbvehicle ������ʣ�������ͳ������۷���
% % V_rff = inst.vehicle(1, :);  %V_rff���� - ��ǰ����V_Nbvehicle����ʱ�仯��ʣ������ ��ʼ��Ϊ������� - 1*V_Nbvehicle ----------------------------
% % V_rff = V_rff(:, jdec);
% % V_cff = zeros(1,m);  % V_cff ��ǰ����V_Nbvehicle����ʱ�仯��cost ��ʼ��Ϊ0�ķ��� 1*V_Nbvehicle ----------------------------
% % %% NEW - 1 ���Ӳ��ֽ�partsol������Ĵ���
% % tmpN = [];
% % if ~isempty(partsol) %���²��ֽ��V_cff �� V_rff; ������¹�������,������n����
% %     wff = partsol.wff(idec, jdec);  %���ֽ�wff˳�������
% % %     M_BCapacity = partsol.M_BCapacity; %Ӧ�ð�M_BCapacity���y,rff��cff,V_maxNbҲҪ���ű�,̫����, �ݲ����ǰ�.
% % %      V_rff=M_BCapacity(1,:);
% %     idxV_UsedVehicle = find(sum(wff,1)>0); %idxV_UsedVehicle:���ֽ���ʹ�õĳ�idx
% %     [idxRow,idxCol,~] = find( wff(:,idxV_UsedVehicle)>0 ); %idxRow: �Ӿ�����item������, idxCol: ������
% %     for tmp = 1: max(idxCol)  %ÿ��used�����������ͳɱ�
% %         tmpN= [tmpN idxRow(idxCol==tmp)']; %�Ѿ����ŵ�item������
% %         tmpCost = max(idecM_cost2( idxV_UsedVehicle(tmp), idxRow(idxCol==tmp) )); %555 ���� �� ���� �������� = �������� %     tmpCost = max(idecM_cost( M_BCapacity(2, idxV_UsedVehicle(tmp)), idxRow(idxCol==tmp) )); %555 ����(��ֵ) �� ����(����) �������� = ��������
% %         V_cff(idxV_UsedVehicle(tmp)) = max(V_cff(idxV_UsedVehicle(tmp)), tmpCost); %���³������۷���
% %                         %         idxRow(idxCol==tmp)
% %                         %         idecV_size(idxRow(idxCol==tmp))
% %                         %         sum(idecV_size(idxRow(idxCol==tmp)))
% %                         %         V_rff(idxV_UsedVehicle(tmp))
% %                         %         sum(wff)
% %         V_rff(idxV_UsedVehicle(tmp)) = V_rff(idxV_UsedVehicle(tmp)) - sum(idecV_size(idxRow(idxCol==tmp))); % 555 ʣ������ = ԭ����-�Ѱ���item������
% %         if length(V_rff(idxV_UsedVehicle(tmp))) ~= length(sum(idecV_size(idxRow(idxCol==tmp)))) error('EEE'); end
% %     end
% % %     V_rff
% % %     V_cff
% % end
% % % V_rff(V_rff<0)
% % NoVehicleFlag = 0;       %- 1�����޳����� �����1 ����ǿ��н�
% % 
% % %%
% % for i=1:n, %��p����˳��������Ź�Ӧ��
% %     if ismember(i, tmpN)
% %         continue;
% %     end
% %     %% Step1 ��ɸѡ��i�ɷ����Used���ĳ��� : idxV_UsedVehicle
% %     idxV_UsedVehicle = find(sum(wff,1)>0); %1.1 ����used��������
% %     if any(idxV_UsedVehicle) %����ǿ�
% %         idxV_UsedVehicle = idxV_UsedVehicle( V_rff(idxV_UsedVehicle) >= idecV_size(i) );  %1.2 Լ��1: ����ʣ������>=i�ߴ�
% %         if any(idxV_UsedVehicle)
% %             idxV_UsedVehicle =  idxV_UsedVehicle( sum(wff(:,idxV_UsedVehicle)) <= jdecV_MaxItemNb_2(idxV_UsedVehicle) );    %1.3 Լ��2: ÿ��������items��<= ��Ӧ����������
% %                  if length(jdecV_MaxItemNb_2(idxV_UsedVehicle)) ~= length(sum(wff(:,idxV_UsedVehicle))) error('1.3'); end
% %             if any(idxV_UsedVehicle)
% %                 [idxRow,idxCol,~] = find( wff(:,idxV_UsedVehicle)>0 ); %idxRow: �Ӿ�����item������, idxCol: ������
% %                 tmpV = [];
% %                 for tmp = 1: max(idxCol)
% %                     if ~all( idecM_depot(i, idxRow(idxCol==tmp)) ==0 )  %�������ȫ����ͨ
% %                         tmpV = [tmpV tmp];
% %                     end
% %                 end
% %                 idxV_UsedVehicle(tmpV) = [];  %1.4 Լ��3: ȥ��������ͨ�ĳ���
% %             end
% %         end
% %     end
% %     %% Step 2: �Ƿ�UsedVehicle�ڴ��ڳ���
% %     if any(idxV_UsedVehicle) %�������: �Ǿ�ѡ����õ�best bin; �������ѡ��һ��bin        
% %         [~,k]=min(abs( V_cff(idxV_UsedVehicle) - idecM_cost2(idxV_UsedVehicle,i)' ) ); %% 555 ����Ҫ������(Best Fit) %         [~,k]=min(abs( V_cff(idxV_UsedVehicle) - idecM_cost(M_BCapacity(2,idxV_UsedVehicle),i)' ) ); %% 555 ����Ҫ������(Best Fit)
% %         k = idxV_UsedVehicle(k);                  %��k������ѡ��
% % %         disp('From Used Bin');
% %     else
% %         %�ҳ�new bin, ���³�
% %         idxV_NewVehicle = find(sum(wff,1)==0); %����: �ճ�����������Լ��
% %         idxV_NewVehicle = idxV_NewVehicle( V_rff(idxV_NewVehicle) >= idecV_size(i) );  %��ֹ�³����� < item�������
% %         if ~any(idxV_NewVehicle), 
% %             NoVehicleFlag = 1;
% %             warning('Not Enough Vehicles');
% %             break;
% %         end  %���Ͳ����ڴ˳���55555
% %         k = randsample(length(idxV_NewVehicle),1);
% %         k = idxV_NewVehicle(k);                                  %��k���³���ѡ��
% % %         disp('From New Bin');
% %     end
% %     
% %     %% Step 3:i����k���в�����    
% %     wff(i, k) = 1;    %- ����item(i)�����bin(idxCol)
% %     V_rff(k)= V_rff(k) - idecV_size(i); %- ��������
% %     V_cff(k) = max( idecM_cost2( k, i),  V_cff(k));  %- ���³������� %      V_cff(k) = max( idecM_cost( M_BCapacity(2,k), i ),  V_cff(k))  %- ���³�������
% %                if find(V_rff<0),        error('Capacity constraint is broken!');     end %-�ƻ�����Լ����Ӧ������
% % %                if M_BCapacity(2,k) > length(V_MaxItemNb), error('M_BCapacity constraint is broken!'); end    
% % end
% % 
% % % wff��ԭ
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
