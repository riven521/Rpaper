function LibDataArray = readLibData(libFolder,fileType)
%readLibData ��ȡLib�����ı�����
% ����
%   libFolder            ��Ŷ�ȡ�ı���Ŀ¼��
%   fileType              ����ȡ�ļ�����
% ���
%   LibDataArray     Lib�ļ��ṹ��
%% function readLibData(libFolder,fileType)
% Nested Functions:
% # getPathNames
% # getLibDates
% Local Functions:
% # VRP_Input
% # CCVRP_Input

pathFileNameArray = getPathNames();
if isempty(pathFileNameArray), error('ERROR : δ��ȡ���κ������ļ���'); end

LibDataArray = getLibDates();
if isempty(LibDataArray), error('ERROR : δ��ȡ���κ������ļ�'); end

fprintf('Reading file with readLibData() ... done \n');

%% % Ƕ�� FUNCTIONS %
% ==========================================================================
% Ƕ�� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Nested function getPathNames()
% 
% 
    function pathFileNameArray = getPathNames()
        %��ȡlibFolderĿ¼�������ļ���
        % pathFileNames - cell array ����ļ���
        % ������ڲ�ͬ���͵��ļ����á�*����ȡ���У������ȡ�ض������ļ���'.'�����ļ����ͣ������á�.jpg��
        
        fileFolder=fullfile(strcat('.\',libFolder));
        dirOutput=dir(fullfile(fileFolder,strcat('*.',fileType)));
        fileNames={dirOutput.name}';
        pathFileNameArray = strcat('.\', libFolder, '\', fileNames);
        
        
    end

%% 2 Nested function getLibDates()
% 
% 
    function libDataArray = getLibDates()
        %������ȡĿ¼��ÿ���ļ�����
        % libDataArray - struct array ���lib���ݽṹ��

        if strcmp('vrp', fileType)
            
%         n = length(pathFileNameArray);
%             libDataArray(n,1) = struct('NAME',[],'BEST_KNOWN',[],'COMMENT',[],'VEHICLE_TYPE',[],...
%                 'DIMENSION',[],'EDGE_WEIGHT_FORMAT',[],'EDGE_WEIGHT_TYPE',[],'VEHICLE_Capacity',[],'VEHICLE_FixCost',[],...
%                 'VEHICLE_VariableCost',[],'VEHICLE_Number',[],'VEHICLE_FeePoint',[],'VEHICLE_MaxPoint',[],...
%                 'Node_Coord',[],'Node_Demand',[],'DEPOT_SECTION',[]);
            
            libDataArray = arrayfun(@VRP_Input,pathFileNameArray);
            
        elseif strcmp('ccvrp', fileType)
            
%         n = length(pathFileNameArray);
%             libDataArray(n,1) = struct('NAME',[],'TYPE',[],'COMMENT',[],'DIMENSION',[],...
%                 'CAPACITY',[],'EDGE_WEIGHT_TYPE',[],'NODE_COORD_TYPE',[],'Node_Coord',[],...
%                 'Demand',[],'Cluster',[],'DEPOT_SECTION',[]);

            libDataArray = arrayfun(@CCVRP_Input,pathFileNameArray);
            
        else
            error('Error. \n �ļ����ʹ���.')
        end
    end
end

%%  % �ֲ� FUNCTIONS %
% ==========================================================================
% �ֲ� FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 2.1 Local function 333 VRP_Input
%
%
function VRPStruct = VRP_Input(thisPathFileName)
%Open the VRP data file and test for errors
[fidtsplib, error_message] = fopen(thisPathFileName{1}, 'rt');
if fidtsplib == -1
    fprintf('error during fopen of data file (%s): %s\n', thisPathFileName, error_message);
else
    % Process data in the data file
    while true
        Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
        if any([~ischar(Line), strcmp('EOF', Line)])
            fprintf('Reading file %s ... done \n',thisPathFileName{1});
            break
        end
        [VRPOpt, LineRem] = strtok(Line, ':'); %�Ա���Line��:���Ž��л���
        if strcmp('NAME', VRPOpt) == 1 %����ַ����
            VRPStruct.NAME = strtok(LineRem(2:end));
        elseif strcmp('BEST_KNOWN', VRPOpt) == 1
            VRPStruct.BEST_KNOWN = strtok(LineRem(2:end));
        elseif strcmp('COMMENT', VRPOpt) == 1
            VRPStruct.COMMENT = strtok(LineRem(2:end));
        elseif strcmp('VEHICLE_TYPE', VRPOpt) == 1
            VRPStruct.VEHICLE_TYPE = sscanf(strtok(LineRem(2:end)), '%g');
        elseif strcmp('DIMENSION', VRPOpt) == 1
            VRPStruct.DIMENSION = sscanf(strtok(LineRem(2:end)), '%g');
        elseif strcmp('EDGE_WEIGHT_FORMAT', VRPOpt) == 1
            VRPStruct.EDGE_WEIGHT_FORMAT = strtok(LineRem(2:end));
        elseif strcmp('EDGE_WEIGHT_TYPE', VRPOpt) == 1
            VRPStruct.EDGE_WEIGHT_TYPE = strtok(LineRem(2:end));
        elseif strcmp('VEHICLE_SECTION', VRPOpt) == 1
            for IxData = 1:VRPStruct.VEHICLE_TYPE
                Line = fgetl(fidtsplib); %��ȡһ��
                Data = sscanf(Line, '%g'); %�ֲ���������
                VRPStruct.VEHICLE_Capacity(IxData, :) = Data(1)'; %����
                VRPStruct.VEHICLE_FixCost(IxData, :) = Data(2)'; %����
                VRPStruct.VEHICLE_VariableCost(IxData, :) = Data(3)'; %����
                VRPStruct.VEHICLE_Number(IxData, :) = Data(4)'; %����
                VRPStruct.VEHICLE_FeePoint(IxData, :) = Data(5)'; %����
                VRPStruct.VEHICLE_MaxPoint(IxData, :) = Data(6)'; %����
            end
            % �������������ϵ�һ��
        elseif strcmp('NODE_COORD_DEMAND_SECTION', VRPOpt) == 1
            for IxData = 1: VRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %��ȡһ��
                Data = sscanf(Line, '%g'); %�ֲ���������
                VRPStruct.Node_Coord(IxData, :) = Data(2:3)'; %��������Node_Coord
                VRPStruct.Node_Demand(IxData, :) = Data(4:end)'; %����������Demand
            end
        elseif strcmp('DEPOT_SECTION', VRPOpt) == 1
            IxData = 1;
            while IxData
                Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
                if any([~ischar(Line), strcmp('EOF', Line)])
                    break
                end
                Data = sscanf(Line, '%g'); %�ֲ���������
                VRPStruct.DEPOT_SECTION(IxData, :) = Data(1:end)';
                IxData = IxData + 1;
            end
        end
    end
    fclose(fidtsplib);
end %�ر��ļ�
end

%% 2.2 Local function 333 CCVRP_Input
%
%
function CCVRPStruct = CCVRP_Input(pathFileName)
% Open the CCVRP data file and test for errors
[fidtsplib, error_message] = fopen(pathFileName{1}, 'rt');
if fidtsplib == -1
    fprintf('error during fopen of data file (%s): %s\n', pathFileName, error_message);
else
    % Process data in the data file
    while true
        Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
        if any([~ischar(Line), strcmp('EOF', Line)])
            break
        end
        [VRPOpt, LineRem] = strtok(Line, ':'); %�Ա���Line��:���Ž��л���
        VRPOpt = strtrim(VRPOpt);
        if strcmp('NAME', VRPOpt) == 1 %����ַ����
            CCVRPStruct.NAME = strtok(LineRem(2:end));
        elseif strcmp('TYPE', VRPOpt) == 1
            CCVRPStruct.TYPE = strtok(LineRem(2:end));
        elseif strcmp('COMMENT', VRPOpt) == 1
            CCVRPStruct.COMMENT = strtok(LineRem(2:end));
        elseif strcmp('DIMENSION', VRPOpt) == 1
            CCVRPStruct.DIMENSION = sscanf(strtok(LineRem(2:end)), '%g');
        elseif strcmp('CAPACITY', VRPOpt) == 1
            CCVRPStruct.CAPACITY = sscanf(strtok(LineRem(2:end)), '%g');
        elseif strcmp('EDGE_WEIGHT_TYPE', VRPOpt) == 1
            CCVRPStruct.EDGE_WEIGHT_TYPE = strtok(LineRem(2:end));
        elseif strcmp('NODE_COORD_TYPE', VRPOpt) == 1
            CCVRPStruct.NODE_COORD_TYPE = strtok(LineRem(2:end));
        elseif strcmp('NODE_COORD_SECTION', VRPOpt) == 1
            for IxData = 1:CCVRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %��ȡһ��
                Data = sscanf(Line, '%g'); %�ֲ���������
                CCVRPStruct.Node_Coord(IxData, :) = Data(2:3)'; %��������Node_Coord
            end
        elseif strcmp('DEMAND_SECTION', VRPOpt) == 1
            for IxData = 1:CCVRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %��ȡһ��
                Data = sscanf(Line, '%g'); %�ֲ���������
                CCVRPStruct.Demand(IxData, :) = Data(2:2)'; %��������Demand
            end
        elseif strcmp('CLUSTER_SECTION', VRPOpt) == 1
            for IxData = 1:CCVRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %��ȡһ��
                Data = sscanf(Line, '%g'); %�ֲ���������
                CCVRPStruct.Cluster(IxData, :) = Data(2:2)'; %��������Cluster
            end
        elseif strcmp('DEPOT_SECTION', VRPOpt) == 1
            IxData = 1;
            while IxData
                Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
                if any([~ischar(Line), strcmp('EOF', Line)])
                    break
                end
                Data = sscanf(Line, '%g'); %�ֲ���������
                CCVRPStruct.DEPOT_SECTION(IxData, :) = Data(1:end)';
                IxData = IxData + 1;
            end
        end
    end
    disp(CCVRPStruct)
    fclose(fidtsplib);
end %�ر��ļ�
end