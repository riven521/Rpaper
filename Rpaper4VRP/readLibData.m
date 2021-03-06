function LibDataArray = readLibData(libFolder,fileType)
%readLibData 读取Lib类型文本数据
% 输入
%   libFolder            存放读取文本的目录名
%   fileType              待读取文件类型
% 输出
%   LibDataArray     Lib文件结构体
%% function readLibData(libFolder,fileType)
% Nested Functions:
% # getPathNames
% # getLibDates
% Local Functions:
% # VRP_Input
% # CCVRP_Input

pathFileNameArray = getPathNames();
if isempty(pathFileNameArray), error('ERROR : 未获取到任何数据文件名'); end

LibDataArray = getLibDates();
if isempty(LibDataArray), error('ERROR : 未读取到任何数据文件'); end

fprintf('Reading file with readLibData() ... done \n');

%% % 嵌套 FUNCTIONS %
% ==========================================================================
% 嵌套 FUNCTIONS  ------------------------------------------------------------=
% ==========================================================================

%% 1 Nested function getPathNames()
% 
% 
    function pathFileNameArray = getPathNames()
        %获取libFolder目录下所有文件名
        % pathFileNames - cell array 存放文件名
        % 如果存在不同类型的文件，用‘*’读取所有，如果读取特定类型文件，'.'加上文件类型，例如用‘.jpg’
        
        fileFolder=fullfile(strcat('.\',libFolder));
        dirOutput=dir(fullfile(fileFolder,strcat('*.',fileType)));
        fileNames={dirOutput.name}';
        pathFileNameArray = strcat('.\', libFolder, '\', fileNames);
        
        
    end

%% 2 Nested function getLibDates()
% 
% 
    function libDataArray = getLibDates()
        %迭代获取目录下每个文件数据
        % libDataArray - struct array 存放lib数据结构体

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
            error('Error. \n 文件类型错误.')
        end
    end
end

%%  % 局部 FUNCTIONS %
% ==========================================================================
% 局部 FUNCTIONS  ------------------------------------------------------------=
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
        [VRPOpt, LineRem] = strtok(Line, ':'); %对变量Line按:符号进行划分
        if strcmp('NAME', VRPOpt) == 1 %如果字符相等
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
                Line = fgetl(fidtsplib); %读取一行
                Data = sscanf(Line, '%g'); %分布读入数据
                VRPStruct.VEHICLE_Capacity(IxData, :) = Data(1)'; %读入
                VRPStruct.VEHICLE_FixCost(IxData, :) = Data(2)'; %读入
                VRPStruct.VEHICLE_VariableCost(IxData, :) = Data(3)'; %读入
                VRPStruct.VEHICLE_Number(IxData, :) = Data(4)'; %读入
                VRPStruct.VEHICLE_FeePoint(IxData, :) = Data(5)'; %读入
                VRPStruct.VEHICLE_MaxPoint(IxData, :) = Data(6)'; %读入
            end
            % 坐标和需求量混合到一起
        elseif strcmp('NODE_COORD_DEMAND_SECTION', VRPOpt) == 1
            for IxData = 1: VRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %读取一行
                Data = sscanf(Line, '%g'); %分布读入数据
                VRPStruct.Node_Coord(IxData, :) = Data(2:3)'; %读入坐标Node_Coord
                VRPStruct.Node_Demand(IxData, :) = Data(4:end)'; %读入需求率Demand
            end
        elseif strcmp('DEPOT_SECTION', VRPOpt) == 1
            IxData = 1;
            while IxData
                Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
                if any([~ischar(Line), strcmp('EOF', Line)])
                    break
                end
                Data = sscanf(Line, '%g'); %分布读入数据
                VRPStruct.DEPOT_SECTION(IxData, :) = Data(1:end)';
                IxData = IxData + 1;
            end
        end
    end
    fclose(fidtsplib);
end %关闭文件
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
        [VRPOpt, LineRem] = strtok(Line, ':'); %对变量Line按:符号进行划分
        VRPOpt = strtrim(VRPOpt);
        if strcmp('NAME', VRPOpt) == 1 %如果字符相等
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
                Line = fgetl(fidtsplib); %读取一行
                Data = sscanf(Line, '%g'); %分布读入数据
                CCVRPStruct.Node_Coord(IxData, :) = Data(2:3)'; %读入坐标Node_Coord
            end
        elseif strcmp('DEMAND_SECTION', VRPOpt) == 1
            for IxData = 1:CCVRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %读取一行
                Data = sscanf(Line, '%g'); %分布读入数据
                CCVRPStruct.Demand(IxData, :) = Data(2:2)'; %读入坐标Demand
            end
        elseif strcmp('CLUSTER_SECTION', VRPOpt) == 1
            for IxData = 1:CCVRPStruct.DIMENSION
                Line = fgetl(fidtsplib); %读取一行
                Data = sscanf(Line, '%g'); %分布读入数据
                CCVRPStruct.Cluster(IxData, :) = Data(2:2)'; %读入坐标Cluster
            end
        elseif strcmp('DEPOT_SECTION', VRPOpt) == 1
            IxData = 1;
            while IxData
                Line = fgetl(fidtsplib);             %%% fprintf(fid,'%c\n',Line);
                if any([~ischar(Line), strcmp('EOF', Line)])
                    break
                end
                Data = sscanf(Line, '%g'); %分布读入数据
                CCVRPStruct.DEPOT_SECTION(IxData, :) = Data(1:end)';
                IxData = IxData + 1;
            end
        end
    end
    disp(CCVRPStruct)
    fclose(fidtsplib);
end %关闭文件
end