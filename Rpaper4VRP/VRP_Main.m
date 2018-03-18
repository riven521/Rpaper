%% 主函数文件
% 今后的模板文件，用于美观Matlab发布。
%% function S = VRP_Main()
%VRP_Main 主函数
% 输出
%   S     
%% function VRP_Main(libFolder,fileType)
% Call Functions:
% # readLibData
% # updateLibData
% # getSortedInsData
% # callHeuBF
% # printSolu
%%
function S = VRP_Main()
%  初始化
clc;  clear; close all;  format long g; format bank; %NOTE 不被MATLAB CODE 支持
rng('default')
rng(1); % NOTE 是否随机的标志

fileType = 'vrp';
Para = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
Para.ins = struct('PHO',0.1,'GAMMA',1.0); %PHO正比聚类数；GAMMA正比连通性
% Para.alg = struct('NAME','BF');
Para.ga = struct('MAX_ITERATION', 10, 'NPOP', 16, 'HEUID', 1);
% 获取lib数据: 从TempInstance中读取后缀为vrp的文件到nested结构体libDatas
LibDataArray = readLibData('TempInstance',fileType); %TODO 后期修改参数二未以提出author命名
% % printstruct(LibDataArray)
% 更新lib数据: 从libDates获取更新后可用的nested结构体insDatas
InsDataArray = updateLibData(LibDataArray,Para.ins,fileType);
% % printstruct(InsDataArray)
% load('5.mat');
      
clear fileType LibDataArray;
% printstruct(insDatas(1))

%%
% 每个文件单独运行 输出cell array: solutions; 每个cell包含一个struct array数组
for iDa = 1:length(InsDataArray)
    % 针对InsDataArray进行画图
%     plotInsData(InsDataArray(iDa)); 

    % 针对InsDataArray进行排序
    SortedInsDataArray = getSortedInsData( InsDataArray (iDa) );
% %     printstruct(SortedInsDataArray)
    % 针对SortedInsDataArray调用 BF 算法
    for iBF = 1:length(SortedInsDataArray)
        SolutionArray(iBF) = callHeuBF(SortedInsDataArray(iBF)); 
    end
% %     printstruct(SolutionArray)
    % 打印BF解的信息
    for iSolu =1:length(SolutionArray), fprintf('BF solu = %d ', iSolu); printSolu(SolutionArray(iSolu)); end   %打印信息 
    % 画出BF解的信息
%     for iSolu =1:length(SolutionArray), plotInsData(InsDataArray(iDa)); hold on; plotSolu(SolutionArray(iSolu)); end   %打印信息 
    % BF解放入Solution中
    clear iBF iSolu SortedInsDataArray;   clear SolutionArray;
    
    % 针对InsDataArray调用 GA 算法
    SolutionMatrix = callHeuGA(InsDataArray,Para.ga);
    % 打印GA解的信息   
    for iIter = 1:size(SolutionMatrix,1)
        [iterBestArray(iIter),~] = findMinIdxOfArray(SolutionMatrix(iIter,:));
%   fprintf('GA Generation = %d ', iIter); printSolu(iterBestArray(iIter));
    end    
    [globalBestSolu,globalBestSoluIdx] = findMinIdxOfArray(iterBestArray);
    fprintf('GA Generation Idx = %d ', globalBestSoluIdx);
    printSolu(globalBestSolu);
        
    clear iIter SolutionMatrix globalBestSoluIdx iterBestArray;   clear globalBestSolu;
    
%     SOLUTIONS{iDa} = globalBestSolu;

end 
clear iDa Para InsDataArray;
 
%     S = SOLUTIONS{1}.instance;
%     save('s.mat','S');
end %结束VRP_Main