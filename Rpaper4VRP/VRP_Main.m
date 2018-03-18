%% �������ļ�
% ����ģ���ļ�����������Matlab������
%% function S = VRP_Main()
%VRP_Main ������
% ���
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
%  ��ʼ��
clc;  clear; close all;  format long g; format bank; %NOTE ����MATLAB CODE ֧��
rng('default')
rng(1); % NOTE �Ƿ�����ı�־

fileType = 'vrp';
Para = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
Para.ins = struct('PHO',0.1,'GAMMA',1.0); %PHO���Ⱦ�������GAMMA������ͨ��
% Para.alg = struct('NAME','BF');
Para.ga = struct('MAX_ITERATION', 10, 'NPOP', 16, 'HEUID', 1);
% ��ȡlib����: ��TempInstance�ж�ȡ��׺Ϊvrp���ļ���nested�ṹ��libDatas
LibDataArray = readLibData('TempInstance',fileType); %TODO �����޸Ĳ�����δ�����author����
% % printstruct(LibDataArray)
% ����lib����: ��libDates��ȡ���º���õ�nested�ṹ��insDatas
InsDataArray = updateLibData(LibDataArray,Para.ins,fileType);
% % printstruct(InsDataArray)
% load('5.mat');
      
clear fileType LibDataArray;
% printstruct(insDatas(1))

%%
% ÿ���ļ��������� ���cell array: solutions; ÿ��cell����һ��struct array����
for iDa = 1:length(InsDataArray)
    % ���InsDataArray���л�ͼ
%     plotInsData(InsDataArray(iDa)); 

    % ���InsDataArray��������
    SortedInsDataArray = getSortedInsData( InsDataArray (iDa) );
% %     printstruct(SortedInsDataArray)
    % ���SortedInsDataArray���� BF �㷨
    for iBF = 1:length(SortedInsDataArray)
        SolutionArray(iBF) = callHeuBF(SortedInsDataArray(iBF)); 
    end
% %     printstruct(SolutionArray)
    % ��ӡBF�����Ϣ
    for iSolu =1:length(SolutionArray), fprintf('BF solu = %d ', iSolu); printSolu(SolutionArray(iSolu)); end   %��ӡ��Ϣ 
    % ����BF�����Ϣ
%     for iSolu =1:length(SolutionArray), plotInsData(InsDataArray(iDa)); hold on; plotSolu(SolutionArray(iSolu)); end   %��ӡ��Ϣ 
    % BF�����Solution��
    clear iBF iSolu SortedInsDataArray;   clear SolutionArray;
    
    % ���InsDataArray���� GA �㷨
    SolutionMatrix = callHeuGA(InsDataArray,Para.ga);
    % ��ӡGA�����Ϣ   
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
end %����VRP_Main