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
par = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
par.ins = struct('pho',{0.2},'gamma',{0.2});

% ��ȡlib����: ��TempInstance�ж�ȡ��׺Ϊvrp���ļ���nested�ṹ��libDatas
libDatas = readLibData('TempInstance',fileType); %TODO �����޸Ĳ�����δ�����author����
    % printstruct(libDatas)

% ����lib����: ��libDates��ȡ���º���õ�nested�ṹ��insDatas
insDatas = updateLibDatas(libDatas,par.ins,fileType);
    % printstruct(insDatas(1));

% ��������
S = insDatas(1);

end



