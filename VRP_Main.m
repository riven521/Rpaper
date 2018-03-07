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
par = struct('file',[],'ins',[]); % par.file = struct('name',{'~'},'load',{0},'gene',{0});
par.ins = struct('pho',{0.2},'gamma',{0.2});

% 获取lib数据: 从TempInstance中读取后缀为vrp的文件到nested结构体libDatas
libDatas = readLibData('TempInstance',fileType); %TODO 后期修改参数二未以提出author命名
    % printstruct(libDatas)

% 更新lib数据: 从libDates获取更新后可用的nested结构体insDatas
insDatas = updateLibDatas(libDatas,par.ins,fileType);
    % printstruct(insDatas(1));

% 返回数据
S = insDatas(1);

end



