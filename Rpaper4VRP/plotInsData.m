function plotInsData(thisInsData)
%PLOTINSDATA 此处显示有关此函数的摘要
%   画图Golden等单车型VRP的Instance
%
da = thisInsData;

x0 = da.Depot.Coord(1);
y0 = da.Depot.Coord(2);
x = da.Node.Coord(1:end,1);  %depot也作为一个聚类
y = da.Node.Coord(1:end,2);

% xc = inst.cluxy(2:end,1)';
% yc = inst.cluxy(2:end,2)';

xmin = min(x); ymin = min(y);
xmax = max(x); ymax = max(y);

nCluster = size(da.Cluster.Coord,1) + 1;  %聚类数(含depot)
Colors = hsv(nCluster);

u =[da.Depot.IdxCluster; da.Customer.IdxCluster];
 for i=1 : nCluster
    X=[x0 x(u==i)' x0];    %X只取聚类为i的点
    Y=[y0 y(u==i)' y0];
    Color=0.8*Colors(i,:);  %不同聚类不同颜色
    
    %画图：不同颜色展示不同聚类的点
    plot(X,Y,'o',...       % 'LineWidth',1,...
        'MarkerSize',5,...
        'Color',Color,...
        'MarkerFaceColor','white');
    hold on;
end

%画图：画Depot点
plot(x0,y0,'ks',...
    'LineWidth',1,...
    'MarkerSize',8,...
    'MarkerFaceColor','blue');
%画图：画聚类质心点
% plot(xc,yc,'ks',...
%     'LineWidth',1,...
%     'MarkerSize',0.05,...
%     'MarkerFaceColor','red');

hold off;
grid on;
axis equal;
xlim([xmin xmax]);
ylim([ymin ymax]);

% p = fig2plotly;
end
