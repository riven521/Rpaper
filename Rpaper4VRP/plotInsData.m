function plotInsData(thisInsData)
%PLOTINSDATA �˴���ʾ�йش˺�����ժҪ
%   ��ͼGolden�ȵ�����VRP��Instance
%
da = thisInsData;

x0 = da.Depot.Coord(1);
y0 = da.Depot.Coord(2);
x = da.Node.Coord(1:end,1);  %depotҲ��Ϊһ������
y = da.Node.Coord(1:end,2);

% xc = inst.cluxy(2:end,1)';
% yc = inst.cluxy(2:end,2)';

xmin = min(x); ymin = min(y);
xmax = max(x); ymax = max(y);

nCluster = size(da.Cluster.Coord,1) + 1;  %������(��depot)
Colors = hsv(nCluster);

u =[da.Depot.IdxCluster; da.Customer.IdxCluster];
 for i=1 : nCluster
    X=[x0 x(u==i)' x0];    %Xֻȡ����Ϊi�ĵ�
    Y=[y0 y(u==i)' y0];
    Color=0.8*Colors(i,:);  %��ͬ���಻ͬ��ɫ
    
    %��ͼ����ͬ��ɫչʾ��ͬ����ĵ�
    plot(X,Y,'o',...       % 'LineWidth',1,...
        'MarkerSize',5,...
        'Color',Color,...
        'MarkerFaceColor','white');
    hold on;
end

%��ͼ����Depot��
plot(x0,y0,'ks',...
    'LineWidth',1,...
    'MarkerSize',8,...
    'MarkerFaceColor','blue');
%��ͼ�����������ĵ�
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
