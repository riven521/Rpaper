function plotSolu(thisSolution)
%plotSolution 此处显示有关此函数的摘要
  

da = thisSolution.instance;
% fieldnames(da.instance)
wff = da.Veh_Cus.wff;

x0 = da.Depot.Coord(1);
y0 = da.Depot.Coord(2);
x = da.Customer.Coord(1:end,1);  %depot不作为点
y = da.Customer.Coord(1:end,2);

% xc = clu.xy(2:end,1)';
% yc = clu.xy(2:end,2)';

xmin = min(x);   ymin = min(y);
xmax = max(x);  ymax = max(y);

% Colors = hsv(clu.num);

%获取每条线路的点数cell
% maxLine  = da.Vehicle.Cost_nCustomer; %最大线路数,>=空线路数
maxLine = size(wff,1);   %最大线路数,>=空线路数
for iLine=1:maxLine    %Line1 - LineX - 最大车辆数     
     idx = find( wff(iLine,:)==1);
     L{iLine} = idx;
end

Colors=hsv(maxLine); %J is ROUTE number

%% DEL
%获取非空线路条数 - 供颜色用
% J = length(find(sum(RecGlobalWff,1)>0)); 

% [a,b]=find(RecGlobalWff==1);
% for i = 1 : length(unique(b))
%     if length(a(find(b==i)))==0, continue, end;
%     L{i} = a(find(b==i));
% end

%%
for j=1:maxLine
    
    if isempty(L{j})
         continue;
    end

    X=[x0 x(L{j})' x0];
    Y=[y0 y(L{j})' y0];
    
    Color=0.8*Colors(j,:);
    
    %画线路
    plot(X,Y,'-',... %           
        'Color',Color,...
        'LineWidth',2.0,... %        'MarkerSize',10,...
        'MarkerFaceColor','blue');
    hold on;
    
end

%画Depot点
plot(x0,y0,'ks',...
    'LineWidth',2,...
    'MarkerSize',2,...
    'MarkerFaceColor','yellow');

hold off;
grid on;
axis equal;

xlim([xmin xmax]);
ylim([ymin ymax]);

% p = fig2plotly;
end