function plotSolu(thisSolution)
%plotSolution �˴���ʾ�йش˺�����ժҪ
  

da = thisSolution.instance;
% fieldnames(da.instance)
wff = da.Veh_Cus.wff;

x0 = da.Depot.Coord(1);
y0 = da.Depot.Coord(2);
x = da.Customer.Coord(1:end,1);  %depot����Ϊ��
y = da.Customer.Coord(1:end,2);

% xc = clu.xy(2:end,1)';
% yc = clu.xy(2:end,2)';

xmin = min(x);   ymin = min(y);
xmax = max(x);  ymax = max(y);

% Colors = hsv(clu.num);

%��ȡÿ����·�ĵ���cell
% maxLine  = da.Vehicle.Cost_nCustomer; %�����·��,>=����·��
maxLine = size(wff,1);   %�����·��,>=����·��
for iLine=1:maxLine    %Line1 - LineX - �������     
     idx = find( wff(iLine,:)==1);
     L{iLine} = idx;
end

Colors=hsv(maxLine); %J is ROUTE number

%% DEL
%��ȡ�ǿ���·���� - ����ɫ��
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
    
    %����·
    plot(X,Y,'-',... %           
        'Color',Color,...
        'LineWidth',2.0,... %        'MarkerSize',10,...
        'MarkerFaceColor','blue');
    hold on;
    
end

%��Depot��
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