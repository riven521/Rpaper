function flag = isCompletedWff(wff)
%ISCOMPLETEWFF 判断wff是否完整

% 1 所有都是1 2 重要 每列之和都为1  3 每行之和>=0 <=车辆最大点数
if ~all(sum(wff))  % 每列之和都为1  % if ~all(sum(iwff)) && (sum(iwff(:)) ~=  nCustomer), 
    flag = 0;
else
    flag = 1;
end

end


