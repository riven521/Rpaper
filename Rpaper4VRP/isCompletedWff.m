function flag = isCompletedWff(wff)
%ISCOMPLETEWFF �ж�wff�Ƿ�����

% 1 ���ж���1 2 ��Ҫ ÿ��֮�Ͷ�Ϊ1  3 ÿ��֮��>=0 <=����������
if ~all(sum(wff))  % ÿ��֮�Ͷ�Ϊ1  % if ~all(sum(iwff)) && (sum(iwff(:)) ~=  nCustomer), 
    flag = 0;
else
    flag = 1;
end

end


