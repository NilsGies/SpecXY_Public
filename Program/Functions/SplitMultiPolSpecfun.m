function [y,status]=SplitMultiPolSpecfun(in,n)
status='';

% if not(size(in,2)==n)
%     status='error wrong n value';
%     y={};
%     return
% end

y=cell(1,n);

for k=1:n
    y{k}=in(:,k:n:end);
end
status='OK';
end

%