

function [endmember_3,id_max_3,endmember_3_norm,id_max_3_norm,abundanceMap,abundanceMap_norm] = MultiSpecEMClassificationFun(CUBE,in,num_EM,selections)
if numel(selections)==1
    selections=1:size(in,2);
end

ENDMEMSIG = nfindr(CUBE, num_EM,'NumIterations',1000,'ReductionMethod','PCA');

for n=1:size(in,2)
    abundanceMap(n,:) = lsqnonneg(ENDMEMSIG,in(:,n))';
end

[maxval_3,id_max_3]=max(abundanceMap);
endmember_3=selections(id_max_3);
% figure;
% for n=1:num_EM
%     nexttile
%     [xx,srt]=sort(abundanceMap(:,n));
%     imagesc(in(:,srt))
% end

ENDMEMSIG = nfindr(normalize(CUBE,3,'range'), num_EM,'NumIterations',1000,'ReductionMethod','PCA');

in = normalize(in,'range');
for n=1:size(in,2)
    abundanceMap_norm(n,:) = lsqnonneg(ENDMEMSIG,in(:,n))';
end

[maxval_3_norm,id_max_3_norm]=max(abundanceMap_norm);
endmember_3_norm=selections(id_max_3_norm);
% for n=1:num_EM
%     nexttile
%     [xx,srt]=sort(abundanceMap(:,n));
%     imagesc(in(:,srt))
% end
end
