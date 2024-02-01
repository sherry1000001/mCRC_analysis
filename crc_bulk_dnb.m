data = importdata('data/bulk_data_0626.csv');
feature = data.textdata(2:end,1);

candi.alpha=[0.2,0.2];
[dnbset,maxCI,resModul,dnbstruct,dnbnet]=dnb(data.data,feature,timeIdx,timeIdx_con,candi,clust);
dnbset_2_2=dnbset;


fid = fopen('dnb_3w_1025_0811.txt','w');
for i = 1:size(dnbset,1)
    a = dnbset(i,:);
    a = cell2mat(a);
    fprintf(fid,'%s\n',a);
end
fclose(fid);




