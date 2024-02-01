% get modules by clustering the filterd candidate genes for each time point
% adjust the cutoff by the distribution of correlation (quantile crit)
% CREATED BY WANWEI ZHANG
function [modules,modIdx,cutoff,Z]=clusterModule_adj(candidate,crit,feature,data,timeIdx,candIdx)
if nargin<6
    candIdx=cell(size(candidate));
    for i=1:numel(candidate)
        if numel(candidate{i})>0
            candIdx{i}=zeros(size(feature,1),1)==1;
            for j=1:numel(candidate{i})
                candIdx{i}=candIdx{i} | strcmp(candidate{i}{j},feature(:,1));
            end
        end
    end
end

if islogical(candIdx{1})
    for i=1:numel(candIdx)
        candIdx{i}=find(candIdx{i});
    end
end

modules=cell(size(candidate));
modIdx=cell(size(candidate));
cutoff=zeros(size(candidate));
for i=1:numel(candidate)
    if numel(candidate{i})>0
%     Z=linkage(expr,'average','correlation');
    expr=data(candIdx{i},timeIdx{i,2});
    n=size(expr,1);
    m=size(expr,2);
    z=(expr-repmat(mean(expr,2),1,m))./repmat(std(expr,1,2)*sqrt(m),1,m);
    y=zeros(1,n*(n-1)/2);
    l=0;
    for k=1:n-1
        for j=k+1:n
            l=l+1;
            for h=1:m
                y(l)=y(l)+z(j,h)*z(k,h);
            end
            %y(l)=sum(z(j,:).*z(k,:));
        end
    end
    y=1-abs(y);
    y(y<0)=0;
    y_s=sort(y);
    l=round(numel(y)*crit);
    while isnan(y_s(l))
        l=l-1;
    end
    while y_s(l)<=0
        l=l+1;
    end
    cutoff(i)=y_s(l);
    %{
    distMat=1-abs(pcc_xy(expr,expr));
    distMat(distMat<0)=0;
    y=zeros(1,size(expr,1)*(size(expr,1)-1)/2);
    l=0;
    for j=1:size(expr,1)-1
        for k=j+1:size(expr,1)
            l=l+1;
            y(l)=distMat(k,j);
        end
    end
    %}
    %{
    y=pdist(expr,@dist_absPCC);
    %}
    Z=linkage(y,'average');
%     Z=linkage(y,'single');
    T=cluster(Z,'cutoff',cutoff(i),'criterion','distance');
%     T=cluster(Z,'cutoff',cutoff);
    modules{i}=cell(max(T),1);
    modIdx{i}=cell(max(T),1);
    for j=1:max(T)
        modules{i}{j}=candidate{i}(T==j);
        modIdx{i}{j}=candIdx{i}(T==j);
    end
%     temp=cellfun(@numel,modules{1});
    else
        warning(['No Candidates found at timepoint: ',num2str(i)])
    end
end
