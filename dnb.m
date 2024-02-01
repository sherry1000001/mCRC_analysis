% function [dnbset,maxCI,resModul,modules,candidate,CI,outIdx]=dnb(data,feature,timeIdx,...
function [dnbset,maxCI,resModul,dnbstruct,dnbnet]=dnb(data,feature,timeIdx,...
    timeIdx_con,candi,clust,out,display,scale,dnbIdx,alpha)
    % Detecting Dynamic Network Biomarker
    % data: feature by sample matrix
    % feature: names of gene or proteins or others
    % timeIdx: indices of time points [{'Time1'},{[1,2,...]};
    %                                  {'Time2'},{[3,4,...]};
    %                                          ...
    %                                  {'TimeN'},{[8,9,...]}]
    % TPC: time slide window (number of time point combined as one)
    % clust: arguments for hierarchical clustering
    % clust.crit: criterion for selecting the cutoffs in cutting modules from
    % clust.adjust: 1 clustering with adjust cutoff
    %               0 clustering with fix dist as cutoff
    % clust.msmin: minimum module size
    % candi: arguments for candidate selection
    % candi.alpha: criterion for selecting candidates
    %              if control samples exist, alpha is the significance level of statistical test,
    %              else alpha is the ratio of features with top SD of CV
    %              hirarchical cluster tree, based on distribution of distance
    % candi.method: method for selecting candidates
    %               0 top alpha(1) ratio of SD, top alpha(2) ratio of CV.
    %               1 differential SD (fold change > alpha(1)),
    %                 differential expression (Pvalue<alpha(2))
    %               2 selecting top alpha(1) ratio SD and differential
    %                 expression (pval<alpha(2))
    %               3 using SD>alpha.
    %               4 SDcurrent/mean(SDremain)>candi.alpha
    % out.crit: [0...1,0...1],criterion for selecting the outside genes.High correlate
    %          with modules in previous time point and low correlated in
    %          current time. if control samples exist, select genes that
    %          are highly correlated with given module in current time.
    % timeIdx_con: indices of control time points.
    % display: [0,1] display the plot
    %           2 plot to specified file
    % scale: 2  scaling data along feature(row)
    %        1  scaling data along sample(column)
    %        0  no scaling
    %        3  block normalization
    %        4  '2' + '3'
    % timepoint combination (time window)
    %
    %
    % DEMO
    % data=randn(3000,100);
    % feature=strcat('feature',cellfun(@num2str,num2cell(1:3000)','unif',0));
    % timeIdx={'t1',(1:10)';
    %     't2',(11:20)';
    %     't3',(21:30)';
    %     't4',(31:40)';
    %     't5',(41:50)'};
    % timeIdx_con={'t1',(51:60)';
    %     't2',(61:70)';
    %     't3',(71:80)';
    %     't4',(81:90)';
    %     't5',(91:100)'};
    % dnbset=dnb(data,feature,timeIdx,timeIdx_con);
    % CREATED BY WANWEI ZHANG

if nargin<11
    alpha=1;
end
if nargin<10
    dnbIdx=[];
end
if nargin<9
    scale=2;
end
if nargin<8
    display=1;
end
if nargin<7
    if ~isempty(timeIdx_con)
        out.method=0;
        out.crit=[1,1];
    else
        out.method=1;
        out.crit=[1,1];
    end
end
if nargin<6
    clust.adjust=0;
    clust.crit=0.2;
    clust.msmin=3;
end
if nargin<5
    if ~isempty(timeIdx_con)
        candi.method=2;
        candi.alpha=[0.2,0.05];
    else
        candi.method=0;
        candi.alpha=[0.2,0.2];
    end
end
if nargin<4
    timeIdx_con=[];
end

TPC=1;

timeIdx_c=timeIdx;
for i=1:size(timeIdx,1)-TPC+1
    for j=1:TPC-1
        timeIdx_c{i,2}=[timeIdx_c{i,2},timeIdx{i+j,2}];
    end
end
for i=size(timeIdx,1)-TPC+2:size(timeIdx,1)
    for j=i+1:size(timeIdx,1)
        timeIdx_c{i,2}=[timeIdx_c{i,2},timeIdx{j,2}];
    end
end

if ~isempty(timeIdx_con)
timeIdx_con_c=timeIdx_con;
for i=1:size(timeIdx_con,1)-TPC+1
    for j=1:TPC-1
        timeIdx_con_c{i,2}=[timeIdx_con_c{i,2},timeIdx_con{i+j,2}];
    end
end
for i=size(timeIdx_con,1)-TPC+2:size(timeIdx_con,1)
    for j=i+1:size(timeIdx_con,1)
        timeIdx_con_c{i,2}=[timeIdx_con_c{i,2},timeIdx_con{j,2}];
    end
end
end

idx=cell2mat(timeIdx(:,2));
if ~isempty(timeIdx_con)
    idx=[idx;cell2mat(timeIdx_con(:,2))];
end
idx=unique(idx);
if scale>0 && scale<3
    data(:,idx)=zscore(data(:,idx),0,scale);
end

if scale==3 % block normalization
    for i=1:size(timeIdx,1)
        temp=data(:,timeIdx{i,2});
        data(:,timeIdx{i,2})=(data(:,timeIdx{i,2})-mean(temp(:)))/std(temp(:));
    end
    if ~isempty(timeIdx_con)
        for i=1:size(timeIdx_con,1)
            temp=data(:,timeIdx_con{i,2});
            data(:,timeIdx_con{i,2})=(data(:,timeIdx_con{i,2})-mean(temp(:)))/std(temp(:));
        end
    end
end

if scale==4 % block normalization + feature normalization
    % feature normalization
    data(:,idx)=zscore(data(:,idx),0,2);
    % block normalization 
    for i=1:size(timeIdx,1)
        temp=data(:,timeIdx{i,2});
        data(:,timeIdx{i,2})=(data(:,timeIdx{i,2})-mean(temp(:)))/std(temp(:));
    end
    if ~isempty(timeIdx_con)
        for i=1:size(timeIdx_con,1)
            temp=data(:,timeIdx_con{i,2});
            data(:,timeIdx_con{i,2})=(data(:,timeIdx_con{i,2})-mean(temp(:)))/std(temp(:));
        end
    end
end

%%
% get candidate genes (candidate)
% get sample std (ssd)
ssd=zeros(size(data,1),size(timeIdx_c,1));
for i=1:size(timeIdx_c,1)
    ssd(:,i)=std(data(:,timeIdx_c{i,2}),0,2);
end

% get control sample std 
if ~isempty(timeIdx_con)
    ssd_con=zeros(size(data,1),size(timeIdx_con_c,1));
    for i=1:size(timeIdx_con_c,1)
        ssd_con(:,i)=std(data(:,timeIdx_con_c{i,2}),0,2);
    end
end

% top alpha(1) ratio of SD, top alpha(2) ratio of CV
if candi.method==0
    %CV
    expr_mean=zeros(size(data,1),size(timeIdx_c,1));
    for i=1:size(timeIdx_c,1)
        expr_mean(:,i)=mean(data(:,timeIdx_c{i,2}),2);
    end
    cv=ssd./expr_mean;
    
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        %CV
        cv_tmp=sort(cv(:,i),'descend');
        temp=floor(candi.alpha(2)*numel(cv_tmp));
        if temp==0
            temp=nan;% threshold for no cv selected
        else
            temp=cv_tmp(temp);
        end
        candIdx{i}=cv(:,i)>=temp;
        
        temp=sort(ssd(:,i),'descend');
        temp=temp(floor(candi.alpha(1)*numel(temp)));
        candIdx{i}=candIdx{i}|ssd(:,i)>=temp;
%         candIdx{i}=ssd(:,i)>=temp;
        %{
        figure
        hold
        plot(expr_mean(:,i),ssd(:,i),'.','markersize',0.5)
        plot(expr_mean(candIdx{i},i),ssd(candIdx{i},i),'.r','markersize',0.5)
        %}
    end
end

% selecting differential SD and differential expression
if candi.method==1
    if isempty(timeIdx_con)
        error('Control data must exist! (Differential SD or expression)')
    end
    foldchange=ssd./ssd_con;
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        pval=mattest(data(:,timeIdx_c{i,2}),data(:,timeIdx_con_c{i,2}));
        candIdx{i}=foldchange(:,i)>candi.alpha(1) | pval<candi.alpha(2);
    end
end

% selecting top alpha(1) ratio SD and differential expression (pval<alpha(2))
if candi.method==2
    if isempty(timeIdx_con)
        error('Control data must exist! (Differential expression)')
    end
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        temp=sort(ssd(:,i),'descend');
        temp=temp(floor(candi.alpha(1)*numel(temp)));
        pval=mattest(data(:,timeIdx_c{i,2}),data(:,timeIdx_con_c{i,2}));
        candIdx{i}=ssd(:,i)>=temp | pval<candi.alpha(2);
    end
end
    
% SD>alpha.
if candi.method==3
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        candIdx{i}=ssd(:,i)>=candi.alpha;
    end
end

% SDcurrent/mean(SDremain)>candi.alpha
if candi.method==4
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        temp=ssd(:,i)./mean(ssd(:,[1:i-1,i+1:end]),2);
        candIdx{i}=temp>candi.alpha;
    end
end

% variance of sd is large
if candi.method==5
    sdsd=std(ssd,0,2);
    [temp,idx]=sort(sdsd,'descend');
    sdIdx=zeros(size(idx));
    sdIdx(idx(1:round(candi.alpha(1)*numel(idx))))=1;
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        pval=mattest(data(:,timeIdx_c{i,2}),data(:,timeIdx_con_c{i,2}));
        candIdx{i}=pval<candi.alpha(2)|sdIdx==1;
    end
end

% SDcurrent significantly higher than SDremain:
% ttest(SDremain,SDcurrent,candi.alpha,'left')==1
if candi.method==6
    candIdx=cell(size(timeIdx_c,1),1);
    for i=1:size(timeIdx_c,1)
        h=myRowfun(@ttest,[],ssd(:,[1:i-1,i+1:end]),ssd(:,i),candi.alpha,'left');
        candIdx{i}=h==1;
    end
end

candidate=cell(size(timeIdx_c,1),1);
for i=1:size(timeIdx_c,1)
    candidate{i}=feature(candIdx{i});
end

%% get modules
if clust.adjust==1
    [raw_mod,raw_modIdx]=clusterModule_adj(candidate,clust.crit,feature,data,timeIdx_c,candIdx);
else
    [raw_mod,raw_modIdx]=clusterModule(candidate,clust.crit,feature,data,timeIdx_c,candIdx);
end
%{
% check largest and smallest module size
msize=zeros(numel(raw_mod),2);
for i=1:numel(raw_mod)
    temp=cellfun(@numel,raw_mod{i});
    msize(i,1)=min(temp);
    msize(i,2)=max(temp);
%     figure
%     plot(sort(temp));
end
%}
% take modules with >= msmin genes
maxsize=zeros(size(raw_mod));
for i=1:numel(raw_mod)
    temp=cellfun(@numel,raw_mod{i});
    maxsize(i)=max(temp);
end
if min(maxsize)<clust.msmin
    warning(['No module of size ',num2str(clust.msmin),' at some timepoints. Lower the size criterion: ',num2str(min(maxsize))]);
    clust.msmin=min(maxsize);
end

modules=raw_mod;
modIdx=raw_modIdx;
for i=1:numel(raw_mod)
    temp=cellfun(@numel,raw_mod{i});
    modules{i}(temp<clust.msmin)=[];
    modIdx{i}(temp<clust.msmin)=[];
    if isempty(modules{i})
%         fprintf(['Error:No module at timepoint:',num2str(i),'\n']);
        error(['No module at timepoint:',num2str(i)]);
    end
end

%% get outIdx
outIdx=cell(size(modIdx));

% select genes that are highly correlated in control samples but low
% correlated in case samples as outside member
if out.method==0
    if isempty(timeIdx_con)
        error('Current out-member selection method is not suitable: no control sample!');
    end
    for t=1:size(timeIdx_c,1)
        outIdx{t}=cell(size(modIdx{t}));
        for i=1:numel(modIdx{t})
            outIdx_tmp=setdiff((1:size(data,1))',modIdx{t}{i});

            PCCcon=pcc_xy(data(modIdx{t}{i},timeIdx_con_c{t,2}),...
                data(outIdx_tmp,timeIdx_con_c{t,2}));
            PCCcon(isnan(PCCcon))=0;
            PCCcase=pcc_xy(data(modIdx{t}{i},timeIdx_c{t,2}),...
                data(outIdx_tmp,timeIdx_c{t,2}));
            PCCcase(isnan(PCCcase))=0;
        
            [temp,idx]=sort(mean(abs(PCCcon),1),'descend');
            idx(isnan(temp))=[];
            outIdx{t}{i}=zeros(size(temp'));
            outIdx{t}{i}(idx(1:round(out.crit(1)*numel(idx))))=1;
            [temp,idx]=sort(mean(abs(PCCcase),1),'ascend');
            idx(isnan(temp))=[];
            outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))=...
                outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))+1;
            outIdx{t}{i}=outIdx_tmp(outIdx{t}{i}==2);
        
            if isempty(outIdx{t}{i})
                error(['a module in the ',num2str(t),' timepoint has no outsiders']);
            end
        end
    end
end

% select genes that are highly correlated in all case timepoints but low
% correlated in current case timepoint as outside member
if out.method==1
    throughTime=cell2mat(timeIdx(:,2));
    for t=1:size(timeIdx_c,1)
        outIdx{t}=cell(size(modIdx{t}));
        for i=1:numel(modIdx{t})
            outIdx_tmp=setdiff((1:size(data,1))',modIdx{t}{i});
            PCCthrough=pcc_xy(data(modIdx{t}{i},throughTime),...
                data(outIdx_tmp,throughTime));
            PCCthrough(isnan(PCCthrough))=0;
            PCCcurrent=pcc_xy(data(modIdx{t}{i},timeIdx_c{t,2}),...
                data(outIdx_tmp,timeIdx_c{t,2}));
            PCCcurrent(isnan(PCCcurrent))=0;
            
            [temp,idx]=sort(mean(abs(PCCthrough),1),'descend');
            idx(isnan(temp))=[];
            outIdx{t}{i}=zeros(size(temp'));
            outIdx{t}{i}(idx(1:round(out.crit(1)*numel(idx))))=1;
            [temp,idx]=sort(mean(abs(PCCcurrent),1),'ascend');
            idx(isnan(temp))=[];
            outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))=...
                outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))+1;
            outIdx{t}{i}=outIdx_tmp(outIdx{t}{i}==2);

            if isempty(outIdx{t}{i})
                error(['a module in the ',num2str(t),' timepoint has no outsiders']);
            end
        end
    end
end

% select those genes that are highly correlated (top 1-75% pcc) with the module during previous timepoints
% but loss the correlation (bottom 25% pcc) in current timepoint
if out.method==2
    for t=1:size(timeIdx_c,1)
        outIdx{t}=cell(size(modIdx{t}));
        for i=1:numel(modIdx{t})
            outIdx_tmp=setdiff((1:size(data,1))',modIdx{t}{i});
            
            if t==1
                prevTime=timeIdx_c{2,2};
            else
                prevTime=timeIdx_c{t-1,2};
            end
            %{
            prevTime=[];
            for j=1:t
                prevTime=[prevTime;timeIdx{j,2}];
            end
            %}
            PCCprev=pcc_xy(data(modIdx{t}{i},prevTime),...
                data(outIdx_tmp,prevTime));
            PCCprev(isnan(PCCprev))=0;
            PCCcurrent=pcc_xy(data(modIdx{t}{i},timeIdx_c{t,2}),...
                data(outIdx_tmp,timeIdx_c{t,2}));
            PCCcurrent(isnan(PCCcurrent))=0;
            
            [temp,idx]=sort(mean(abs(PCCprev),1),'descend');
            idx(isnan(temp))=[];
            outIdx{t}{i}=zeros(size(temp'));
            outIdx{t}{i}(idx(1:round(out.crit(1)*numel(idx))))=1;
            [temp,idx]=sort(mean(abs(PCCcurrent),1),'ascend');
            idx(isnan(temp))=[];
            outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))=...
                outIdx{t}{i}(idx(1:round(out.crit(2)*numel(idx))))+1;
            outIdx{t}{i}=outIdx_tmp(outIdx{t}{i}==2);
            %{
            outIdx{t}{i}=find(mean(abs(PCCprev),1)>out.crit(1) & ...
                mean(abs(PCCcurrent),1)<out.crit(2))';
            %}

            if isempty(outIdx{t}{i})
                error(['a module in the ',num2str(t),' timepoint has no outsiders']);
            end
        end
    end
end


%% get CI
% tic;
%{
alpha=1/2;
%}

CI=cell(size(modIdx));
for t=1:size(timeIdx_c,1)
    CI{t}=zeros(numel(modIdx{t}),4);
    for i=1:numel(modIdx{t})
        
        midx=zeros(size(data,1),1)>0;
        for j=1:numel(modules{t}{i})
            midx=midx | strcmp(modules{t}{i}{j},feature);
        end
        expr_m=data(midx,timeIdx_c{t,2});
        %{
        expr_m=data(modIdx{t}{i},timeIdx_c{t,2});
        %}
        
        % SD_in
        CI{t}(i,2)=sum(std(expr_m,0,2))/(size(expr_m,1)^alpha);
        
        % PCC_in
        temp=abs(corrcoef(expr_m'));
        temp(isnan(temp))=0;
        CI{t}(i,3)=sum(sum(temp,2)-diag(temp))/(size(expr_m,1)*(size(expr_m,1)-1));
        
        % PCC_out
        expr_o=data(outIdx{t}{i},timeIdx_c{t,2});
        temp=std(expr_o,0,2);
        expr_o(temp==0,:)=[];
        temp=pcc_xy(expr_m,expr_o);
        CI{t}(i,4)=mean(abs(temp(~isnan(temp))));
        
        % Critical Indicator
        CI{t}(i,1)=sqrt(size(expr_m,1))*CI{t}(i,2)*CI{t}(i,3)/CI{t}(i,4);
    end
end

% get maximum Critical Indicators in each timepoint
maxCI=nan*ones(size(CI,1),5);
for i=1:numel(CI)
    if ~isempty(CI{i})
    [temp,maxCI(i,5)]=max(abs(CI{i}(:,1)));
    for j=1:4
        maxCI(i,j)=CI{i}(maxCI(i,5),j);
    end
    end
end
%{
%}
% get responsive module
resModul=cell(size(maxCI,1),1);
for i=1:size(maxCI,1)
    if ~isempty(modules{i})
    resModul{i}=modules{i}{maxCI(i,5)};
    end
end
if isempty(dnbIdx)
    [temp,dnbIdx]=max(abs(maxCI(:,1)));
end

%% DNB behavior
dnbset=resModul{dnbIdx};
idx=zeros(size(feature))>0;
for i=1:numel(dnbset)
    idx=idx | strcmp(dnbset{i},feature);
end
idx=find(idx);
%{
idx=modIdx{dnbIdx}{maxCI(dnbIdx,5)};
%}
oidx=outIdx{dnbIdx}{maxCI(dnbIdx,5)};

dnbCI=zeros(size(timeIdx_c,1),4);
for i=1:size(timeIdx_c,1)
    expr_m=data(idx,timeIdx_c{i,2});
    temp=std(expr_m,0,2);
    expr_m(temp==0,:)=[];
    
    dnbCI(i,2)=sum(std(expr_m,[],2))/(numel(idx)^alpha);
    
    temp=abs(corr(expr_m'));
    dnbCI(i,3)=(sum(temp(~isnan(temp)))-size(temp,1))/(size(temp,1)*(size(temp,1)-1));
    
    expr_o=data(oidx,timeIdx_c{i,2});
    temp=std(expr_o,0,2);
    expr_o(temp==0,:)=[];
    
    temp=pcc_xy(expr_m,expr_o);
    dnbCI(i,4)=mean(abs(temp(~isnan(temp))));
    
    dnbCI(i,1)=dnbCI(i,2)*dnbCI(i,3)/dnbCI(i,4);
end
if ~isempty(timeIdx_con)
    dnbCI_con=zeros(size(timeIdx_con_c,1),4);
    for i=1:size(timeIdx_con_c,1)
        expr_m=data(idx,timeIdx_con_c{i,2});
        temp=std(expr_m,0,2);
        expr_m(temp==0,:)=[];
        
        dnbCI_con(i,2)=mean(std(expr_m,[],2));
        
        temp=abs(corr(expr_m'));
        dnbCI_con(i,3)=(sum(temp(~isnan(temp)))-size(temp,1))/(size(temp,1)*(size(temp,1)-1));
        
        expr_o=data(oidx,timeIdx_con_c{i,2});
        temp=std(expr_o,0,2);
        expr_o(temp==0,:)=[];
        
        temp=pcc_xy(expr_m,expr_o);
        dnbCI_con(i,4)=mean(abs(temp(~isnan(temp))));
        
        dnbCI_con(i,1)=dnbCI_con(i,2)*dnbCI_con(i,3)/dnbCI_con(i,4);
    end
end

dnbstruct.CI=CI;
dnbstruct.modules=modules;
dnbstruct.modIdx=modIdx;
dnbstruct.outIdx=outIdx;
dnbstruct.candidate=candidate;
dnbstruct.maxCI=maxCI;
dnbstruct.dnbIdx=dnbIdx;
dnbstruct.data=data;
dnbstruct.feature=feature;
dnbstruct.timeIdx=timeIdx;
dnbstruct.timeIdx_con=timeIdx_con;
dnbstruct.candi=candi;
dnbstruct.alpha=alpha;

%% plot 
%{
%}
if display==1
    ylabl={'CI','SD','PCCin','PCCout'};
    lw = 1.2;
    mk = '.';
    ms = 16;
    fs = 14;
    % DNB selection
    figure('Name','DNB selection')
    for i=1:4
        subplot(2,2,i);hold;
        plot(maxCI(:,i),'-k','linewidth',lw,'marker',mk,'MarkerSize',ms)
        plot(dnbIdx,maxCI(dnbIdx,i),'.r','MarkerSize',ms)
        ylabel(ylabl{i})
        set(gca,'xlim',[1,size(timeIdx,1)],'xtick',1:size(timeIdx,1))
        if i == 1
            ht = title('DNB selection');
            pos = get(ht, 'position');
            xlim = get(gca, 'xlim');
            set(ht, 'pos', pos + [xlim(2)-pos(1)+(xlim(2)-xlim(1))/6, 0, 0], 'HorizontalAlignment', 'center', 'FontSize', fs);
        end
        box on;
        grid on;
    end
    
    % DNB behavior
    % plot 4 components
    figure('Name','DNB behavior');
    for i=1:4
        subplot(2,2,i);hold;
        plot(dnbCI(:,i),'-k','linewidth',lw,'marker',mk,'MarkerSize',ms)
        %{
        if ~isempty(timeIdx_con)
            plot([dnbCI(:,i),dnbCI_con(:,i)])
            plot(dnbCI(:,i),'k','linewidth',1.5)
        end
        %}
        ylabel(ylabl{i})
        plot(dnbIdx,dnbCI(dnbIdx,i),'.r','MarkerSize',ms)
        set(gca,'xlim',[1,size(timeIdx,1)],'xtick',1:size(timeIdx,1))
        if i == 1
            ht = title('DNB behavior');
            pos = get(ht, 'position');
            xlim = get(gca, 'xlim');
            set(ht, 'pos', pos + [xlim(2)-pos(1)+(xlim(2)-xlim(1))/6, 0, 0], 'HorizontalAlignment', 'center', 'FontSize', fs);
        end
        box on;
        grid on;
    end
    
    %{
    % only plot CI
    figure;hold on
    plot(dnbCI(:,1),'-*k','linewidth',1,'MarkerSize',8)
    ylabel('CI')
    plot(dnbIdx,dnbCI(dnbIdx,1),'*r','linewidth',1,'MarkerSize',8)
    box on;
    set(gca,'xlim',[1,size(timeIdx,1)],'xtick',1:size(timeIdx,1))
    xlabel('Progression \rightarrow');
    %}
end

%% get dnbnet information
if nargout>=5
    dnbnet.node=zeros(size(data,1),size(timeIdx_c,1));
    eidx=tril(ones(size(data,1)),-1)==1;
    dnbnet.edge=zeros(sum(eidx(:)),size(timeIdx_c,1));
    
    for i=1:size(timeIdx_c,1)
        sd_tmp=std(data(:,timeIdx_c{i,2}),0,2);
        pcc_tmp=abs(corrcoef(data(:,timeIdx_c{i,2})'));
        pcc_tmp(isnan(pcc_tmp))=0;
        %{
    midx=modIdx{dnbIdx}{maxCI(dnbIdx,5)};
    oidx=outIdx{dnbIdx}{maxCI(dnbIdx,5)};
    temp=sd_tmp.*mean(pcc_tmp(:,midx),2)./mean(pcc_tmp(:,oidx),2);
    temp(isnan(temp))=0;
        %}
        temp=sd_tmp;
        dnbnet.node(:,i)=temp;
        dnbnet.edge(:,i)=pcc_tmp(eidx);
    end
    
    % normalize sd
    % dnbnet.node=quantilenorm(dnbnet.node);
    % dnbnet.node=zscore(dnbnet.node,0,2);
    
    dnbnet.node=[(1:numel(feature))',zeros(size(data,1),1),dnbnet.node];
    dnbnet.nodeName=feature;
    dnbnet.node(idx,2)=1;
    
    % normalize pcc
    dnbnet.edge=quantilenorm(dnbnet.edge);
    
    temp=sort(dnbnet.edge(:,1),'descend');
    edgeCrit=temp(3*size(data,1));
    idx=any(dnbnet.edge>=edgeCrit,2);
    
    % do not show outside connections
    temp=eidx;
    temp(dnbnet.node(:,2)==0,:)=false;
    temp=temp(eidx);
    idx=idx & temp;
    
    [x,y]=find(eidx);
    dnbnet.edge=[x,y,dnbnet.edge];
    dnbnet.edge=dnbnet.edge(idx,:);
    dnbnet.edgeCrit=edgeCrit;
    
    dnbnet.dnbCI=dnbCI;
    if exist('dnbCI_con','var')
        dnbnet.dnbCI_con=dnbCI_con;
    end
end
