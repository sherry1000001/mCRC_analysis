% compute the pcc matrix between variables from x and y
% rows of x and y are variables; columns are observations
% normalized by N-1, where N is the column numbers of x and y
% mode:0(default)---return size(x,1)*size(y,1) pcc matrix
%      1---return size(x,1)*1 pcc matrix. size(x)==size(y)
function pccMat=pcc_xy(x,y,mode)
if nargin<3
    mode=0;
end
N=size(x,2);
if mode==0;
    pccMat=(x-repmat(mean(x,2),1,N))*(y-repmat(mean(y,2),1,N))'/...
        (N-1)./(std(x,0,2)*std(y,0,2)');
end
if mode==1
    pccMat=sum((x-repmat(mean(x,2),1,N)).*(y-repmat(mean(y,2),1,N)),2)./...
        (std(x,0,2).*std(y,0,2))/(N-1);
end