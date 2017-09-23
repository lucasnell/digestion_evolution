% runs bootstrap if bootflag==1
bootflag=0;

% number of bootstrap data sets (2000 recommended for 95% confidence
% intervals
reps=200;

%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load BauwensDataRaw.txt /ascii
data=BauwensDataRaw;
% untransformed data
%1 body weight
%2 sprint speed
%3 hind limb

n=length(data);

%selection of traits
tr1=1;
tr2=3;
traits=[tr1 tr2];

load BauwensV.txt /ascii
SP=BauwensV;
SP=SP/det(SP)^(1/n);

% separate covariance matrices for different traits
CXbase=SP;
CYbase=SP;

%%%%%%%%%%%%%%%%%%%%%%%
% note assumes CX=CY; if CX!=CY, CXY=DX^-1*(DY')^-1 where
% DX*CX*DX'=eye(n)
%%%%%%%%%%%%%%%%%%%%%%%
CXYbase=SP;
	
load BauwensMERaw.txt /ascii
%untransformed SEs
ME=BauwensMERaw;

% log transform data with measurement error
% if X is lognormally distributed with mean M and variance V, X=exp(Y),
% then the mean m and variance v of Y are
%
%	v=log(V/M^2 +1)
%	m=log(M)-v/2

ME=log((ME./data).^2 + 1).^.5;

%%%%%%%%%%%%%%%%%
ME=.0001*ME
%%%%%%%%%%%%%%%%

MEX=ME(:,tr1);
MEY=ME(:,tr2);

MX=diag(MEX.^2);
MY=diag(MEY.^2);

W=log(data)-ME.^2/2;
X=W(:,tr1);
Y=W(:,tr2);

data=[Y X];