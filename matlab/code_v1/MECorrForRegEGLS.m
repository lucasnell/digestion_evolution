function [aX, COVX]=MECorrForRegEGLS(X,MEX,Vall)% MECorrForRegEGLS.m% Based on MECorrEGLS.m but used to calcualte the pairwise correlation between% independent variables for MERegEGLS.m% Anthony R. Ives% 9 Aug 06% number of tips and variables[n p]=size(X);U=ones(n,1);aX=[];s2X=[];COVX=zeros(p,p);for i=1:p-1	for j=i+1:p		X1=X(:,i);		X2=X(:,j);		CX1=Vall(:,:,i);		CX2=Vall(:,:,j);				if sum(sum((CX1-CX2).^2))<10^-8			sameCXflag=1;		else			sameCXflag=0;		end				MX1=diag(MEX(:,i).^2);		MX2=diag(MEX(:,j).^2);		% EGLS estimates with measurement error		VX1=CX1;		VX2=CX2;		s2X1=0;		s2X2=0;		dSS=1;		ii=1;		while dSS>10^-16 & ii<50,			ii=ii+1;			olds2X1=s2X1;			olds2X2=s2X2;						invVX1=VX1\eye(n);			invVX2=VX2\eye(n);			a1=(U'*invVX1*U)\(U'*invVX1*X1);			a2=(U'*invVX2*U)\(U'*invVX2*X2);			s2X1=(X1-U*a1)'*invVX1*(X1-U*a1)/(n-1);			s2X2=(X2-U*a2)'*invVX2*(X2-U*a2)/(n-1);			VX1=CX1+(n/(n-1))*MX1/s2X1;			VX2=CX2+(n/(n-1))*MX2/s2X2;						dSS=(s2X1-olds2X1)^2 + (s2X2-olds2X2)^2;		end				sX1=s2X1^.5;		sX2=s2X2^.5;		% calculate correlation			invCX1=CX1\eye(n);		invCX2=CX2\eye(n);		if sameCXflag			CX12=CX1;		else			iDX1=chol(CX1)';			iDX2=chol(CX2)';			DX1=iDX1\eye(n);			DX2=iDX2\eye(n);			CX12=iDX1*iDX2';		end		invCX12=CX12\eye(n);				c=(X2-a2)'*invCX12*(X1-a1)/(n-1);		r=c/(sX1*sX2);		% out of bounds estimate of r, so correct with Method of Moments		if abs(r)>1,			r=sign(r);			sX1n=sX1*c/(sX1*sX2);			sX2n=sX2*c/(sX1*sX2);			sX1=sX1n;			sX2=sX2n;		end				COVX(i,j)=r*sX1*sX2;	end	aX=[aX;a1];	s2X=[s2X;s2X1];endaX=[aX;a2];s2X=[s2X;s2X2];COVX=COVX+COVX'+diag(s2X);