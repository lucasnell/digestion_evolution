% Version 13 Nov, 2007% updated 1Jun09function MERegEGLS(data,ME,Vall,bootflag,alpha,reps)% number of tips and variables[n p]=size(data);Y=data(:,1);X=data(:,2:end);MEY=ME(:,1);MEX=ME(:,2:end);MY=diag(MEY.^2);MX=diag(MEX.^2);CY=Vall(:,:,1);for i=2:p	for j=2:p						CXi=Vall(:,:,i);		CXj=Vall(:,:,j);				if sum(sum((CXi-CXj).^2)) > 10^-8			fprintf('Warning: covariance matrices for X variables are not identical;\n')			fprintf('    EGLS invalid so use ML or REML\n')		end	endend%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CX is assumed to be the same for all independent variablesCX=Vall(:,:,2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLS estimates without measurement errorU=ones(n,1);invCX=CX\eye(n);aX=[];s2X=[];for i=2:p	a=(U'*invCX*U)\(U'*invCX*X(:,i-1));	aX=[aX;a];	s2X=[s2X; (X(:,i-1)-a)'*invCX*(X(:,i-1)-a)/(n-1)];endU2=[ones(n,1) X];invCY=CY\eye(n);b=(U2'*invCY*U2)\(U2'*invCY*Y);s2E=(Y-U2*b)'*invCY*(Y-U2*b)/(n-p);% EGLS estimates with measurement error% function that returns the means and covariance matrix from the independent% variables using pairwise EGLSif p-1>1	[aX, COVX]=MECorrForRegEGLS(X,MEX,Vall(:,:,2:end));	R=COVX./((diag(COVX).^.5)*(diag(COVX).^.5)');else	%find mean using EGLS (from MEUnivarEGLS.m)	V=CX;	olds2X=0;	dSS=1;	t=0;	while dSS>10^-16 & t<50,		t=t+1;		invV=V\eye(n);		aX=(U'*invV*U)\(U'*invV*X);		s2X=(X-aX)'*invV*(X-aX)./(n-1);		V=CX + (n/(n-1))*MX/s2X;		dSS=(s2X-olds2X)^2;		olds2X=s2X;	end	aX=aX;	COVX=s2X;	R=1;end% EGLS for regression coefficientsVY=CY;t=0;oldb=ones(p-1,1);dSS=1;while dSS>10^-8 & t<10^3,	t=t+1;	invVY=VY\eye(n);	aY=(U'*invVY*U)\(U'*invVY*Y);	s2E=(Y-aY)'*invVY*(Y-aY)/(n-1);	c=[];	for i=1:p-1		c=[c;(X(:,i)-ones(n,1)*aX(i)')'*invCX*(Y-aY)/(n-1)];	end	b=(COVX^-1)*c;	CCX=(b'*COVX*b)*CX;	VY=CY + CCX/s2E + (n/(n-1))*MY/s2E;			dSS=(b-oldb)'*(b-oldb)/length(b);	oldb=b;endb0=aY-b'*aX;sE=s2E^.5;sX=s2X.^.5;est=[b0 b' s2E aX' s2X'];estR=R;aXtrue=aX;btrue=[b0;b];% approximate standard errorsif ~bootflag	b=(U2'*invVY*U2)\(U2'*invVY*Y);	s2EY=(Y-U2*b)'*invVY*(Y-U2*b)/(n-p);	s2b=s2EY*(U2'*invVY*U2)^-1;	bSE=diag(s2b).^.5;endif bootflag	% note that VX is created assuming CX is the same for all independent	% variables	VX=kron(R.*(sX'*sX),CX);	[TX LX TT]=svd(VX);	DX=TX*LX.^(.5);		VY=s2E*CY;	[TY LY TT]=svd(VY);	DY=TY*LY.^(.5);	bootlist=[];	for t=1:reps		EX=DX*randn(n*(p-1),1);		EY=DY*randn(n,1);					a=ones(n,1)*aXtrue';		rX=a(:)+EX;		nX=rX+MEX(:).*randn(n*(p-1),1);			rX=reshape(rX,n,p-1);		nX=reshape(nX,n,p-1);		nY=[ones(n,1) rX]*btrue+EY+MEY.*randn(n,1);		% GLS estimates without measurement error		invCX=CX\eye(n);		s2X=[];		for i=2:p			a=(U'*invCX*U)\(U'*invCX*nX(:,i-1));			s2X=[s2X; (nX(:,i-1)-a)'*invCX*(nX(:,i-1)-a)/(n-1)];		end		U2=[ones(n,1) nX];		invCY=CY\eye(n);		b=(U2'*invCY*U2)\(U2'*invCY*nY);		s2E=(nY-U2*b)'*invCY*(nY-U2*b)/(n-p);		% EGLS estimates with measurement error		% function that returns the means and covariance matrix from the independent		% variables using pairwise EGLS		if p-1>1			[aX, COVX]=MECorrForRegEGLS(nX,MEX,Vall(:,:,2:end));			R=COVX./((diag(COVX).^.5)*(diag(COVX).^.5)');		else			%find mean using EGLS (from MEUnivarEGLS.m)			V=CX;			olds2X=0;			dSS=1;			t=0;			while dSS>10^-16 & t<10^3,				t=t+1;				invV=V\eye(n);				aX=(U'*invV*U)\(U'*invV*nX);				s2X=(nX-U*aX)'*invV*(nX-U*aX)./(n-1);				V=CX + (n/(n-1))*MX/s2X;				dSS=(s2X-olds2X)^2;				olds2X=s2X;			end			aX=aX;			COVX=s2X;		end		% EGLS for regression coefficients		VY=CY;		t=0;		oldb=ones(p-1,1);		dSS=1;		while dSS>10^-16 & t<50,			t=t+1;			invVY=VY\eye(n);			aY=(U'*invVY*U)\(U'*invVY*nY);			s2E=(nY-aY)'*invVY*(nY-aY)/(n-1);			c=[];			for i=1:p-1				c=[c;(nX(:,i)-ones(n,1)*aX(i)')'*invCX*(nY-aY)/(n-1)];			end			b=COVX\c;			CCX=(b'*COVX*b)*CX;			VY=CY + CCX/s2E + (n/(n-1))*MY/s2E;			dSS=(b-oldb)'*(b-oldb)/length(b);			oldb=b;		end		b0=aY-b'*aX;		bootlist=[bootlist;b0 b' s2E aX' s2X'];	end		% obtain alpha% confidence intervals	[bn, bc]=size(bootlist);	bootconf=[];	for jj=1:bc		d=bootlist(:,jj);		d=sort(d);		m=length(d);		ub=d(ceil((1-.5*alpha)*m));		lb=d(floor(.5*alpha*m));		meand=mean(d);		bootconf=[bootconf; lb meand ub];	endend		% print final valuesif bootflag	fprintf('Coefficients with (lb,mean,ub) from simulation\n')	fprintf(['b0 (intercept) = ',num2str(est(1)),' (',...		num2str(bootconf(1,1)),', ',num2str(bootconf(1,2)),', ',num2str(bootconf(1,3)),')\n']);	for i=1:p-1		fprintf(['b',num2str(i)',' = ',num2str(est(i+1)),' (',...			num2str(bootconf(i+1,1)),', ',num2str(bootconf(i+1,2)),', ',num2str(bootconf(i+1,3)),')\n']);	end		fprintf(['\nsigma2 = ',num2str(est(p+1)),' (',...		num2str(bootconf(p+1,1)),', ',num2str(bootconf(p+1,2)),', ',num2str(bootconf(p+1,3)),')\n\n']);else	fprintf('Coefficients with approximate SE\n')	fprintf(['b0 (intercept) = ',num2str(est(1)),' +- ',num2str(bSE(1)),'\n']);	for i=1:p-1		fprintf(['b',num2str(i)',' = ',num2str(est(i+1)),' +- ',num2str(bSE(i+1)),'\n']);	end		fprintf(['\nsigma2 = ',num2str(est(p+1)),'\n\n']);endfprintf('Independent variable means, variances, and correlations\n')for i=1:p-1	fprintf(['aX',num2str(i)',' = ',num2str(est(p+1+i)),'\n']);endfor i=1:p-1	fprintf(['s2X',num2str(i)',' = ',num2str(est(2*p+i)),'\n']);endfor i=1:p-1	for j=i+1:p-1		fprintf(['r(',num2str(i),',',num2str(j)',') = ',num2str(estR(i,j)),'\n']);	endend