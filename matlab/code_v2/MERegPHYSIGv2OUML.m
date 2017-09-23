function [lnlike AIC esttrue boot]=MERegPHYSIGv2OUML(data,ME,Vall,bootflag,alpha,reps)global CX CY MX MY nX nY p tipsum ib jb mb nb counter col% number of tips and variables[n p]=size(data);Y=data(:,1);X=data(:,2:end);MEY=ME(:,1);MEX=ME(:,2:end);	MY=diag(MEY.^2);MX=diag(MEX(:).^2);CY=Vall(:,:,1);CX=[];for i=2:p	Crow=[];	for j=2:p						CXi=Vall(:,:,i);		CXj=Vall(:,:,j);				if sum(sum((CXi-CXj).^2)) < 10^-8			Crow=[Crow CXi];		else			invCXi=CXi\eye(n);			invCXj=CXj\eye(n);			iDXi=chol(CXi)';			iDXj=chol(CXj)';			CXij=iDXi*iDXj';			Crow=[Crow CXij];		end	end	CX=[CX;Crow];end% GLS estimates without measurement errorU=ones(n,1);aX=[];for i=2:p	invCX=Vall(:,:,i)\eye(n);	a=(U'*invCX*U)\(U'*invCX*X(:,i-1));	aX=[aX;a];ends2X=[];for i=2:p	ss2X=[];	for j=2:p		invCXY=CX(n*(i-2)+1:n*(i-1),n*(j-2)+1:n*(j-1))\eye(n);		ss2X=[ss2X (X(:,i-1)-aX(i-1))'*invCXY*(X(:,j-1)-aX(j-1))/(n-1)];	end	s2X=[s2X;ss2X];ends2X=s2X^.5;U2=[ones(n,1) X];invCY=CY\eye(n);b=(U2'*invCY*U2)\(U2'*invCY*Y);s2E=(Y-U2*b)'*invCY*(Y-U2*b)/(n-p);sE=s2E^.5;s2=[];for i=1:p-1									s2=[s2 s2X(i,i:end)];endif isempty(s2)	s2=0;end% REML estimates with measurement errortipsum=diag(CY)*ones(1,length(CY)) + ones(length(CY),1)*diag(CY)'-2*CY;nX=X;	nY=Y;	startd=.5;startp=[b(2:end)' sE aX' s2(2:end) startd];		%%%%%%%%%%%%% fixed[mb,nb]=size(startp);[ib,jb,startp]=find(startp);ii=0;while MERegPHYSIGv2OUMLfunct(startp)>=10^10 & ii <50	ii=ii+1;	s2Xtemp=s2X^2;		dsX=diag(diag(s2Xtemp).^-.5);	R=dsX*s2Xtemp*dsX;	R=R.*(.01*ii*eye(p-1)+(1-.01*ii)*ones(p-1,p-1));		dsX=diag(diag(s2Xtemp).^.5);	temps2X=dsX*R*dsX;		temps2X=temps2X^.5;	s2=[];	for i=1:p-1		s2=[s2 temps2X(i,i:end)];	end	startp=[b(2:end)' sE aX' s2(2:end) startd];		%%%%%%%%%%%%%% fixedendcounter=0;annealingopts=struct('StopTemp',10^-4,'Verbosity',0);fhandle=@MERegPHYSIGv2OUMLfunct;[startpAnneal LL]=anneal(fhandle,startp,annealingopts);[startpAnneal LL]=anneal(fhandle,startpAnneal,annealingopts);[startpAnneal LL]=anneal(fhandle,startpAnneal,annealingopts);options=optimset('tolX',10^-6,'tolFun',10^-10,'MaxFunEvals',10^4,'MaxIter',10^4);[startpFminsearch LL]=fminsearch('MERegPHYSIGv2OUMLfunct',startpAnneal,options);[est LL]=fminsearch('MERegPHYSIGv2OUMLfunct',startpFminsearch,options);% MLlnlike=-LL - p*n*log(2*pi)/2;% REML% lnlike=-LL - (n-length(b))*log(2*pi)/2 + log(det(n^p))/2est=full(sparse(ib,jb,est,mb,nb));b=est(1:p-1)';	sE=abs(est(p));s2=est(p+1:end-1);									%%%%%%%%%%%%% fixedd=abs(est(end));CYd=real(d.^tipsum.*(1-d.^(2*CY))/(1-d^2));startboot=est;if startboot(end)<.1	startboot(end)=.1;end% extract covariance matrixs2X=zeros(p-1,p-1);ii=1;for i=1:p-1	for j=i:p-1		s2X(i,j)=s2(ii);		ii=ii+1;	endends2X=s2X+s2X'-diag(diag(s2X));s2X=s2X^2;sX=diag(s2X).^.5;ds2X=diag(diag(s2X).^(-.5));R=ds2X*s2X*ds2X;% extract estimates of the meansvX=R.*(sX*sX');VX=CX.*kron(vX,ones(n,n));VX=VX+MX;vXY=vX.*(b*ones(1,p-1));vXY=CX.*kron(vXY,ones(n,n));vYY=vX.*(b*b');vYY=CX.*kron(vYY,ones(n,n));VXY=zeros(n*(p-1),n);VYY=zeros(n,n);for i=1:p-1	VXY=VXY + vXY(:,n*(i-1)+1:n*i);	for j=1:p-1		VYY=VYY + vYY(n*(i-1)+1:n*i,n*(j-1)+1:n*j);	endendVYY=VYY+sE^2*CYd+MY;VV=[VX VXY;VXY' VYY];invVV=VV\eye(p*n);Z=kron(eye(p),ones(n,1));w=kron(eye(p),[nX nY]);W=w(:,1:p+1:end);a=(Z'*invVV*Z)\(Z'*invVV*W);a=diag(a);if p==1	aY=a(p);	b0=aY;	est=[b0 sE^2 d];else	aX=a(1:p-1);	aY=a(p);	b0=aY-b'*aX;	est=[b0 b' sE^2 aX' s2 d];endesttrue=est;estR=R;aXtrue=aX;btrue=[b0;b];dtrue=d;% GLS approximate standard errorsU2=[ones(n,1) X];invVYY=VYY\eye(n);bGLS=(U2'*invVYY*U2)\(U2'*invVYY*Y);s2EY=(Y-U2*bGLS)'*invVYY*(Y-U2*bGLS)/(n-p);s2b=s2EY*(U2'*invVYY*U2)^-1;bSE=diag(s2b).^.5;bGLS_approx_conf_int_for_b=[[b0;b] bGLS-2*bSE bGLS bGLS+2*bSE]% if dtrue<.01% 	bootflag=0;% endif ~bootflag	bootconf=[];endif bootflag		if p>1		VX=kron(R.*(sX*sX'),ones(n,n)).*CX;    %%%%% fixed 30Jun17		[TX LX TT]=svd(VX);		DX=TX*LX.^(.5);	end		VY=s2E*CYd;	[TY LY TT]=svd(VY);	DY=TY*LY.^(.5);	bootlist=[];	for t=1:reps				if p==1			EY=DY*randn(n,1);			nY=ones(n,1)*btrue+EY+MEY.*randn(n,1);		else			EX=DX*randn(n*(p-1),1);			EY=DY*randn(n,1);			a=ones(n,1)*aXtrue';			rX=a(:)+EX;			nX=rX+MEX(:).*randn(n*(p-1),1);				rX=reshape(rX,n,p-1);			nX=reshape(nX,n,p-1);			nY=[ones(n,1) rX]*btrue+EY+MEY.*randn(n,1);		end		[ib,jb,startpboot]=find(startboot);		annealingopts=struct('StopTemp',10^-3,'Verbosity',0);		fhandle=@MERegPHYSIGv2OUMLfunct;		[startpAnneal LL]=anneal(fhandle,startpboot,annealingopts);		options=optimset('tolX',10^-6,'tolFun',10^-10,'MaxFunEvals',10^4,'MaxIter',10^4);		[estboot LL]=fminsearch('MERegPHYSIGv2OUMLfunct',startpAnneal,options);		[estboot LL]=fminsearch('MERegPHYSIGv2OUMLfunct',estboot,options);		%[estboot LL]=fminsearch('MERegPHYSIGv2OUMLfunct',startpboot,options);				estboot=full(sparse(ib,jb,estboot,mb,nb));		b=estboot(1:p-1)';			sE=abs(estboot(p));		s2=estboot(p+1:end-1);						%%%%%%%%%%%%% fixed		d=abs(estboot(end));		bootCYd=real(d.^tipsum.*(1-d.^(2*CY))/(1-d^2));				% extract covariance matrix		s2X=zeros(p-1,p-1);		ii=1;		for i=1:p-1			for j=i:p-1				s2X(i,j)=s2(ii);				ii=ii+1;			end		end		s2X=s2X+s2X'-diag(diag(s2X));		s2X=s2X^2;		sX=diag(s2X).^.5;		ds2X=diag(diag(s2X).^(-.5));		R=ds2X*s2X*ds2X;		% extract estimates of the means		vX=R.*(sX*sX');		VX=CX.*kron(vX,ones(n,n));		VX=VX+MX;		vXY=vX.*(b*ones(1,p-1));		vXY=CX.*kron(vXY,ones(n,n));		vYY=vX.*(b*b');		vYY=CX.*kron(vYY,ones(n,n));		VXY=zeros(n*(p-1),n);		VYY=zeros(n,n);		for i=1:p-1			VXY=VXY + vXY(:,n*(i-1)+1:n*i);			for j=1:p-1				VYY=VYY + vYY(n*(i-1)+1:n*i,n*(j-1)+1:n*j);			end		end		VYY=VYY+sE^2*bootCYd+MY;		VV=[VX VXY;VXY' VYY];		invVV=VV\eye(p*n);		Z=kron(eye(p),ones(n,1));		w=kron(eye(p),[nX nY]);		W=w(:,1:p+1:end);		a=(Z'*invVV*Z)\(Z'*invVV*W);		a=diag(a);		if p==1			aY=a(p);			b0=aY;			bootlist=[bootlist;b0 sE^2 d];		else			aX=a(1:p-1);			aY=a(p);			b0=aY-b'*aX;			bootlist=[bootlist;b0 b' sE^2 aX' s2 d];		end				%[t bootlist(t,:)]	end		figure(100)	hist(bootlist(:,2))	hold on	plot([btrue(2) btrue(2)],[0 300],'r')	hold off		drawnow		% obtain alpha% confidence intervals	[bn, bc]=size(bootlist);	bootconf=[];	for jj=1:bc		w=bootlist(:,jj);		w=sort(w);		m=length(w);		ub=w(ceil((1-.5*alpha)*m));		lb=w(floor(.5*alpha*m));		meanw=mean(w);		bootconf=[bootconf; lb meanw ub];	endend		% print final valuesif bootflag	fprintf('Parameters with (lb,mean,ub) from simulation\n\n')	fprintf(['b0 (intercept) = ',num2str(est(1)),' (',...		num2str(bootconf(1,1)),', ',num2str(bootconf(1,2)),', ',num2str(bootconf(1,3)),')\n']);	for i=1:p-1		fprintf(['            b',num2str(i)',' = ',num2str(est(i+1)),' (',...			num2str(bootconf(i+1,1)),', ',num2str(bootconf(i+1,2)),', ',num2str(bootconf(i+1,3)),')\n']);	end		fprintf(['\n        sigma2 = ',num2str(est(p+1)),' (',...		num2str(bootconf(p+1,1)),', ',num2str(bootconf(p+1,2)),', ',num2str(bootconf(p+1,3)),')\n']);	fprintf(['\n             d = ',num2str(est(end)),' (',...		num2str(bootconf(end,1)),', ',num2str(bootconf(end,2)),', ',num2str(bootconf(end,3)),')\n\n']);else	fprintf('Parameters\n\n')	fprintf(['b0 (intercept) = ',num2str(est(1)),'\n']);	for i=1:p-1		fprintf(['            b',num2str(i)',' = ',num2str(est(i+1)),'\n']);	end		fprintf(['\n        sigma2 = ',num2str(est(p+1)),'\n']);	fprintf(['\n             d = ',num2str(est(end)),'\n\n']);endfprintf('Independent variable means, variances, and correlations\n')for i=1:p-1	fprintf(['aX',num2str(i)',' = ',num2str(est(p+1+i)),'\n']);endfor i=1:p-1	fprintf(['s2X',num2str(i)',' = ',num2str(est(2*p+i)),'\n']);endfor i=1:p-1	for j=i+1:p-1		fprintf(['r(',num2str(i),',',num2str(j)',') = ',num2str(estR(i,j)),'\n']);	endendfprintf(['\nLnLikelihood = ',num2str(lnlike),'\n'])fprintf(['-2LL = ',num2str(-2*lnlike),'\n'])par=length(est)-1-length(b)-length(aX');AIC=-2*lnlike + 2*par;fprintf(['AIC(par=',num2str(par),') = ',num2str(AIC),'\n\n'])if bootflag	boot=bootconf;else	boot=zeros(length(startp)+1,3);end