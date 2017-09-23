% "MERegPHYSIGv2.m" Matlab program.% Version of 21 September 2013.% This is part of the PHYSIG package.% Calls programs	% RegLS.m	% RegGLS.m	% RegOU.m		% RegOUfunct		% RegOUMLfunct	% RegACDC.m		% RegACDCfunct		% RegACDCMLfunct	% RegGrafen.m		% RegGrafenfunct		% RegGrafenMLfunct	% RegPagel.m		% RegPagelfunct		% RegPagelMLfunctfprintf('MERegPHYSIGv2.m\n');fprintf('(c) Anthony R. Ives and Ted Garland Jr. - 9 Aug, 2006\n');% fprintf('See Johnson, Ives, Ahern, and Salminen ()\n');% fprintf('Version updated 21Sep2013\n\n');%format longlogon = input('Do you wish to log the session (Turn diary on)? (Y/N=<return>) ', 's');if lower(logon) == 'y'	fname = input('What would you like to call the diary file? (default: MERegPHYSIGv2.doc) ','s');	if isempty(fname)		fname = 'MERegPHYSIGv2.doc';	end;	diary(fname);end;fprintf('The date is: %s\n', date);c=clock;c=c(4:5);fprintf('Time (24 hour clock): %i:%.2i\n', c(1), c(2));%%%%%%%%%%%%%%%%%%%% input data file %%%%%%%%%%%%%%%%%%%%%%%fprintf('Data files should contain values corresponding to species (rows)\n');fprintf('    and variables (columns)\n');input('Hit Return to choose the data file.');[filename, path] = uigetfile('*.*', 'Open a data file!');fileinp=strcat(path, filename);fprintf('You chose data file: %s\n', fileinp);columns = input('Input the number of columns in the file, including tip names if present: ');fid1 = fopen(fileinp, 'r');header = input('Does the data file have a header row? (Y/N) ','s');if lower(header) == 'y'	fgetl(fid1);end;tips = input('Does the file include tip names in the first column? (Y/N) ', 's');% Read in the tip data.  Dynamically construct the format string.formstring = ' ';if lower(tips) == 'y'	columns = columns - 1;	for i = 1:columns		formstring=strcat('%e', formstring);	end;	[f] = fscanf(fid1, strcat('%*s', formstring), [columns,inf]); %skip over tip nameselse	for i = 1:columns		formstring = strcat('%e', formstring);	end	[f] = fscanf(fid1, formstring, [columns, inf]);endfclose(fid1);dataFile=f';[numtips numvars]=size(dataFile);if lower(tips) == 'y'	colcount=numvars+1;else	colcount=numvars;end;if lower(header) == 'y'	rowcount=numtips+1;else	rowcount=numtips;end;fprintf(['Your tip data file contains ',num2str(rowcount),' rows and ',num2str(colcount),' columns\n']);%%%%%%%%%%%%%%%%%%%% input measurement error file %%%%%%%%%%%%%%%%%%%%%%%fprintf('\nThe measurement error file contains SEs in columns corresponding to the data\n');fprintf('    file; if there is no measurement error, the column should contain zeros.\n');input('Hit return to choose the measurement error file containing standard errors.');[filename, path] = uigetfile('*.*', 'Open a measurement error file!');fileinp=strcat(path, filename);fprintf('Chose your measurement error file: %s\n', fileinp);fid1 = fopen(fileinp, 'r');MEheader = input('Does the measurement error file have a header row? (Y/N) ','s');if lower(MEheader) == 'y'	fgetl(fid1);end;MEtips = input('Does the file include tip names in the first column? (Y/N) ', 's');% Read in the tip data.  Dynamically construct the format string.formstring = ' ';if lower(MEtips) == 'y'	for i = 1:columns		formstring=strcat('%e', formstring);	end	%skip over tip names	[f] = fscanf(fid1, strcat('%*s', formstring), [columns,inf]); else	for i = 1:columns		formstring = strcat('%e', formstring);	end	[f] = fscanf(fid1, formstring, [columns, inf]);endfclose(fid1);MEfile=f';[MEnumtips MEnumvars]=size(MEfile);if lower(MEtips) == 'y'	colcount=MEnumvars+1;else	colcount=MEnumvars;end;if lower(MEheader) == 'y'	rowcount=MEnumtips+1;else	rowcount=MEnumtips;end;fprintf(['Your measurement error file contains ',num2str(rowcount),' rows and ', ...	num2str(colcount),' columns\n']);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		LOOP OVER SELECTION OF VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%continueflag='y';while lower(continueflag)=='y'		columnnum = input('\nHow many independent variables do you want to analyze? ');	columnnum=columnnum+1;		Vall=zeros(numtips,numtips,columnnum);	ME=[];	colnum=1;	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	%%%%%%%%%%%%%%%%%%%% select dependent variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		fprintf(['\nNote: If there are tip names in the first column, the numbering is determined\n']);	fprintf(['      after this column is removed (which is different from REGRESSIONv2.m)\n\n']);	columninp = input(['Which column contains the dependent variable? ']);% 	if lower(header)=='y'% 		columninp=columninp-1;% 	end	data=dataFile(:,columninp);	MEcol=MEfile(:,columninp);	tf=input('Do you want to log-transform this variable? (Y/N) ','s');	if lower(tf)=='y'		% log transform data with measurement error		% if X is lognormally distributed with mean M and variance V, X=exp(Y),		% then the mean m and variance v of Y are		%		%	v=log(V/M^2 +1)		%	m=log(M)-v/2		ME=[ME log((MEcol./data(:,colnum)).^2 + 1).^.5];		data(:,colnum)=log(data(:,colnum))-MEcol.^2/2;	else		ME=[ME MEcol];	end	%%%%%%%%%%%%%%%%%%%% input covariance matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%	fprintf('You can select either a star phylogeny or input your own:\n');	fprintf('             Star (S)\n');	fprintf('            Input (I)\n');	phylo = input('Select a phylogeny: ','s');	switch lower(phylo)		case {'s'}			Vall(:,:,colnum)=eye(numtips);		case {'i'}									input('Hit return to choose the covariance matrix file.');			[file, path] = uigetfile('*.*', 'Open a matrix file!');			matinp = strcat(path, file);			fprintf('You choose matrix file: %s\n', matinp);			fid2 = fopen(matinp, 'r');			Vheader = input('Does the matrix have a header row? (Y/N) ','s');			if lower(Vheader) == 'y'				fgetl(fid2);			end			mat = fscanf(fid2, '%e', [numtips, inf]);			%standardize size of elements of mat			mat=mat/mat(1,1);			mat=mat/det(mat)^(1/length(mat));			Vall(:,:,colnum)=mat;	end	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	%%%%%%%%%%%%%%%%%%%% select independent variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	for colnum=2:columnnum		fprintf('\n');		columninp = input(['Which column contains independent variable ',num2str(colnum-1),'? ']);% 		if lower(header)=='y'% 			columninp=columninp-1;% 		end		data=[data dataFile(:,columninp)];		MEcol=MEfile(:,columninp);		tf=input('Do you want to log-transform this variable? (Y/N) ','s');		if lower(tf)=='y'			% log transform data with measurement error			% if X is lognormally distributed with mean M and variance V, X=exp(Y),			% then the mean m and variance v of Y are			%			%	v=log(V/M^2 +1)			%	m=log(M)-v/2			ME=[ME log((MEcol./data(:,colnum)).^2 + 1).^.5];			data(:,colnum)=log(data(:,colnum))-MEcol.^2/2;		else			ME=[ME MEcol];		end				%%%%%%%%%%%%%%%%%%%% input covariance matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%		fprintf('You can select either a star phylogeny, choose the one entered for\n');		fprintf('        the dependent variable, or input a new one:\n');		fprintf('	                   Star (S)\n');		fprintf('	                  Input (I)\n');		fprintf('    Same as dependent variable (V)\n');		phylo = input('Select a phylogeny: ','s');		switch lower(phylo)			case {'s'}				Vall(:,:,colnum)=eye(numtips);			case {'v'}				Vall(:,:,colnum)=Vall(:,:,1);			case {'i'}									input('Hit return to choose the covariance matrix file.');				[file, path] = uigetfile('*.*', 'Open a matrix file!');				matinp = strcat(path, file);				fprintf('You chose matrix file: %s\n', matinp);				fid2 = fopen(matinp, 'r');				Vheader = input('Does the matrix have a header row? (Y/N) ','s');				if lower(Vheader) == 'y'					fgetl(fid2);				end;				mat = fscanf(fid2, '%e', [numtips, inf]);				%standardize size of elements of mat				mat=mat/mat(1,1);				mat=mat/det(mat)^(1/length(mat));				Vall(:,:,colnum)=mat;		end	end	%%%%%%%%%%%%%%%%%%%% Select a method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	fprintf('Methods:\n');	fprintf('              GLS (no ME) (G)\n');	fprintf('                     EGLS (E) [In multiple regression, this requires the covariance\n');	fprintf('                         matrices of all X variables to be the same]\n');	fprintf('		                ML (M)\n');	fprintf('                     REML (R)\n');	fprintf('              OU with ML (OM)\n');	fprintf('            OU with REML (OR)\n');	fprintf('   Pagels lambda with ML (LM)\n');	fprintf(' Pagels lambda with REML (LR)\n');	method = input('Select a method: ','s');	switch lower(method)		case {'g'}			methodstr='GLS';		case {'e'}			methodstr='EGLS';		case {'m'}			methodstr='ML';		case {'om'}			methodstr='OUML';		case {'or'}			methodstr='OUREML';		case {'lm'}			methodstr='PLML';		case {'lr'}			methodstr='PLREML';	end	bootflag=0;	reps=[];	alpha=[];	boot=input('Do you want obtain confidence intervals by simulation? (Y/N) ','s');	switch lower(boot)		case {'y'}			bootflag=1;			reps = input('Input the number of simulations you want to run (default = 2000): ');			if isempty(reps)				reps = 2000;			end			alpha=input('Select an alpha value for confidence intervals (default = .05): ');			if isempty(alpha)				alpha = .05;			end    end        fprintf('\nStarted at %s\n\n', datestr(now, 'HH:MM:SS'));	fprintf('\nOUTPUT FROM REGRESSIONv3.m\n');	fprintf(['METHOD = ',methodstr,'\n']);		fprintf('\nData file and ME file saved as ''workingDataFile.txt'' and ''workingMEfile.txt''\n\n');	save 'workingDataFile.txt' data -ascii	save 'workingMEfile.txt' ME -ascii	pgm=strcat('MERegPHYSIGv2',methodstr,'(data,ME,Vall,bootflag,alpha,reps)');	eval(pgm);        fprintf('\n\nEnded at %s\n', datestr(now, 'HH:MM:SS'));		continueflag=input('\nDo you want to perform additional analyses? (Y/N) ','s');endfprintf('\nEnd REGRESSIONv3.m\n');diary off