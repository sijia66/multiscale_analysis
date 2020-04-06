function Session = VisEP_Database()
%
%Session = VisEP_Database()
%
%database of sessions with full-field visual flashes to calculate evoked-potentials
%for Goose_Multiscale_M1_ECOG244
%
%output: Session - cell (#sessions x 1). Within each cell, format is:
%                  {DAY, {REC1, REC2, REC3, ....}, session#}


Session = cell(0,0);
ind = 1;

Session{ind} = {'180216',{'003', '006'},ind}; ind = ind+1;
Session{ind} = {'180217',{'004'},ind}; ind = ind+1;
Session{ind} = {'180218',{'005', '006', '030'},ind}; ind = ind+1;
Session{ind} = {'180219',{'041'},ind}; ind = ind+1; %probably more--missing page of notebook
Session{ind} = {'180220',{'004', '008'},ind}; ind = ind+1;
Session{ind} = {'180221',{'003', '008'},ind}; ind = ind+1;
Session{ind} = {'180222',{'003', '011'},ind}; ind = ind+1;
Session{ind} = {'180223',{'003','008'},ind}; ind = ind+1;

%Session{ind} = {'180321',{'005'},ind}; ind = ind+1;
%Session{ind} = {'180322',{'006', '011', '012'},ind}; ind = ind+1; %%issue w/ ecog connection
%Session{ind} = {'180323',{'005'},ind}; ind = ind+1; %why are there 212
%electrodes?
Session{ind} = {'180324',{'006', '009'},ind}; ind = ind+1;
Session{ind} = {'180325',{'003', '006'},ind}; ind = ind+1;
Session{ind} = {'180326',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180327',{'003', '005', '011'},ind}; ind = ind+1;
Session{ind} = {'180328',{'003', '005', '013'},ind}; ind = ind+1;
Session{ind} = {'180329',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180330',{'004', '006'},ind}; ind = ind+1;
Session{ind} = {'180331',{'003', '005', '011'},ind}; ind = ind+1;
Session{ind} = {'180401',{'003', '005', '011'},ind}; ind = ind+1;
Session{ind} = {'180403',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180404',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180405',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180406',{'003', '005', '011'},ind}; ind = ind+1;
%Session{ind} = {'180407',{'003', '005'},ind}; ind = ind+1; %task code is
%all 9
Session{ind} = {'180408',{'003', '005'},ind}; ind = ind+1;
%Session{ind} = {'180409',{'003', '005'},ind}; ind = ind+1; %task code is
%all 9
Session{ind} = {'180410',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180411',{'003', '005', '009'},ind}; ind = ind+1;
Session{ind} = {'180412',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180413',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180414',{'003', '005'},ind}; ind = ind+1;
Session{ind} = {'180415',{'003'},ind}; ind = ind+1;

%Session{ind} = {'180608',{},ind}; ind = ind+1;
%Session{ind} = {'180609',{},ind}; ind = ind+1;
%Session{ind} = {'180610',{},ind}; ind = ind+1;


%Session{ind} = {'180724',{},ind}; ind = ind+1;
Session{ind} = {'180725',{'003'},ind}; ind = ind+1;
Session{ind} = {'180726',{'003', '004'},ind}; ind = ind+1;
Session{ind} = {'180727',{'003', '004'},ind}; ind = ind+1;
Session{ind} = {'180728',{'003'},ind}; ind = ind+1;
Session{ind} = {'180729',{'003'},ind}; ind = ind+1;
Session{ind} = {'180730',{'004'},ind}; ind = ind+1;
Session{ind} = {'180731',{'003'},ind}; ind = ind+1;
Session{ind} = {'180801',{'003'},ind}; ind = ind+1;
Session{ind} = {'180802',{'021'},ind}; ind = ind+1;
Session{ind} = {'180803',{'002'},ind}; ind = ind+1;

Session{ind} = {'180807',{'004'},ind}; ind = ind+1;
Session{ind} = {'180808',{'004'},ind}; ind = ind+1;
Session{ind} = {'180809',{'003'},ind}; ind = ind+1;
Session{ind} = {'180810',{'003'},ind}; ind = ind+1;
Session{ind} = {'180813',{'003'},ind}; ind = ind+1;
Session{ind} = {'180814',{'005'},ind}; ind = ind+1;
Session{ind} = {'180815',{'004'},ind}; ind = ind+1;
Session{ind} = {'180816',{'003'},ind}; ind = ind+1;
Session{ind} = {'180817',{'002'},ind}; ind = ind+1;
