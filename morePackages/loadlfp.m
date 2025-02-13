function lfp = loadlfp(file, Events, trials, field, bn, CH, FS)
%  LOADLFP loads lfp data 
%
%  LFP = LOADLFP(FILE, EVENTS, TRIALS, FIELD, BN, CH, FS)
%
%  Inputs:  FILE    = String.  Lfp data file prefix.
%           EVENTS  = Structure.  Trial events data structure.
%           TRIALS  = Vector.  Trials to load data for.
%           FIELD   = Scalar.  Event to align data to.
%           BN      = Vector.  Time to start and stop loading data.
%           CH      = Scalar.  Number of recording channels
%           FS
%           Note:  Times should be ms.
%
%   Outputs:    LFP  = [TRIAL,CH,TIME] Array. LFP data.
%

%  Written by:  Bijan Pesaran
%


if ~ischar(field); error('FIELD needs to be a string'); end
if nargin < 6; CH = 2; end
if nargin < 7; FS = 1e3; end

fid = fopen([file '.lfp.dat']);
ss = 4;
format = 'float=>single';
bn = double(bn);
ntr = length(trials);
%disp(num2str(ntr));
N1 = round(bn(1)*FS./1e3);
N2 = round(bn(2)*FS./1e3);
dN=N2-N1;
lfp = zeros(ntr,CH,dN,'single');

%Events.StartOn(trials)
if strcmp(field,'PulseStarts')
    pulse = 0;
    for i = 1:ntr
        pulseTimes = Events.PulseStarts{trials(i)};
        for iPulse = 1:length(pulseTimes)
            pulse = pulse + 1;
            at = round(pulseTimes(iPulse));
            start = at+bn(1);
            pos = round(start.*FS./1e3).*ss.*CH;
            status = fseek(fid,pos,'bof');
            if status < 0; ferror(fid); end
            h = fread(fid,[CH,dN],format);
            lfp(pulse,:,1:size(h,2)) = h;
        end
    end
else
    for i = 1:ntr
        at = getfield(Events,field,{trials(i),1});
        start = at + bn(1);
        pos = round(start.*FS./1e3).*ss.*CH;
        status = fseek(fid,pos,'bof');

        if status < 0
            ferror(fid)
            disp(['error in reading: ' file ' position: ' num2str(pos)])
        %else
            %disp(['ok in reading: ' file ' position: ' num2str(pos)])
        end
            
        h = fread(fid, [CH,dN], format);
        
        lfp(i,:,1:size(h,2)) = h;
    end
end
fclose(fid);
% 
% if CH==1 || ntr == 1
%     mlfp = squeeze(mlfp);
% end
