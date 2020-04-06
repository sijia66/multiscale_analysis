function experiment = expmap_mapChannels(experiment)

%experiment = expmap_mapChannels(experiment)
%uses specified bank definition and daq mapping to compute channel mapping 
%i.e. updates channelid and multiunitid for each electrode given the set acquisition configuration. 

muidoffset = experiment.hardware.acquisition(1).muidOffset;

nBank = size(experiment.hardware.acquisition.bank,2);


%get daq mapping from pin -> acquisitionid
banklist = [experiment.hardware.acquisition(1).channel.bankid];
pinlist  = [experiment.hardware.acquisition(1).channel.pinid];
acqChlist = 1:length(banklist);
nCh = length(banklist);


%also get info on where banks are plugged in
pcbs  = [experiment.hardware.acquisition(1).bank(:).pcbid];
drives = [experiment.hardware.acquisition(1).bank(:).driveid];

%also get headstage gain info for each electrode (if headstage gains in acquisition are defined)
if isfield(experiment.hardware.acquisition(1).bank(1).headstage, 'gain') 
    tmp = [experiment.hardware.acquisition(1).bank(:).headstage];
    gains = [tmp(:).gain];
    
    %if 'amp' also included, factor that in as well. total gain = headstage_gain*amp_gain
    if isfield(experiment.hardware.acquisition(1).bank(1), 'amp')
        tmp = [experiment.hardware.acquisition(1).bank(:).amp];
        gains = gains.* [tmp(:).gain];
    end
    reassignGains = 1; %flag to use these gains. 
else
    gains = nan(size(pcbs));
    reassignGains = 0;
end


%convert to list of where each channel is plugged in
pcblist   = nan(1,nCh);
drivelist = nan(1,nCh);
gainlist  = nan(1,nCh);
for b=1:nBank
    ind = banklist == b;
    
    pcblist(ind)   = pcbs(b);
    drivelist(ind) = drives(b);
    gainlist(ind)  = gains(b); 
end



%loop through drive electrodes and assign ids accordingly

nDrive = size(experiment.hardware.microdrive,2);


for iD=1:nDrive
    
    nE = size(experiment.hardware.microdrive(iD).electrodes,2);
    
    for iE=1:nE
        
        %pcb and pin assignments for the electrode
        e_pcb = experiment.hardware.microdrive(iD).electrodes(iE).pcbid;
        e_pin = experiment.hardware.microdrive(iD).electrodes(iE).pinid;

        %find where e_pin and e_pcb are plugged in
        ind = pcblist==e_pcb & pinlist==e_pin & drivelist==iD;
        
        if sum(ind==1) %electrode connected
            
            experiment.hardware.microdrive(iD).electrodes(iE).channelid = acqChlist(ind);
            experiment.hardware.microdrive(iD).electrodes(iE).multiunitid = acqChlist(ind) + muidoffset;

            %only assign if gains for headstages have been defined
            if reassignGains
                experiment.hardware.microdrive(iD).electrodes(iE).gain = gainlist(ind);
            end
            
        elseif sum(ind==0) %electrode not connected
            experiment.hardware.microdrive(iD).electrodes(iE).channelid = nan;
            experiment.hardware.microdrive(iD).electrodes(iE).multiunitid = nan;
        
        else %found multiple channels for a given electrode
            error('problem with channel mapping')
        end
    end
end

