function [def] = drmap_layoutElectrodes(def, nD)

%def = drmap_layoutElectrodes(def)
%Reads a drive layout file and populates electrodes.position accordingly. 
%input: def - experiment definition structure (standard lab format assummed)
%       nD  - drive to build (optional. If not defined, loops through all
%       defined drives)
%output: def - updated experiment definition structure

if exist('nD', 'var')
    drives = nD;
else
    nD = size(def.hardware.microdrive,2);
    drives = 1:1:nD;
end

for d=drives
    
    %for increased readability...
    tmp = def.hardware.microdrive(d);
    
    %read in position information from def file
    [pos_ids, labels] = drmap_readMap(tmp.layout.file, 1, tmp.layout.in, tmp.layout.out);
    
    

    %if 'electrodes' already defined, check that things seem to match up.
    %Warn accordingly. 
    
    if isfield(tmp, 'electrodes')
        nE_def = size(tmp.electrodes,2);
        nE_layout = size(pos_ids,1);
        
        if ~isequal(nE_def, nE_layout)
            
            if nE_def> nE_layout
                fprintf('!!!!!!WARNING!!!!!! More electrodes appear to be defined in drive %s than are included in the layout file. Please check configurations.', tmp.name)
            else
                fprintf('!!!!!!WARNING!!!!!! More electrodes appear to be defined in layout file than are defined in drive %s. Please check configurations.', tmp.name)
            end
        end
    end
    
    
    %populate 'electrodes.position'
    nE = size(pos_ids,1);
    E = pos_ids(:,ismember(labels, 'electrode'));
    
    x_col = ~cellfun('isempty', strfind(labels, 'x'));
    y_col = ~cellfun('isempty', strfind(labels, 'y'));
    row_col = ~cellfun('isempty', strfind(labels, 'row'));
    col_col = ~cellfun('isempty', strfind(labels, 'col'));
    
       
    for e=1:nE
        
        tmp.electrodes(E(e)).position.x   = pos_ids(e,x_col);
        tmp.electrodes(E(e)).position.y   = pos_ids(e,y_col);
        tmp.electrodes(E(e)).position.row = pos_ids(e,row_col);
        tmp.electrodes(E(e)).position.col = pos_ids(e,col_col);
        
    end
    
    
    %update def
    def.hardware.microdrive(d).electrodes = tmp.electrodes;
    
end
