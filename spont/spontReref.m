function [data_reref, reref_pairs] = reref(data, type, pos, params)

nC = size(data,1);
if size(pos,1) ~= nC
    if size(pos,2) == nC
        pos = pos';
    else
        error('pos and data must have same # electrodes.')
    end
end


switch lower(type)
    
    case 'grandavg'
        
        %re-reference all electrodes to grand average
        m = mean(data,1);
        d = ndims(data);
        data_reref = data - permute( repmat(m, [ones(d-1,1) nC]), [d 1:(d-1)] );
        
        %set reref_pairs
        reref_pairs = zeros(nC,1);
        
    case 'nn_mean'
                
        %rereference each electrode to the average of it's nearest
        %neighbors
        
        data_reref = data;
        reref_pairs = nan(nC,4);
        for iC = 1:nC
            
            dist = sqrt( sum( (pos - repmat(pos(iC,:),[nC 1])).^2,2) );
            dist(dist==0) = nan;
            nn   = dist==min(dist);
            
            m_nn = mean(data(nn,:,:,:,:),1);
            data_reref(iC,:,:,:,:) = data(iC,:,:,:,:) - m_nn;
            
            reref_pairs(iC,1:sum(nn)) = find(nn==1);
        end
        
    case 'nn_rand'
        
        %rereference each electrode to one of the nearest neighbors
        %(randomly chosen)
        
        data_reref = data;
        reref_pairs = nan(nC,1);
        for iC = 1:nC
            
            dist = sqrt( sum( (pos - pos(iC,:)).^2,2) );
            nn   = find(dist==min(dist));
            
            pick = ceil(rand*length(nn)+1);
            data_reref(iC,:) = data(iC,:,:,:,:) - data(nn(pick),:,:,:,:);
            
            reref_pairs(iC,1) = nn(pick);
        end
        
    case 'radialdist'
        
        %rereference each electrode to average of electrodes w/in a certain
        %distance (specified in 'params')
        data_reref = data;
        for iC = 1:nC
            
            dist = sqrt( sum( (pos - pos(iC,:)).^2,2) );
            nn   = dist<=params.dist;
            
            m_nn = mean(data(nn,:,:,:,:),1);
            data_reref(iC,:) = data(iC,:,:,:,:) - m_nn;
            
            reref_pairs(iC,1:sum(nn)) = find(nn==1);
            
        end
        
        reref_pairs(reref_pairs==0) = nan;
        
    otherwise
        error('Type not supported')
end

        
        
        
        