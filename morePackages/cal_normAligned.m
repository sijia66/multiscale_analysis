function mERP_norm = cal_normAligned(data)

    %z-score each channel
    m = mean(data,3);
    sd = std(data, [],3);
    data_norm = (data - repmat(m, [1 1 size(data,3)]))./ repmat(sd, [1 1 size(data,3)]);
    
    mERP = sq(mean(data,1));
    %subtract out dc differences across channels
    m = sq(mean(mERP,2));
    mERP = mERP - repmat(m, [1 size(mERP,2)]);
    
    mERP_norm = sq(mean(data_norm,1));
    
end 