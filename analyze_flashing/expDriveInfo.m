function [experiment, monkeydir, day_impDate_list, day_id_list, day_theta_list, day_x_list,day_y_list] = expDriveInfo(day, recs, drive_base, monkeydir)

for r = 1:length(recs)
    rec_drive = recs{r};
    try
        load([monkeydir '/' day '/' rec_drive '/rec' rec_drive '.experiment.mat'], 'experiment')
        
    catch
        warning(['No experiment definition file found for rec' rec_drive])
        return
    end
    
    drive_idx = find(contains({experiment.hardware.microdrive.name},drive_base));
    
    day_impDate_list(r) = {experiment.hardware.microdrive(drive_idx).implantDate};
    day_id_list(r)      = {experiment.hardware.microdrive(drive_idx).ID};
    day_theta_list(r)  = experiment.hardware.microdrive(drive_idx).coordinate.theta;
    day_x_list(r)       = experiment.hardware.microdrive(drive_idx).coordinate.x;
    day_y_list(r)       = experiment.hardware.microdrive(drive_idx).coordinate.y;
    
end
end