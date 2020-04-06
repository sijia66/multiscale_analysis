function P_chamber = drmap_coordTransform_microdrive2chamber(experiment, driveID)

% P_chamber = drmap_coordTransform_microdrive2chamber(experiment, driveID
%function that takes translates electrode positions from microdrive-based
%coordinates (defined in
%experiment.hardware.microdrive.electrodes.position) into chamber-based
%coordinates. Uses experiment.hardware.microdrive.coordinate for
%transformation of 1) a translation (by [.coordinate.x .coordinate.y] and
%2) a rotation of coordinate system by [.coordinate.theta] where theta is
%defined as the angle relative to the +x axis (+theta = CCW)
%
%input: experiment - experiment structure
%       driveID    - index of drive to return. Only supports one drive at a time.
%                    driveID can be a number index or the name of the microdrive (string)
%output: P_chamber - matrix (2 x #electrodes) of coordinates for each drive
%                    electrode in chamber-based coordinates. 
%
%A.Orsborn, 1/12/2018



if ischar(driveID) %if driveID is a string, find # of drive w/ matching name
    drive_names = {experiment.hardware.microdrive(:).name};
    drive_num   = find( ismember(drive_names, driveID) );
    
    if ~isempty(drive_num)
        driveID = drive_num;
    else
        driveID %#ok<NOPRT>
        drive_names %#ok<NOPRT>
        error('Drive name input not found in experiment definition file.')
    end
end
%otherwise, driveID should be a number. Do error checking. 
if length(driveID)>1
    error('drmap_CoordTransform_microdrive2chamber does not support multiple drive inputs. Use for one drive at a time.')
elseif driveID>length(experiment.hardware.microdrive)
    driveID %#ok<NOPRT>
    size(experiment.hardware.microdrive)
    error('driveID provided is outside the range of provided experiment definition file.')
end


%get microdrive-coordinate (x,y) position
pos = [experiment.hardware.microdrive(driveID).electrodes(:).position];
P = [pos(:).x; pos(:).y];

%get translation, rotation of microdrive relative to chamber from
%microdrive.coordinate
coord = [experiment.hardware.microdrive(driveID).coordinate];
offset = [coord.x; coord.y];
theta  = [coord.theta]; %rotation angle, defined as  angle relative to +x axis, with +theta --> counter clockwise


%translate microdrive to chamber-centered coordinates
P_shift = P + repmat(offset, [1 size(P,2)]);

%rotate electrode positions by theta
ROT = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
P_chamber = ROT*P_shift;
