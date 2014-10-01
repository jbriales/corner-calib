% path = fullfile( pwd, 
path = '/home/jesus/Corner_Datasets/GT1';
scans = loadLidarData('Rawlog', path, '_LASER_LMS');

XY = [scans.xy];
XY = reshape( XY, 2, 361, [] );
XY = mean( XY, 3 );

