function sc = uav_scenarios(name)

start = [0 0 0];
goal  = [100 100 100];

switch name
    case 'S1_sparse'
        obs = [
            20,25,15, 7;  50,20,30, 8;  78,28,18, 7;
            15,55,42, 8;  45,48,22, 7;  72,55,48, 8;
            30,78,28, 7;  62,74,62, 8;  86,70,32, 7;  40,33,62, 7];
        label = 'S1: Sparse obstacles (10)';

    case 'S2_medium'
        obs = [
            20,25,15, 7;  50,20,30, 8;  78,28,18, 7;  35,18,55, 7;
            15,55,42, 8;  45,48,22, 7;  72,55,48, 8;  58,42,70, 7;
            30,78,28, 7;  62,74,62, 8;  86,70,32, 7;  40,33,62, 7;
            12,35,75, 7;  88,45,65, 7;  50,85,45, 8;  25,60,12, 7];
        label = 'S2: Medium-density obstacles (16)';

    case 'S3_dense'
     
        obs = [
            15,15,10,6; 30,10,20,7; 45,20,5,5;  60,10,35,8; 85,15,15,6;
            10,35,25,7; 25,40,10,6; 40,35,40,9; 55,45,15,7; 80,35,30,8;
            90,50,10,6; 15,60,30,8; 35,65,50,10;55,55,35,9; 75,60,10,6;
            10,80,20,7; 30,85,40,8; 50,80,60,10;70,85,30,7; 90,75,50,9;
            20,95,10,5; 45,90,25,6; 65,95,45,8; 85,90,70,10;95,80,20,6];
        label = 'S3: Dense obstacles (25)';

    case 'S4_corridor'
       
        r = 12; zs = [10 30 50 70 90];
        obs = zeros(0,4);
        for z = zs
            for y = [0 24 48 96]
                obs(end+1,:) = [28, y, z, r]; 
            end
            for y = [0 48 72 96]
                obs(end+1,:) = [64, y, z, r]; 
            end
        end
        label = 'S4: S-shaped corridor (two gapped walls)';

    otherwise
        error('Unknown scenario：%s', name);
end

sc = struct('name',name,'label',label,'start',start,'goal',goal,'obs',obs);
end
