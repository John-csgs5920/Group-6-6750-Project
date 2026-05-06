function simulate_A3_handover()



rng(5);                        %seed trials 
numUsers = 2000;
numCells = 8;
simTime = 200;                 % seconds
dt = 0.1;                     %Time step(s)
steps = round(simTime/dt);
%Hexagonal Layout
%CellRadius in meters.
cellRadius = 500;             
sitePositions = generate_cell_positions(numCells, cellRadius);
%TX power dbm
txPower_dBm = 46;              
freqGHz = 3.5;                 
bwHz = 100e6;                 
noiseFigure = 7;               
noiseTemp = 290;
kBoltz = 1.38064852e-23;

% Event A3 parameters (configurable)
hysteresis_dB = 3;           
margin_dB = 1;               % Offset
TTT =0.256;                   % TTT(s)


shadowStd_dB = 8;              
corrDist = 50;                 


minSpeed = 0.5; maxSpeed = 1.5;  %(m/s)
pauseProb = 0.1;

%RSRP pathloss 
plModel = @(d) tr38901_pathloss(d, freqGHz);

%Compute noise power.
noisePowerW = kBoltz * noiseTemp * bwHz * 10^(noiseFigure/10);
noisePower_dBm = 10*log10(noisePowerW*1000);

%%
users(numUsers) = struct();
for u=1:numUsers
    users(u).pos = [ (rand-0.5)*2*cellRadius, (rand-0.5)*2*cellRadius ]; 
    users(u).dest = rand_dest_in_area(cellRadius);
    users(u).speed = minSpeed + (maxSpeed-minSpeed)*rand;
    users(u).pause = 0;
    users(u).serving = randi(numCells);    %first serving cell
    users(u).rscp = -Inf;
    users(u).sinr = -Inf;
    users(u).handoverState = struct('candidate',0,'tttTimer',0);
    users(u).handoverCount = 0;
    users(u).handoverFailed = 0;
    users(u).totalThroughput = 0;
    users(u).totalLatency = 0;
    users(u).numPackets = 0;
end

gridResolution = 10; % meters
[xg, yg, shadowField] = generate_correlated_shadowing_field(sitePositions, cellRadius, gridResolution, shadowStd_dB, corrDist, numCells);

%Sim Records 
handoverEvents = [];
servingHistory = zeros(numUsers, steps);

%% Main simulation loop
for tstep = 1:steps
    t = (tstep-1)*dt;
    %User Movement.
    for u=1:numUsers
        users(u).pos = move_user(users(u), dt, cellRadius);
    end
    
    %Compute RSRP for each user 
    
    posMat = reshape([users.pos],[2,numUsers])'; 
    cellPosMat = sitePositions; 
    dists = pdist2(posMat, cellPosMat); 
    
    
    PL_dB = plModel(dists);
   
    shadow_dB = interp_shadow(shadowField, xg, yg, posMat, numCells);
    
   
    rsrp_dBm = txPower_dBm - PL_dB + shadow_dB;
    
    
    for u=1:numUsers
        
        s = users(u).serving;
        rsrp_serv = rsrp_dBm(u,s);
        
        
        [bestVal, bestCell] = max(rsrp_dBm(u,:));
        if bestCell==s
           
            users(u).handoverState.candidate = 0;
            users(u).handoverState.tttTimer = 0;
        else
            
            if bestVal - rsrp_serv > margin_dB + hysteresis_dB
                
                if users(u).handoverState.candidate ~= bestCell
                    users(u).handoverState.candidate = bestCell;
                    users(u).handoverState.tttTimer = 0;
                else
                    users(u).handoverState.tttTimer = users(u).handoverState.tttTimer + dt;
                    if users(u).handoverState.tttTimer >= TTT

                        users(u).handoverCount = users(u).handoverCount + 1;
                        success = attempt_handover(u, s, bestCell, rsrp_dBm, posMat, sitePositions, cellRadius);
                        if success
                            users(u).serving = bestCell;
                            users(u).handoverState.candidate = 0;
                            users(u).handoverState.tttTimer = 0;
                            handoverEvents = [handoverEvents; t, u, s, bestCell, 1]; %#ok<AGROW>
                        else
                            users(u).handoverFailed = users(u).handoverFailed + 1;
                            users(u).handoverState.candidate = 0;
                            users(u).handoverState.tttTimer = 0;
                            handoverEvents = [handoverEvents; t, u, s, bestCell, 0]; %#ok<AGROW>
                        end
                    end
                end
            else
                
                users(u).handoverState.candidate = 0;
                users(u).handoverState.tttTimer = 0;
            end
        end
        
        
        servingCell = users(u).serving;
        sinr_lin = db2pow(rsrp_dBm(u,servingCell) - noisePower_dBm);
        
    end
    
    usersPerCell = zeros(numCells,1);
    for c=1:numCells
        usersPerCell(c) = sum([users.serving]==c);
    end
   
    for u=1:numUsers
        c = users(u).serving;
        if usersPerCell(c)<=0
            allocBw = 0;
        else
            allocBw = bwHz / usersPerCell(c);
        end
        rsrp_u = rsrp_dBm(u,c);
        sinr_lin = db2pow(rsrp_u - noisePower_dBm);
        cap_bps = allocBw * log2(1 + sinr_lin);
        users(u).totalThroughput = users(u).totalThroughput + cap_bps * dt;
        pktSize = 1500*8; 
        if cap_bps>0
            estLat = pktSize / cap_bps; 
        else
            estLat = Inf;
        end
        users(u).totalLatency = users(u).totalLatency + estLat * dt;
        users(u).numPackets = users(u).numPackets + dt; % accumulate time for averaging later
    end
    
    
    servingHistory(:,tstep) = [users.serving]';
end

%% Display metrics.
totalHandovers = sum([users.handoverCount]);
totalFailures = sum([users.handoverFailed]);
handoverFailureRatio = totalFailures / max(totalHandovers,1);

throughputPerUser_bps = arrayfun(@(u) users(u).totalThroughput/simTime, 1:numUsers);
avgLatencyPerUser_s = arrayfun(@(u) users(u).totalLatency / users(u).numPackets, 1:numUsers);
avgThroughput_bps = mean(throughputPerUser_bps);
medianThroughput_bps = median(throughputPerUser_bps);


fprintf('Simulation complete.\n');
fprintf('Total handovers attempted: %d\n', totalHandovers);
fprintf('Handover failures: %d\n', totalFailures);
fprintf('Handover failure ratio: %.4f\n', handoverFailureRatio);
fprintf('Average throughput per user: %.2f kbps\n', avgThroughput_bps/1e3);
fprintf('Median throughput per user: %.2f kbps\n', medianThroughput_bps/1e3);
fprintf('Average latency per user: %.3f s\n', mean(avgLatencyPerUser_s(~isinf(avgLatencyPerUser_s))));


figure;
histogram(throughputPerUser_bps/1e3,50);
xlabel('Throughput per user (kbps)');
ylabel('Number of users');
title('Per-user throughput distribution');

figure;
imagesc(servingHistory(1:200,:)); 
xlabel('Time step');
ylabel('User index');
title('Serving cell over time (first 200 users)');
colorbar;

end

%Handover functions.

function pos = rand_dest_in_area(R)
    pos = [ (rand-0.5)*2*R, (rand-0.5)*2*R ];
end

function pos = move_user(user, dt, R)
    
    if user.pause>0
        user.pause = user.pause - dt;
        pos = user.pos;
        return;
    end
    dest = user.dest;
    vec = dest - user.pos;
    dist = norm(vec);
    if dist < 1
       
        user.dest = rand_dest_in_area(R);
        if rand < 0.1
            user.pause = 2 + 3*rand;
        end
        pos = user.pos;
        return;
    end
    dir = vec / dist;
    step = min(user.speed * dt, dist);
    pos = user.pos + dir * step;
end

function positions = generate_cell_positions(N, R)
    %place cells in network. 
    theta = (0:N-1)' * 2*pi/N;
    positions = [R*cos(theta), R*sin(theta)];
end

function [xg, yg, shadowField] = generate_correlated_shadowing_field(sitePositions, R, res, sigma, corrDist, numCells)
  
    margin = 2*R;
    gx = -R-margin:res:R+margin;
    gy = gx;
    [X,Y] = meshgrid(gx,gy);
    xg = gx; yg = gy;
    shadowField = zeros([size(X), numCells]);
 
    for c=1:numCells
        rng(100+c);
        sz = size(X);
       
        w = randn(sz);
       
        [kx, ky] = meshgrid( (-floor(sz(2)/2)):(ceil(sz(2)/2)-1), (-floor(sz(1)/2)):(ceil(sz(1)/2)-1) );
        dx = res;
        dist = sqrt((kx*dx).^2 + (ky*dx).^2);
        kernel = exp(-dist/corrDist);
        kernel = fftshift(kernel);
       
        field = real(ifft2( fft2(w) .* fft2(kernel) ));
       
        field = field - mean(field(:));
        field = field / std(field(:)) * sigma;
        shadowField(:,:,c) = field;
    end
end

function S = interp_shadow(shadowField, xg, yg, posMat, numCells)
    
    numUsers = size(posMat,1);
    S = zeros(numUsers, numCells);
    [Xg, Yg] = meshgrid(xg, yg);
    for c=1:numCells
        
        S(:,c) = interp2(Xg, Yg, shadowField(:,:,c), posMat(:,1), posMat(:,2), 'linear', 0);
    end
end

function PL = tr38901_pathloss(d, freqGHz)
    
    %TR 38.901 pathloss.
    d0 = 1;
    c = 3e8;
    lambda = c/(freqGHz*1e9);
    PLfs_d0 = 20*log10(4*pi*d0/lambda);
    PL = zeros(size(d));
    
    d(d<1) = 1;
    n = 3.5; 
    PL = PLfs_d0 + 10*n*log10(d./d0);
    
    PL(d>1000) = PL(d>1000) + 20;
end

function ok = attempt_handover(u, s, target, rsrp_dBm, posMat, sitePositions, cellRadius)
    
    diff = rsrp_dBm(u,target) - rsrp_dBm(u,s);
    
    if diff > 5
        p = 0.98;
    elseif diff > 2
        p = 0.9;
    elseif diff > 0
        p = 0.7;
    else
        p = 0.2;
    end
    ok = rand < p;
end

function val = db2pow(x)
    val = 10.^(x/10);
end