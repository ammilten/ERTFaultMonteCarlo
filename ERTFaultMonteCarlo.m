function simtime = ERTFaultMonteCarlo(fold,type,ns)
% fold: string with the path to the Monte Carlo output
% type: string, only works with 'object' right now
% ns:   number of samples to generate

addpath('rockphysics','priors','objectmodel','sampling','io','utils')

%% Load the Rest of the Parameters and Sample

run('ParameterFileNames.m')
run('GridParams.m')
% run('TProGSPrior.m')
run('GeneralMixingPrior.m')

switch type
    case 'object'
        run('FractureZonePrior.m')
end

PARAMS = samplePrior(FRACTURE,ns);
PRES = samplePrior(RESISTIVITY,ns);

setup_folders(fold);
save([fold,'Priors.mat'],'FRACTURE','GRID','RESISTIVITY')
save([fold,'Parameters.mat'],'PARAMS','PRES')


%% Make sure variables get passed to parallel pool
strike = strike;
dip = dip;
mofb = mofb;
tofb = tofb;
mivfb = mivfb;
tivfb = tivfb;
mcb = mcb;
tcb = tcb;
m = m;
t = t;
GRID = GRID;

%% Monte Carlo (will be gradually fixed)
tic
for i=1:ns
    disp(['Simulating ',num2str(i),' of ',num2str(ns)])
    switch type
        case 'object'
            [hw, fw] = FractureModel(PARAMS{i}); % hanging wall & foot wall
    end
    
    
end
simtime = toc;
disp(['Time to create ',num2str(ns),' lithology models: ',num2str(simtime/60),' minutes'])

