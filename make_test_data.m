function testdata = make_test_data(scenario,n)
% This function makes test data to recover for the inversion for the
% different scenarios

% location of random sample
testdata.lat= repmat(30,n,1);
testdata.lon= repmat(10,n,1);
testdata.altitude=repmat(500,n,1);

switch scenario
    case 'step'
        testdata.t = [1500];                         % step change timing
        % testdata.e = [20,50,100,300,500,4000,5000,...    % old erosion rates of different catchments mm/ka
        %               20,50,100,300,500,4000,5000];% new erosion rates of catchments
        testdata.e = [20,50,100,300,50,400,50,...    % old erosion rates of different catchments mm/ka
              1000,2000,4000,200,50,100,500];% new erosion rates of catchments
        testdata.changeVariable = [];                % this is variable for convinience, so that I dont need a swith between loss/change or no change factor later on
    case 'samestep'
        testdata.t = [1500];                        % step change timing
        testdata.e = [20,50,100,40,60,80,90];       % old erosion rates of different catchments
        testdata.chg = [20];                        % change factor of erosion rate at time t

    case 'samebackground_step'
        testdata.t = [1500];                         % step change timing
        testdata.e = 50;                             % old erosion rates of different catchments
        testdata.chg = [10,20,5,40,2,25,30];       % change factor of erosion rate at time t
    case 'samebackground_samestep'
        testdata.t = 1000;
        testdata.e = 50;
        testdata.chg= 10;
    case 'spike'
        testdata.t =  [1500];                           % soil loss timing
        testdata.e =  [20,50,100,300,30,400,50];        % background erosion rates of different catchments
        testdata.loss =[50,20,1,  30, 10, 20, 40];       % soil loss in cm
    case 'samespike'
        testdata.t = [1500];                            % soil loss timing
        testdata.e = [20,50,100,300,50,400,50];         % background erosion rates of different catchments
        testdata.loss = [50];                            % soil loss in cm
    case 'samebackground_spike'
        testdata.t =   [1500];                          % soil loss timing
        testdata.e =   [50];                            % background erosion rates of different catchments
        testdata.loss = [50,10,100,30,20,50,5];         % soil loss in cm
    case 'samebackground_samespike'
        testdata.t =   [1500,500];                      % soil loss timing
        testdata.e =   [50];                            % background erosion rates of different catchments
        testdata.loss = [10,20];                        % soil loss in cm

    case 'curve'
        % apply pollen data as test
        curvedata = load('./data/pollen.mat');
        pollen = curvedata.pollen;

        timebreaks = [10000, 6200, 700, 0];
        for i = 1:length(timebreaks) - 1
            timeRange = pollen.yearsBP >= timebreaks(i+1) & pollen.yearsBP < timebreaks(i); % Define the time range for this period
            meanPercTree(i) = mean(pollen.percTree(timeRange)); % Calculate the mean percTree for this time range
        end
        noTreePerc = 100-meanPercTree; 
        curvechanges    = (noTreePerc(2:end) ./ noTreePerc(1) -1);

        scaleFactor     = [10, 1, 4, 5, 0.1, 100, 20];                           

        testdata.t = timebreaks(2:end-1);
        testdata.e = [10, 500, 100, 50, 20, 80, 150];   % background erosion rate - erosion rate at start of curve
        testdata.chg = scaleFactor;                     % this is the scaling factor of the pollen curve. E_increase = scaleFactor * curvechanges
        testdata.curvechange = curvechanges;            % these are the base changes that will be scaled later on
end

testdata.steps = length(testdata.t);
if or(isfield(testdata,'chg'), isfield(testdata,'curve'))
    testdata.changeVariable = testdata.chg;
elseif isfield(testdata,'loss')
    testdata.changeVariable = testdata.loss;
else
    testdata.changeVariable = [];
end


end