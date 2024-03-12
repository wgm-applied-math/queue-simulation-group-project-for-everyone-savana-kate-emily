% Script that runs a ServiceQueue simulation many times and plots a
% histogram

%% Set up

% Set up to run 100 samples of the queue.
n_samples = 100;

% Each sample is run up to a maximum time of 1000.
max_time = 1000;

% Record how many customers are in the system at the end of each sample.
NInSystemSamples = cell([1, n_samples]);
waitTimeSamples = cell([1, n_samples]);
serviceTimeSamples = cell([1, n_samples]);
totalTimeInSystemSamples = cell([1, n_samples]);

%% Run the queue simulation

% The statistics seem to come out a little weird if the log interval is too
% short, apparently because the log entries are not independent enough.  So
% the log interval should be long enough for several arrival and departure
% events happen.
for sample_num = 1:n_samples
    q = ServiceQueue(DepartureRate=1/1.5, DepartureRateWithHelper=1, LogInterval=100);
    % mu = 1/lambda = departure rate
    q.schedule_event(Arrival(1, Customer(1)));
    run_until(q, max_time);
    % Pull out samples of the number of customers in the queue system. Each
    % sample run of the queue results in a column of samples of customer
    % counts, because tables like q.Log allow easy extraction of whole
    % columns like this.
    NInSystemSamples{sample_num} = q.Log.NWaiting + q.Log.NInService;

    waitTime = zeros([length(q.Served)]);
    serviceTime = zeros([length(q.Served)]); 
    totalTimeInSystem = zeros([length(q.Served)]);

    for j = 1:length(q.Served)
        waitTime = q.Served{j}.BeginServiceTime - q.Served{j}.ArrivalTime;
    end
    waitTimeSamples{sample_num} = waitTime;

    for j = 1:length(q.Served)
        serviceTime = q.Served{j}.DepartureTime - q.Served{j}.BeginServiceTime;
    end
    serviceTimeSamples{sample_num} = serviceTime;

    for j = 1:length(q.Served)
        totalTimeInSystem = q.Served{j}.DepartureTime - q.Served{j}.ArrivalTime;
    end
    totalTimeInSystemSamples{sample_num} = totalTimeInSystem;

end

% Join all the samples. "vertcat" is short for "vertical concatenate",
% meaning it joins a bunch of arrays vertically, which in this case results
% in one tall column.

% use these variable names to call on histogram

% the probability of finding n customers in the system
NInSystem = vertcat(NInSystemSamples{:});
% the time customers spend waiting in the queue
waitTimeSum = vertcat(waitTimeSamples{:});
% the time customers spend being served
serviceTimeSum = vertcat(serviceTimeSamples{:});
% the total time customers spend in the system
totalTimeInSystemSum = vertcat(totalTimeInSystemSamples{:});

%% NInSystem
fig = figure();
t = tiledlayout(fig, 1,1);
ax = nexttile(t);
hold(ax, 'on');
h = histogram(ax, NInSystem, Normalization="probability");
title(ax, "Number of Customers in System");
xlabel(ax, "Time in System (minutes)");
ylabel(ax, "Probability");


% MATLAB-ism: When you pull multiple items from a cell array, the result is
% a "comma-separated list" rather than some kind of array.  Thus, the above
% means
%
%    NInSystem = horzcat(NInSystemSamples{1}, NInSystemSamples{2}, ...)
%
% which horizontally concatenates all the lists of numbers in
% NInSystemSamples.
%
% This is roughly equivalent to "splatting" in Python, which looks like
% f(*args).

%% waitTimeSum
fig = figure();
t = tiledlayout(fig, 1,1);
ax = nexttile(t);
hold(ax, 'on');
h = histogram(ax, waitTimeSum, Normalization="probability");
title(ax, "Time Customers Spend Waiting");
xlabel(ax, "Wait Time (minutes)");
ylabel(ax, "Probability");

%% serviceTimeSum
fig = figure();
t = tiledlayout(fig, 1,1);
ax = nexttile(t);
hold(ax, 'on');
h = histogram(ax, serviceTimeSum, Normalization="probability");
title(ax, "Time Customers Spend Being Served");
xlabel(ax, "Service Time (minutes)");
ylabel(ax, "Probability");

%% totalTimeInSystemSum
fig = figure();
t = tiledlayout(fig, 1,1);
ax = nexttile(t);
hold(ax, 'on');
h = histogram(ax, totalTimeInSystemSum, Normalization="probability");
title(ax, "Total Time Customers Spend in System");
xlabel(ax, "Time (minutes)");
ylabel(ax, "Probability");

%% Make a picture

% Start with a histogram.  The result is an empirical PDF, that is, the
% area of the bar at horizontal index n is proportional to the fraction of
% samples for which there were n customers in the system.
fig = figure();
t = tiledlayout(fig, 1, 1);
ax = nexttile(t);
h = histogram(ax, NInSystem, Normalization="probability", BinMethod="integers");

% MATLAB-ism: Once you've created a picture, you can use "hold on" to cause
% further plotting function to work with the same picture rather than
% create a new one.
hold(ax, "on");

% For comparison, plot the theoretical results for a M/M/1 queue.
% The agreement isn't all that good unless you run for a long time, say
% max_time = 10,000 units, and LogInterval is large, say 10.
rho1 = q.ArrivalRate / q.DepartureRate;
rho2 = q.ArrivalRate / q.DepartureRateWithHelper;

P0 = 1 - rho1;
P20 = 1/(1 + rho1 / (1-rho2));
nMax = 10;
ns = 0:nMax;
wh = 0:nMax;
P = zeros([1, nMax+1]);
P2 = zeros([1, nMax+1]);
P(1) = P0;
P2(1) = P20;

% red dots are WITHOUT helper
for n = 1:nMax
    P(1+n) = P0 * rho1^n;
end
plot(ax, ns, P, 'o', MarkerEdgeColor='k', MarkerFaceColor='r');

% blue dots are WITH helper
for n = 1:nMax
        P2(1+n) = P20 * rho1 * rho2^(n-1);
end
plot(ax, wh, P2, 'o', MarkerEdgeColor='k', MarkerFaceColor='b');

exportgraphics(fig, "histogram file name goes here.pdf");
title(ax, "Theoretical v.s. Actual Customers In System");
xlabel(ax, "Number of Customers in System");
ylabel(ax, "Probability");