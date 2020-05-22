%% Program information %%

%This program calculates different features from a blood pressure curve.
%There is also functions for determining the quality of the signal although
%more work could be needed for improved results.

%Most of the functions are dependent on either cycleSetter or timeSetter.
%cycleSetter is often preferred since it gives a more accurate way of finding
%the cycle's "true" indices. 

%%
%Initializing variables

[time,signal,sample_freq,signalinfo] = rdmat('3326836_0002m');
period = 1/sample_freq;

plotATM('3326836_0002m');

%[time_noise,signal_noise, ~,signalinfo_noise] = rdmat('3326836_0002m_noise');
%plotATM('3326836_0002m_noise');

%[~,signal_noise2, ~,signalinfo_noise2] = rdmat('3326836_0002m_noise2');
%plotATM('3326836_0002m_noise2');

%  [time, signal, sample_freq, signalinfo] = rdmat('3399897_0009m');
%  period = 1/sample_freq;
%
%  plotATM('3399897_0009m');

plot(time, signal);
xlabel('Time (s)', 'FontSize', 20);
ylabel('Blood pressure (mmHg)', 'FontSize', 20)
title('Blood pressure vs. Time', 'FontSize', 20);
%%
%Script 
%Testing functions.

%[start, stop] = timeSetter();
[start, stop] = cycleSetter(period, signal);
%pr = find_pr(start, stop, period, signal);
%ppv = find_ppv(start, stop, period, signal)
%sbp = findSystolic(start, stop, period, signal)
%dbp = findDiastolic(start, stop, period, signal)                       
%MAP = findMAP(start, stop, period, signal)
%cycles = findCycles(period, period, signal, 0.7)

% t_sys = sys_duration(start, stop, period, signal)
% t_dia = dia_duration(start, stop, period, signal)
% t_sys_rise = sys_rise_duration(start, stop, period, signal)
% t_sys_decay = sys_decay_duration(start, stop, period, signal)
% t_decay = decay_duration(start, stop, period, signal)

%sys_area = sys_phase_area(start, stop, period, signal)
%dia_area = dia_phase_area(start, stop, time, period, signal)
%sys_rise_area = sys_rise_phase_area(start, stop, period, signal)
%sys_dec_area = sys_dec_phase_area(start, stop, period, signal)
%decay_area = decay_phase_area(start, stop, period, signal)

%MAP_plotter(start, stop, period, signal)

%signal_control(start, stop, period, signal, sample_freq); 

%[ACF, ACF_score] = autocorr_norm(start, stop, period, signal);
%subplot(2,1,1);
%[ACF, ACF_score] = autocor(start, stop, period, signal);
%subplot(2,1,2);
%[ACF_2, ACF_2_score, ~] = autocor_2(period, signal);

yf = transform(period, signal, sample_freq);

%%
%Stroke Volume %Not working correctly. Needs additional constant derived
%from multiple patients data. Could be re-written as function.

% [~, cycle_time] = find_pr(start, stop, period, signal);
% [~, dbp_index] = findMin(start, stop, period, signal);
% [~, d_notch_index] = findMin((dbp_index+6)*period, period*(dbp_index + (cycle_time/(period*2))), period, signal) ;
% bp_array = zeros(1, (d_notch_index - dbp_index));
% for i = 1 : (d_notch_index - dbp_index)
%     bp_array(1,i) = signal((i+dbp_index),1) - dbp;
% end
%
% dP_dt = diff(bp_array); %Differential
% sv_array = zeros(1, (d_notch_index - dbp_index) - 1);
%
% for i = 1 : ((d_notch_index - dbp_index) - 1) %One time less since dP_dt is 1 index smaller. The first index of bp_array is always 0.
%     sv_array(1,i) = dP_dt(1,i) * bp_array(1,i+1);
% end
%
% sv = trapz(sv_array) * ((d_notch_index - dbp_index) * period) %Integrating over time

%%
% MAP plotter, supposed to plot all the MAPs of the entire signal. 
% Due to changes in the calculation of pulse rate, this function does not 
% work correctly anymore. 

% function [MAP_array] = MAP_plotter(start, stop, period, signal)
% [~,cycle_time] = find_pr(start, stop, period, signal);
% cycles = findCycles(period, period, signal, cycle_time);
% cyclesCounter = 0;
% MAP_array = zeros(1,cycles);
% cycle_array = zeros(1,cycles);
% for i = 1 : cycles
%     %Adding 2 periods to avoid getting out of bounds
%     MAP_array(1,i) = findMAP(period, signal);
%     cycle_array(1,i) = i*cycle_time;
%     
%     cyclesCounter = cyclesCounter + 1;
% end
% 
% plot(MAP_array);
% end
%%
%Systolic blood pressure

function [sbp] = findSystolic(start, stop, period, signal)
sbp = findMax(start, stop, period, signal);
end

%%
%Diastolic blood pressure

function [dbp] = findDiastolic(start, stop, period, signal)
dbp = findMin(start, stop, period, signal);
end

%%
%Mean Arterial Pressure

function [MAP] = findMAP(start, stop, period, signal)
start_index = round(start/period); %Fixing problem with incorrect rounding in matlab
stop_index = round(stop/period);
bp_array = zeros(1, stop_index - start_index);
for i=1 : (stop_index - start_index)
    bp_array(1,i) = signal(i + start_index-1, 1);
end
MAP = mean(bp_array);

%More simplistic ways of calculating MAP, previously used before the
%current function was developed.
%MAP = (sys_phase_area(start, stop, period, signal) + dia_phase_area(start, stop, period, signal))/2;
%MAP = (2*findDiastolic(start, stop, period,signal) + findSystolic(start, stop, period, signal))/3;
end
%%
%Pulse rate
%Not applicable for all cycles, sensitive for false maxes in distorted
%signals.

function [pr,cycle_time] = find_pr(start, stop, period, signal)
[~,startindex] = findMax(start, stop, period, signal);
reference = 120; %Includes a cycle in array steps.
max_count = 1;
cycle_time = 0;
max_array = zeros(1,8); %8 positions just for buffer

for i=startindex:(startindex + reference)
    if (max_count == 1)
        max_array(1,1) = startindex;
        max_count = max_count + 1;
    else %3 step jumping to avoid plateaus
        if (signal(i+3,1) < signal(i,1) && signal(i-3,1) < signal(i,1))
            max_array(1, max_count) = i;
        end %To include plateaus, if maxes
        if ((signal(i+1,1)==signal(i,1) && signal(i-1,1)==signal(i,1)) && (signal(i+7,1) < signal(i,1) && signal(i-7,1) < signal(i,1)))
            max_array(1, max_count) = i;
        end
        if (max_array(1, max_count) - max_array(1, max_count-1) > 10) %Compares indices, not signal values
            max_count = max_count+1;
        end
    end
    if (max_count == 4 && cycle_time == 0) %To find the cycle time for the first cycle of interest
        cycle_time = (i-startindex)*period;
    end
end
pr = 60/cycle_time;
end
%%
%Pulse Pressure Variation
%Calculates the variation of all the pulse pressures in the signal.
%The function has only been tested for smaller time intervals, therefore,
%not suitable for extensive periods of time since the pulse rate probably
%can fluctuate too much.

function [ppv] = find_ppv(start, stop, period, signal)
[~,cycle_time] = find_pr(start, stop, period, signal);
cycles = findCycles(period, period, signal, cycle_time);
cyclesCounter = 0;
dbp_array = zeros(1,cycles);
sbp_array = zeros(1,cycles);
pp_array = zeros(1,cycles);
for i = 1 : cycles
    %Adding 2 periods to avoid getting out of bounds
    dbp_array(1,i) = findMin(cycle_time*(cyclesCounter + 1) - cycle_time + (2*period), cycle_time*(cyclesCounter+1), period, signal);
    sbp_array(1,i) = findMax(cycle_time*(cyclesCounter + 1), cycle_time*(cyclesCounter+1) + cycle_time, period, signal);
    pp_array(1,i) = sbp_array(1,i) - dbp_array(1,i);
    
    cyclesCounter = cyclesCounter + 1;
end
ppv = 100*(max(pp_array) - min(pp_array)) / mean(pp_array);
end

%%
%Time of the systolic phase, from diastole to dicrotic notch

function [t_sys, d_notch_index, dbp_index] = sys_duration(start, stop, period, signal)
[~, cycle_time] = find_pr(start, stop, period, signal);
[~, dbp_index] = findMin(start, stop, period, signal);
[~, d_notch_index] = findMin((dbp_index + 10)*period, (dbp_index*period + (cycle_time/2)), period, signal);
t_sys = (d_notch_index - dbp_index)*period;
end

%%
%Time of the diastolic phase, from dicrotic to next diastole.
%(Dependent on sys_duration).

function [t_dia, d_notch_index, dbp_index] = dia_duration(start, stop, period, signal)
[~,d_notch_index,~] = sys_duration(start, stop, period, signal);
[~, cycle_time] = find_pr(start, stop, period, signal);
[~, dbp_index] = findMin(d_notch_index*period, (d_notch_index*period + cycle_time), period, signal);
t_dia = (dbp_index - d_notch_index)*period;
end

%%
%Time of the systolic rise, from diastole to systole.

function [t_sys_rise, dbp_index, sbp_index] = sys_rise_duration(start, stop, period, signal)
[~, ~, dbp_index] = sys_duration(start, stop, period, signal);
[~, cycle_time] = find_pr(start, stop, period, signal);
[~, sbp_index] = findMax(dbp_index*period, (dbp_index*period + cycle_time), period, signal);
t_sys_rise = (sbp_index - dbp_index) * period;
end

%%
%Time of the systolic decay, from systole to dicrotic notch.

function [t_sys_decay] = sys_decay_duration(start, stop, period, signal)
t_sys_rise = sys_rise_duration(start, stop, period, signal);
[t_sys, ~, ~] = sys_duration(start, stop, period, signal);
t_sys_decay = (t_sys - t_sys_rise);
end

%%
%Time of the overall decay phase, from systole to next diastole.

function [t_decay] = decay_duration(start, stop, period, signal)
t_dia = dia_duration(start, stop, period, signal);
t_sys_decay = sys_decay_duration(start, stop, period, signal);
t_decay = t_dia + t_sys_decay;

end

%%
%Area from diastole to dicrotic notch.

function [sys_area] = sys_phase_area(start, stop, period, signal)
[~, d_notch_index, dbp_index] = sys_duration(start, stop, period, signal);
bp_array = zeros(1,(d_notch_index - dbp_index)); 
for i=1 : (d_notch_index - dbp_index)
    bp_array(1,i) = signal((i + dbp_index - 1),1); %-1 to include the whole interval
end
sys_area = mean(bp_array);
end

%%
%Area from the dicrotic notch to next diastole.

function [dia_area] = dia_phase_area(start, stop, period, signal)
[~, d_notch_index, dbp_index] = dia_duration(start, stop, period, signal);
bp_array = zeros(1, (dbp_index - d_notch_index)); 
for i=1 : (dbp_index - d_notch_index)
    bp_array(1,i) = signal((i + d_notch_index - 1),1); %-1 to include the whole interval
end
dia_area = mean(bp_array);
end

%%
%Area from diastole to systole.

function [sys_rise_area] = sys_rise_phase_area(start, stop, period, signal)
[~, dbp_index, sbp_index] = sys_rise_duration(start, stop, period, signal);
bp_array = zeros(1, (sbp_index - dbp_index+1)); 
for i=1 : (sbp_index - dbp_index+1)
    bp_array(1,i) = signal((i + dbp_index - 1),1); %-1 to include the whole interval
end
sys_rise_area = mean(bp_array);
end

%%
%Area from systole to dicrotic notch.

function [sys_dec_area] = sys_dec_phase_area(start, stop, period, signal)
[~, ~, sbp_index] = sys_rise_duration(start, stop, period, signal);
[~, d_notch_index, ~] = sys_duration(start, stop, period, signal);
bp_array = zeros(1, (d_notch_index - sbp_index)); 
for i=1 : (d_notch_index - sbp_index)
    bp_array(1,i) = signal((i + sbp_index - 1),1); %-1 to include the whole interval
end
sys_dec_area = mean(bp_array);
end

%%
%Area from systole to next diastole.

function [decay_area] = decay_phase_area(start, stop, period, signal)
[~, ~, sbp_index] = sys_rise_duration(start, stop, period, signal);
[~, ~, dbp_index] = dia_duration(start, stop, period, signal);
bp_array = zeros(1, dbp_index -sbp_index);
for i=1 : (dbp_index - sbp_index)
    bp_array(1,i) = signal((i + sbp_index -1),1); %-1 to include the whole interval
end
decay_area = mean(bp_array);
end

%%
%Time setter
%Lets the user enter a start and stop time for analysis. Needs to be a
%factor of "period".

function [start, stop] = timeSetter()
correct = 0;
while(correct == 0)
    start = input('\nStart time:');
    stop = input('End time:');
    if(start < stop && start > 0)
        correct = 1;
    else
        fprintf("\nInvalid inputs");
    end
end

end

%%
%Cycle setter
%Lets the user enter a cycle for analysis. May skip the first cycle of the
%signal due to robustness.

function [start, stop] = cycleSetter(period, signal)
correct = 0;
[~, cycle_time] = find_pr(period, period+1, period, signal); %Assuming that PR is roughly the same throughout the record
cycles = findCycles(period, period, signal, cycle_time);
while(correct == 0)
    cycle = input('\nCycle:');
    if(cycle > 0 && cycle <= cycles)
        correct = 1;
    else
        fprintf("\nInvalid input");
    end
end
start_approx = cycle_time*cycle;
start_approx2 = cycle_time*cycle + cycle_time*0.7;
[magnitude, start] = findMin(start_approx - 10*period, start_approx + 10*period, period, signal);
[magnitude2, start2] = findMin(start_approx2 - 10*period, start_approx2 + 10*period, period, signal);
if(magnitude2 < magnitude)
    start = start2;
end
start = start*period; %Converting to time
stop_approx = start + cycle_time;
[~, stop] = findMin(stop_approx - 10*period, stop_approx + 10*period, period, signal);
stop = stop*period; %Converting to time

end

%%
%Cycles of interest
%Lets the user enter multiple cycles for analysis.

function[interval, noOfCycles] = cycles_of_interest(period, signal)
start = 0; stop = 0;
while(start >= stop)
[start, ~] = cycleSetter(period, signal);
[~, stop] = cycleSetter(period, signal);
end
[~, cycle_time] = find_pr(start, stop, period, signal);
noOfCycles = round((stop-start)/cycle_time);

start = round(start/period); %Converting + fixing problem with incorrect rounding in matlab
stop = round(stop/period);
interval = zeros(1, stop-start);
for i=1 : stop-start
    interval(1,i) = signal(i-1+start,1);
end
end
%%
%Finding the maximum value for an interval.
function [max,index] = findMax(start, stop, period, signal)
max = 0;
index = 0;
loop_starter = fix(start/period); %Rounds downwards to nearest integer
loop_stopper = fix(stop/period);
for i=(loop_starter):(loop_stopper)
    if(signal(i,1) > max)
        max = signal(i,1);
        index = i;
    end
end

end

%%
%Finding the minimum value for an interval.

function [min,index] = findMin(start, stop, period, signal)
min = findMax(start, stop, period, signal);
index = 0;
loop_starter = fix(start/period); %Rounds downwards to nearest integer
loop_stopper = fix(stop/period);
for i=(loop_starter):loop_stopper
    if(signal(i,1) < min)
        min = signal(i,1);
        index = i;
    end
end

end

%%
%Number of cycles
%Returns the number of complete cycles from the entire signal.

function [noOfCycles] = findCycles(max_start, period, signal, cycle_time)
noOfCycles = 1;
max_reference = 40*period; %The max value should be within following interval
cycle_reference = 120*period; %Includes a cycle in array steps
[~, max_index] = findMax(max_start, max_start + max_reference, period, signal);
[~,min_start_index] = findMin(max_index*period, max_index*period + cycle_reference, period, signal);
[elements,~] = size(signal);
min_array = zeros(1,1);

for i = min_start_index : (cycle_time/period) : (elements-rem(elements-min_start_index,(cycle_time/period))-(cycle_time/period))
    %Making sure the loop does not get "out of bounds"
    [~,min_array(1,noOfCycles+1)] = findMin((i+10)*period, (i*period)+cycle_time, period, signal);
    %i+10 to prohibit missing the next min, if greater than previous min
    
    if(min_array(1,noOfCycles+1) ~= min_array(1,noOfCycles))
        noOfCycles = noOfCycles + 1;
    end
end
noOfCycles = noOfCycles - 1; %Removing one cycle, since noOfCycles starts at 1 for the comparison to work

end


%% Signal processing %%

%%
%Checks if the signal is good enough to analyze. Returns 1 if the signal is
%sufficient, and 0 if the signal isn't good enough. Due to robustness it is 
%recommended to use timeSetter instead of cycleSetter. 
%When using autocorrelation for several cycles and implementing timeSetter
%the score might not work properly, since the number of cycles could be
%false.

function [controller] = signal_control(start, stop, period, signal, sample_freq)
controller = 1;
while(controller == 1)
if(findMax(start, stop, period, signal) > 200 || findMax(start, stop, period, signal) < 90)
    controller = 0;
    [~, index] = findMax(start, stop, period, signal);
    start_time = index*period - sys_rise_duration(start, stop, period, signal); %Finding start of broken cycle
    end_time = index*period + decay_duration(start, stop, period, signal); %Finding end of broken cycle
    output = ['Invalid cycle: ', num2str(start_time), '-', num2str(end_time), ' s '];
    disp(output);
    break;
end
if(findMin(start, stop, period, signal) < 30 || findMin(start, stop, period, signal) > 120)
    controller = 0;
    [~, index] = findMin(start, stop, period, signal);
    start_time = index*period;
    end_time = start_time + sys_duration(start, stop, period, signal) + dia_duration(start, stop, period, signal);
    output = ['Invalid cycle: ', num2str(start_time), '-', num2str(end_time), ' s '];
    disp(output);
    break;
end
if(findMAP(start, stop, period, signal) > 150 || findMAP(start, stop, period, signal) < 40)
    controller = 0;
    output = 'Invalid cycle';
    disp(output);
    break;
end

PR = find_pr(start, stop, period, signal);
for i = round(start/period) : round(stop/period) %0.25 and -0.12 seems to be reasonable considering the cycles we've analysed.
    if((signal(i+1,1) - signal(i,1)) / PR > 0.25 || (signal(i+1,1) - signal(i,1)) / PR < -0.12)
        controller = 0;
        output = ('Too steep slope. Wait for a better signal');
        disp(output);
        break;
    end
end
[~, ACF_score] = autocor(start, stop, period, signal); 
if(ACF_score < 21 || ACF_score > 30) %Seems to be reasonable considering the cycles we've analysed.
    controller = 0;
    output = ('Not a valid cycle');
    disp(output);
    break;
end
[~, ACF2_score, noOfCycles] = autocor_2(period, signal); 
if(ACF2_score < (20*noOfCycles) || ACF_score > 30*noOfCycles) %Seems to be reasonable considering the cycles we've analysed.
    controller = 0;
    output = ('Wait for a better signal');
    disp(output);
    break;
end
[~, energyRatio] = transform(period, signal, sample_freq);
if(energyRatio < 0.9 || energyRatio > 0.985) %Below 0.9: Not fixable due to irregular noise. Above 0.985: Barely improves the signal
    controller = 0;
    output = ('Wait for a better signal');
    disp(output);
    break;
end

end
end
%%
%Fourier transform
%Computes and plots the Fourier transform of the signal and returns a
%filtered signal.
%Filter is currently only robust for one cycle. (Enter same cycle twice when 
%giving input)

function [yf, energyRatio] = transform(period, signal, sample_freq)
fprintf("Choose your cycle (only one) for transformation");
[interval, noOfCycles] = cycles_of_interest(period, signal);
frequency_signal = (1/length(interval)) * fft(interval); 
freq_vector = (0:length(frequency_signal)-1)*sample_freq/length(frequency_signal);
energy = sum(abs(frequency_signal(1 : 10))) + sum(abs(frequency_signal(end-9 : end))); %10 indices seems optimal for the signals we have tried.
energy_tot = sum(abs(frequency_signal)); %The ratio is the important variable
energyRatio = energy / energy_tot;

b = complex(zeros(1, length(frequency_signal)-20), 0);
yf = [frequency_signal(1:10), b, frequency_signal(end-9:end)]; %Similar to a low pass filter

signal2 = length(interval)*ifft(yf);
x = (0:length(signal2)-1)*noOfCycles/length(signal2);
subplot(2,1,1);
plot(x, interval);
xlabel('Cycle', 'FontSize', 20);
ylabel('Blood pressure [mmHg]', 'FontSize', 20);
title('Unfiltered signal', 'FontSize', 20);
subplot(2,1,2);

plot(x, signal2)
xlabel('Cycle', 'FontSize', 20);
ylabel('Blood pressure [mmHg]', 'FontSize', 20);
title('Filtered signal', 'FontSize', 20);

figure(2);
plot(freq_vector, abs(frequency_signal), 'r')
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Blood pressure [mmHg]', 'FontSize', 20);
title('Unfiltered Fourier Spectrum', 'FontSize', 20);

mean(interval); %These should be same. (MAP)
mean(signal2);
yf(1,1);

end

%%
%Autocorrelation function (within one cycle)
%Be aware, if you are using timeSetter more than one cycle could be
%included in the correlation and cause eventual problems.

function [ACF, ACF_score] = autocor(start, stop, period, signal)

[~,start_index] = findMin(start, stop, period, signal);
stop_index = start_index + (sys_duration(start, stop, period, signal)/period + dia_duration(start, stop, period, signal)/period);
cycle_of_intrest = zeros(1, (stop_index-start_index));
for i = 1 : stop_index-start_index
    cycle_of_intrest(1,i) = signal(start_index+i-1,1);
end
duration = stop_index-start_index-1;
y = cycle_of_intrest;
y_avg = mean(y);
y_std = std(y);

ACF = zeros(1, duration);
for i = 1 : duration
    y2 = circshift(y, i-1);
    ACF(1,i) = sum(((y-y_avg).*(y2-y_avg)) / (y_std*y_std*duration));
end

lags = 0 : duration-1;
stem(lags, ACF);
xlabel('Lags', 'FontSize', 20);
ylabel('Autocorrelation coefficient', 'FontSize', 20);
title('Autocorrelation within one cycle', 'FontSize', 20);

ACF_score = sum(abs(ACF));
end

%%
%Autocorrelation function (multiple cycles)

function [ACF_2, ACF_2_score, noOfCycles] = autocor_2(period, signal)
fprintf("Choose your cycles for autocorrelation");
[interval, noOfCycles] = cycles_of_interest(period, signal);
duration = length(interval);
y = interval;
y_avg = mean(y);
y_std = std(y);

ACF_2 = zeros(1, duration);

for i = 1 : duration
    y2 = circshift(y, i-1);
    ACF_2(1, i) = sum(((y-y_avg).*(y2-y_avg)) / (y_std*y_std*duration));
end
lags = 0 : duration-1;
subplot(2,1,1)
stem(lags, ACF_2);
xlabel('Lags');
ylabel('Autocorrelation coefficient');

x = noOfCycles*(0 : duration-1)/duration;
subplot(2,1,2)
stem(x, ACF_2);
xlabel('Cycles');
ylabel('Autocorrelation coefficient');

ACF_2_score = sum(abs(ACF_2));
end

%%
%Normaliezed autocorrelation function (within one cycle)
%Be aware, if you are using timeSetter more than one cycle could be
%included in the correlation and cause eventual problems.

function [ACF, ACF_score] = autocorr_norm(start, stop, period, signal)
start_index = round(start/period); %Fixing matlab rounding error
stop_index = round(stop/period);
cycle_of_intrest = zeros(1, (stop_index-start_index));
for i = 1 : stop_index-start_index
    cycle_of_intrest(1,i) = signal(start_index+i,1);
end
duration = stop_index-start_index-1;
y = cycle_of_intrest;
y_avg = mean(y);
y_std = std(y);

ACF = zeros(1, duration);
for i = 1 : duration
    y2 = circshift(y, i-1);
    ACF(1,i) = sum(((y-y_avg).*(y2-y_avg)) / (y_std*y_std*duration));
end

x = (0 : duration-1)/duration;
stem(x, ACF);
xlabel('Cycle');
ylabel('Autocorrelation coefficient');

ACF_score = sum(abs(ACF));
end
