% Train Gauss
%generates a train of bellshaped pulses
function y = TrainGauss(SamplingTimes, PulseTimes, BellSpread)
    nPts = length(SamplingTimes);
    y = zeros(1,nPts);
        for i =1: length(PulseTimes)
            amp =1./(BellSpread.*sqrt(2.*pi))
            y (i+1,:) = y(i,:) + BellCurve(SamplingTimes, amp, PulseTimes(1,i), BellSpread)
        end
        