function [sector_losses] = CalcSectoralLosses(event_data, daytimemult, industrymult)

% each outage event has a loss depending on day/time/season, industry, and
% length of event
% total loss to each sector is sum over all events

[T,~] = size(event_data);
[S,L] = size(industrymult);

sector_losses = zeros(S,1);

for t=1:T
    for l=1:L
        for s=1:S
            sector_losses(s) = sector_losses(s) + sum(event_data{t,l} .* daytimemult(t) .* industrymult(s,l));
        end
    end
end




end

