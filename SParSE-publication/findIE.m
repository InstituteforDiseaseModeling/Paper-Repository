function [ies, prcs, rhois] = findIE(localExt, Kce, distance, counter, event_type)

counter_thrh = 100;

% if sign(distance) == event_type
switch event_type
	case -1
		localExt  = sort(localExt, 'ascend');
	case 1
		localExt = sort(localExt, 'descend');
	otherwise
		error('unknown rare event type. exiting...');
end

if sign(distance)==event_type
    if abs(distance) > 0.4 || counter < counter_thrh
        rhos = [0.005 0.01];
    elseif abs(distance) > 0.2
        rhos = [0.01 0.05 0.1];
    else
        rhos = [0.05 0.1 0.2];
    end 
else
    if abs(distance) > 0.4 || counter < counter_thrh
        rhos = [0.01 0.015];
    elseif abs(distance) > 0.2
        rhos = [0.05 0.1 0.15];
    else
        rhos = [0.1 0.15 0.2];
    end 
end

ie_inds = ceil(Kce * rhos);
ies  = localExt(ie_inds);
prcs = zeros(1,length(ies));

for i=1:length(ies)
	if event_type == -1
		prcs(i) = sum(localExt<=ies(i));
	else
		prcs(i) = sum(localExt>=ies(i));
	end
end

[c, ia, ~]   = unique(ies, 'stable');
ies           = c;
prcs          = prcs(ia);
rhois         = rhos(ia);