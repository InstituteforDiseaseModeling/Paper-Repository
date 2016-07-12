% function to update rate
function  [num_out, denom_out, counter_out] = BD_solveOnce_gamma(x0, t0, tf, k, ies)
ielen        = length(ies);
counter_out  = zeros(1,ielen);
num_out      = zeros(ielen, 2);
denom_out    = zeros(ielen, 2);

x      = x0;
t      = t0;
n      = zeros(1,2);
lambda = zeros(1,2);

% ies has the rarest one in the first element
ie_ind  = ielen;
while(t < tf)
	% check for a rare event
	if (x == ies(ie_ind))
		counter_out(ie_ind) = 1;
	
        num_out(ie_ind,:)   = n;
        denom_out(ie_ind,:) = lambda;
		
		if ie_ind == 1
			break;
		else
			ie_ind = ie_ind - 1;
		end
	end
    
    % compute propensity
    a  = [k(1) k(2)*x];
    a0 = a(1) + a(2);

	% time to the next reaction
	tau = log(1/rand) / a0;
	t   = t + tau;
	if t < tf
        if (a(1)/a0) > rand
            x = x + 1;
            n(1) = n(1) + 1;
        else
            x = x - 1;
            n(2) = n(2) + 1;
        end
		lambda = lambda + a*tau;
	end
end
