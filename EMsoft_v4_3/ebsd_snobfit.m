function [xbest, fbest, exitflag] = ebsd_snobfit(F, dx, lb, ub, options)
%EBSD_SNOBFIT
% Run SNOBFIT to optimize any defined function
% Make sure snobfit is added to the path before running
% INPUTS:
% F: function to minimize
% dx: resolution vector (minimum meaningful step size for each variable), n x 1 vector
% lb: lower bound for each variable, n x 1 vector
% ub: upper bound for each variable, n x 1 vector
% 'options': a struct containing the following fields:
% display: how much info to print to screen (notify, final, off, iter)[same as fminsearch]
% ncall = 100;   % limit on the number of function calls
% iterstop = 3;   % stop after this number of consecutive iterations without improvement in fbest
% npoint = n+6;   % number of random start points to be generated
% nreq = n+6;     % no. of points to be generated in each call to SNOBFIT
% p = 0.5;        % probability of generating random x to explore more of space 
%
% Adapted from snobtest.m from snobfit_v2.1
% Last updated: 7/24/19 (Edward Pang, MIT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get info from struct 'options'
display = options.display;
ncall = options.ncall;
iterstop = options.iterstop;
npoint = options.npoint;
nreq = options.nreq;
p = options.p;

% define some parameters
file = 'test';  % temp file to store intermediate data (after each call to SNOBFIT)
n = length(lb);     % figure out dimension of the problem


% computation of initial function values
x = rand(npoint,n)*diag(ub-lb) + ones(npoint,1)*lb'; % generation of npoint random starting points in [lb,ub]
for j=1:npoint
	f(j,:) = [F(x(j,:)) 0];
end

ncall0 = npoint;   % function call counter
fbestcounter = 0;   % initialize counter for # consecutive calls of snobfit without improvement of fbest
params = struct('bounds',{lb,ub},'nreq',nreq,'p',p); % input structure

exitflag = 0;   % exitflag=0 (max iter exceeded) unless program terminates early

%%% repeated calls to Snobfit
% print header
s = strcat('funevals',pad('[xbest]',length(lb)*10,'left'), '     fbest');
if strcmp(display, 'iter')
    fprintf('%s\n',s);    % header
end

while ncall0 < ncall+1 % repeat till ncall function values are reached
    % (if the stopping criterion is not fulfilled first)
	if ncall0 == npoint  % initial call
        [request,xbest,fbest] = snobfit(file,x,f,params,dx);

        % print results of initial iteration
        if strcmp(display, 'iter')
            fprintf('%8.0f',ncall0);    % more compact printing (Eddie edit)
            for ii=1:length(xbest)
                fprintf('%10.6f',xbest(ii));
            end
            fprintf('%10.6f\n',fbest);
        end
	else                 % continuation call
        [request,xbest,fbest] = snobfit(file,x,f,params);
    end

    clear x
    clear f
  
    % computation of the (perturbed) function values at the suggested points
	for j=1:size(request,1)
        x(j,:) = request(j,1:n);
        f(j,:) = [F(x(j,:)) 0];
    end
    
    ncall0 = ncall0 + size(f,1); % update function call counter
    [fbestn,jbest] = min(f(:,1)); % best function value from this round
    
    % update global min if better than any previous point
    if fbestn < fbest
        fbest = fbestn;
        xbest = x(jbest,:);
        fbestcounter = 0;   % reset counter
        
        % print
        if strcmp(display, 'iter')
            fprintf('%8.0f',ncall0);    % more compact printing (Eddie edit)
            for ii=1:length(xbest)
                fprintf('%10.6f',xbest(ii));
            end
            fprintf('%10.6f\n',fbest);
        end
    else
        fbestcounter = fbestcounter + 1;    % count number of consecutive calls of
                                            % snobfit where there is no improvement of fbest
    end
    
    % check stopping criterion (# of consecutive calls of snobfit without
    % improvement of fbest)
    if fbestcounter > iterstop-1
        exitflag = 1;   % objective not improving
        break
    end
end


%%% Print final results
if strcmp(display, 'final') || strcmp(display, 'iter')
    if exitflag == 0
        fprintf('\nExiting: Maximum number of iterations has been exceeded\n\n');
    else
        fprintf('\nOptimization terminated:\n objective not improving anymore\n\n');
    end

    fprintf('Final results:\n');
    fprintf('%s\n',s);    % header
    fprintf('%8.0f',ncall0);    % more compact printing (Eddie edit)
    for ii=1:length(xbest)
        fprintf('%10.6f',xbest(ii));
    end
    fprintf('%10.6f\n\n',fbest);
elseif strcmp(display, 'notify')
    % print only if does not converge
    if exitflag == 0
        fprintf('\nExiting: Maximum number of iterations has been exceeded\n\n');
    end
end



