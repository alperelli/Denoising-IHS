% compute approximate divergence

function xvar = div_denoi(xhat, rhat, rvar, denoi)

%epsilon = obj.changeFactor*0.1*max(abs(rhat),[],1) + eps;
%epsilon = obj.changeFactor*0.1*mean(abs(rhat),1) + eps;
xhat = xhat(:); rhat = rhat(:);
dv = length(xhat);

epsilon = .1*min( sqrt(mean(rvar,1)), mean(abs(rhat),1) ) + eps;

n     = length(xhat);
n_avg = 2;

div = nan(n_avg,1);
for i = 1:n_avg % for each monte-carlo average
    eta = sign(randn(n,1)); % random direction
    rhat_perturbed = rhat + bsxfun(@times,epsilon,eta);

    xhat_perturbed = reshape(denoi(rhat_perturbed), dv, 1);
    div(i) = mean( eta.*(xhat_perturbed-xhat(:)) , 1)./epsilon;
end

div = mean(div); % max(min(div,obj.divMax),obj.divMin); % enforce limits

% compute posterior variance
xvar = rvar.*div
