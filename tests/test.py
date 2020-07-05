# for Travis
import numpy as np
import rcr

np.random.seed(18318) # get consistent random results

N = 1000 # total measurement count
frac_contaminated = 0.85 # fraction of sample that will be contaminated

# symmetric, uncontaminated distribution
mu = 0 
sigma_uncontaminated = 1
uncontaminated_samples = np.random.normal(mu, sigma_uncontaminated, 
    int(N * (1 - frac_contaminated)))

# one-sided contaminants
sigma_contaminated = 10
contaminated_samples = np.abs(np.random.normal(mu, sigma_contaminated, 
    int(N * frac_contaminated)))

# create whole dataset
data = np.concatenate((uncontaminated_samples, contaminated_samples))
np.random.shuffle(data)

# perform RCR
# initialize RCR with rejection technique:
# (chosen from shape of uncontaminated + contaminated distribution)
r = rcr.RCR(rcr.LS_MODE_68)
r.performBulkRejection(data) # perform outlier rejection

# View results
cleaned_data = r.result.cleanY
cleaned_mu = r.result.mu
cleaned_sigma = r.result.stDev

assert cleaned_sigma >= 0 