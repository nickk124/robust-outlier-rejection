import numpy as np

np.random.seed(18318) # get consistent random results

N = 500                  # total measurement count
frac_contaminated = 0.5  # fraction of sample that will be contaminated

# symmetric, uncontaminated distribution
mu = 0 
sigma_uncontaminated = 1
uncontaminated_samples = np.random.normal(mu, sigma_uncontaminated, 
    int(N * (1 - frac_contaminated)))

# symmetric, contaminated distribution
sigma_contaminated = 5
contaminated_samples = np.random.normal(mu, sigma_contaminated, 
    int(N * frac_contaminated))

# create whole dataset
data = np.concatenate((uncontaminated_samples, contaminated_samples))
np.random.shuffle(data)

# plot data
import matplotlib.pyplot as plt

plt.figure(figsize=(8,5))
ax = plt.subplot(111)

ydata = np.random.uniform(0, 1, N) # project randomly into 2D for better visualization

ax.plot(contaminated_samples, ydata[:int(N * frac_contaminated)], "k.", 
    label="Pre-RCR dataset", alpha=0.75, ms=4)
ax.plot(uncontaminated_samples, ydata[int(N * frac_contaminated):], "k.", 
    alpha=0.75, ms=4)

plt.xlim(-15, 15)
plt.ylim(0, 1)
plt.xlabel("data")
plt.yticks([])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()

# get results pre-RCR
contaminated_mu = np.mean(data)
contaminated_sigma = np.std(data)
print(contaminated_mu, contaminated_sigma)

# perform RCR
import rcr

# initialize RCR with rejection technique:
# (chosen from shape of uncontaminated + contaminated distribution)
r = rcr.RCR(rcr.SS_MEDIAN_DL)
r.performBulkRejection(data) # perform outlier rejection

# View results post-RCR
cleaned_mu = r.result.mu
cleaned_sigma = r.result.stDev
print(cleaned_mu, cleaned_sigma)

# plot rejections
cleaned_data = r.result.cleanY

flags = r.result.flags 
# list of booleans corresponding to the original dataset, 
# true if the corresponding datapoint is not an outlier.

cleaned_data_indices = r.result.indices 
# indices of data in original dataset that are not outliers

plt.figure(figsize=(8,5))
ax = plt.subplot(111)
ax.plot(data[cleaned_data_indices], ydata[cleaned_data_indices], "b.", 
    label="RCR-accepted points", alpha=0.75, ms=4)

plt.xlim(-15, 15)
plt.ylim(0, 1)
plt.xlabel("data")
plt.yticks([])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()


# Weighting Data:

from scipy.stats import norm

def weight_data(datapoint):
    return norm.pdf(datapoint, loc=mu, scale=sigma_uncontaminated)

weights = weight_data(data)

# perform RCR
r = rcr.RCR(rcr.SS_MEDIAN_DL)
r.performBulkRejection(weights, data) # perform outlier rejection, now with weights


# View results post-RCR
cleaned_mu = r.result.mu
cleaned_sigma = r.result.stDev
print(cleaned_mu, cleaned_sigma)


# plot rejections
cleaned_data = r.result.cleanY
cleaned_data_indices = r.result.indices

plt.figure(figsize=(8,5))
ax = plt.subplot(111)
ax.plot(data[cleaned_data_indices], ydata[cleaned_data_indices], "b.", 
    label="RCR-accepted points,\nwith weights applied to data", alpha=0.75, ms=4)

plt.xlim(-15, 15)
plt.ylim(0, 1)
plt.xlabel("data")
plt.yticks([])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()