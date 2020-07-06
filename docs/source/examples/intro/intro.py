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


# plot data
import matplotlib.pyplot as plt

ydata = np.random.uniform(0, 1, N) # project randomly into 2D for better visualization
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
ax.plot(data, ydata, "k.", label="Data pre-RCR", alpha=0.75, ms=4)
ax.plot(cleaned_data, ydata[r.result.indices], "bo", 
    label="Data post-RCR", alpha=0.4, ms=4)

# plot results
cont_mean = np.mean(data)
cont_sigma = np.std(data)

ax.axvspan(mu - sigma_uncontaminated, mu + sigma_uncontaminated, color='g', 
    alpha=0.25, label="1-$\sigma$ region of true\nuncontaminated distribution")
ax.axvline(x=cont_mean, c='r', lw=3, ls="-", alpha=0.75, 
    label="Pre-RCR sample mean of data")
ax.axvspan(cont_mean - cont_sigma, cont_mean + cont_sigma,  color='r', 
    fill=False, alpha=0.75, hatch="/", label="1-$\sigma$ region of data, pre-RCR")

ax.axvline(x=cleaned_mu, c='b', lw=3, ls="-", alpha=0.75, 
    label="RCR-recovered $\mu$ of\nuncontaminated distribution")
ax.axvspan(cleaned_mu - cleaned_sigma, cleaned_mu + cleaned_sigma, color='b', 
    fill=False, alpha=0.75, hatch= '\\', 
    label="1-$\sigma$ region of uncontaminated\ndistribution, after RCR")

plt.xlim(-5, 30)
plt.ylim(0, 1)
plt.xlabel("data")
plt.title("Results of RCR being used on an" + 
    " {}% contaminated dataset".format(frac_contaminated*100))
plt.yticks([])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()