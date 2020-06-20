import numpy as np
import rcr

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


# plot data
import matplotlib.pyplot as plt

ydata = np.random.uniform(0, 1, N)
plt.figure(figsize=(8,5))
ax = plt.subplot(111)
ax.plot(data, ydata, "k.", label="Data pre-RCR", alpha=0.75, ms=4)

plt.xlim(-5, 30)
plt.ylim(0, 1)
plt.xlabel("data")
plt.yticks([])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()