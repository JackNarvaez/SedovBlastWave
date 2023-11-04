import matplotlib.pyplot as plt
import numpy as np

r, rho = np.loadtxt("results.txt", unpack=True)

fig, ax = plt.subplots(figsize=(7,5))
ax.plot(r, rho, "--", color="crimson")
ax.set_xlabel("r")
ax.set_ylabel(r"$\rho$")
plt.show()