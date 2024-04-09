import matplotlib.pyplot as plt
import numpy as np

with open("res_rel.txt", "r") as data_file :
    data = data_file.readlines()

data = [x.split("\t") for x in data]
for x in data :
    x.remove("\n")

cols = len(data[0])
plot_data = np.zeros( (int(cols/2), len(data), 2) )

for i in range(len(data)) :
    for j in range(cols) :
        plot_data[int(np.floor(j/2))][i][j%2] = data[i][j]

fig, ax = plt.subplots(figsize=(10,10))

ax.set_xlabel("H[s]", fontsize = 25)
ax.set_ylabel("H[k]", fontsize = 25)
ax.tick_params(axis='both', which='major', labelsize=18)

for i in range(plot_data.shape[0]) :
    ax.scatter(plot_data[i,:,0], plot_data[i,:,1], label = f"N = {i}")

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=6, fancybox=True, shadow=True, fontsize=15)

plt.show()

