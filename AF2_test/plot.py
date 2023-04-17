import matplotlib.pyplot as plt
import seaborn as sns
import json
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Cambria'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2

# Read JSON data from a file
with open("ASN-NAG_reduced.json", "r") as file:
    data = json.load(file)

phi_values = [float(d["Phi"]) for d in data]
psi_values = [float(d["Psi"]) for d in data]
print(len(phi_values))

# Calculate standard deviations
phi_std = np.std(phi_values)
psi_std = np.std(psi_values)

# Calculate mean
phi_mean = np.mean(phi_values)
psi_mean = np.mean(psi_values)

# Plot KDE for Phi and Psi on the same axis
sns.kdeplot(phi_values, label=f"Phi [CG-ND2-C1-O5] (mean: {phi_mean:.2f}, std: {phi_std:.2f})", shade=True, color="blue", alpha=0.5)
sns.kdeplot(psi_values, label=f"Psi [CB-CG-ND2-C1] (mean: {psi_mean:.2f}, std: {psi_std:.2f})", shade=True, color="red", alpha=0.5)

plt.xlabel("Angle")
plt.ylabel("Density")
plt.title("ASN-NAG\nPhi and Psi Distributions (KDE)")
plt.legend()
plt.savefig("dis.png",dpi=450)
