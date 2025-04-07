import pandas as pd
import matplotlib.pyplot as plt

# Load cell read counts from the file
file_path = "read1_read1_20_80_K562_PBMC_cell_counts.txt"
data = pd.read_csv(file_path, sep="\t")

# Sort the data by read counts in descending order
data_sorted = data.sort_values(by="ReadCount", ascending=False)

# Plot the knee plot with log scale for the x-axis
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(data_sorted) + 1), data_sorted["ReadCount"], marker="o", linestyle="-")
plt.xscale("log")  # Set x-axis to log scale
plt.yscale("log")  # Set y-axis to log scale
plt.xlabel("Cell Rank (Log Scale)", fontsize=14)
plt.ylabel("Read Count (Log Scale)", fontsize=14)
plt.title("Knee Plot of Cell Read Counts (Log Scale)", fontsize=16)
plt.grid(alpha=0.5)
plt.tight_layout()

# Save the plot to a file
output_plot_path = "20_80_K562_PBMC_counts_log_scale.png"
plt.savefig(output_plot_path, dpi=300)
plt.close()

print(f"Knee plot saved as {output_plot_path}")
