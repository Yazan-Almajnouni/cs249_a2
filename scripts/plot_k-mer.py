import pandas as pd
import matplotlib.pyplot as plt

# 1. Read the .hist file
df = pd.read_csv(
    "/ibex/user/majnouym/cs249_a2/lizard_asembly/merqury_2/lizard.hist",
    sep="\t",
    header=None,
    names=["coverage", "count"],
    comment="#",
    engine="python",
)

# 2. Sanity‐check shapes using df["count"]
print("coverage shape:", df["coverage"].shape)
print("   count shape:", df["count"].shape)
print(df.head())

# 3. Plot using df["count"]
plt.figure(figsize=(6, 4))
plt.plot(df["coverage"], df["count"], lw=1, color="steelblue")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("k-mer coverage")
plt.ylabel("Frequency")
plt.title("k-mer histogram")
plt.grid(True, which="both", ls="--", c="gray", alpha=0.3)
plt.tight_layout()
plt.savefig("kmer_histogram.png", dpi=300)
print("Saved kmer_histogram.png")
