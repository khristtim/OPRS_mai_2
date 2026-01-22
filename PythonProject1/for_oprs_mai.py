import pandas as pd
import matplotlib.pyplot as plt

def plot_orbit(path, title):
    df = pd.read_csv(path, header=None, names=["t","y1","y2","v1","v2"])
    plt.figure()
    plt.plot(df["y1"], df["y2"])
    plt.xlabel("y1")
    plt.ylabel("y2")
    plt.title(title)
    plt.axis("equal")
    plt.grid(True)

plot_orbit("D:/d/repositories/for_OPRS_mai_2/build/Desktop_Qt_6_10_0_MinGW_64_bit-Debug/arenstorf_1_dt001.csv", "Arenstorf #1, dt_out=0.01")
plot_orbit("D:/d/repositories/for_OPRS_mai_2/build/Desktop_Qt_6_10_0_MinGW_64_bit-Debug/arenstorf_1_dt1.csv",   "Arenstorf #1, dt_out=1.0 (dense output)")
plt.show()