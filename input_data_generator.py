import numpy as np
import numpy.random as random

N = 10_000
OUTPUT_FILE = "gen_data.txt"
with open(OUTPUT_FILE, "w") as f:
    print(N, file=f)
    for i in range(N):
        m = 8810324116.227 
        start = random.normal(size=(3,))
        vel = random.normal(size=(3,))
        print(f"{m} {" ".join(map(str, start.tolist()))} {" ".join(map(str, vel.tolist()))}", file=f)