import os
import numpy as np

dirs = sorted(os.listdir())
print(dirs[1])
x = np.load(dirs[1])
print(x)
