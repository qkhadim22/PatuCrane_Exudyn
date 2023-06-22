
#%%
import os
import numpy as np
import matplotlib.pyplot as plt
import csv



data = []
all_data = []

FILE_NAME = "ExData/LiftFlexible/Angle.txt"

with open(FILE_NAME, 'r') as f:
    data = list(csv.reader(f, delimiter="," ))
    data = np.array(data)
    
print(data.shape)

n_samples, n_features = data.shape
n_features -=1
X = data[:,0:n_features]
Y = data[:,n_features]

print(X.shape, Y.shape)
 

