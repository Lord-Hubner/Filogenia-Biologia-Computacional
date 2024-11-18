import numpy as np
import pandas as pd
import numpy as np

#Dados de input
data = [
    [0.0000, 0.1890, 0.1100, 0.1130, 0.2150],
    [0.1890, 0.0000, 0.1790, 0.1920, 0.2110],
    [0.1100, 0.1790, 0.0000, 0.0941, 0.2050],
    [0.1130, 0.1920, 0.0941, 0.0000, 0.2140],
    [0.2150, 0.2110, 0.2050, 0.2140, 0.0000]
]

labels = ["Gorila", "Orangotango", "Humano", "Chimpanzé", "Gibão"]

distanceMatrix = pd.DataFrame(data, index=labels, columns=labels)
nOTUs = len(distanceMatrix)
###################################################################
###################################################################
###################################################################


itemsTotalDistances = calculateTotalDistances(distanceMatrix, nOTUs)

