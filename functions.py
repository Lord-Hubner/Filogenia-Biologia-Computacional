import numpy as np
import pandas as pd
import numpy as np
import math;


class NeighborJoiningTree:
    _distanceMatrix:pd.DataFrame
    _nOTUs:int

    def __init__(self, distanceMatrix, nOTUs) -> None:
        _distanceMatrix = distanceMatrix

    def calculateTotalDistances(self) -> list:
        itemsTotalDistances = np.empty(self._nOTUs, dtype=float)

        for i in range(len(itemsTotalDistances)):
            currentRow = self._distanceMatrix.iloc[i, :]
            print(currentRow)
            itemsTotalDistances[i] = sum(currentRow)/(self._nOTUs-2)

        return itemsTotalDistances

    def getSmallestDistancePair(self, itemsTotalDistances: int):
        #(math.factorial(self._nOTUs)/math.factorial(self._nOTUs-2)*math.factorial(2)) #Combinação de nOTUs elementos dois a dois
        distanceMatrix:pd.DataFrame = self._distanceMatrix

        for i in range(self._nOTUs):
            if(i == 0): continue;

            currentFrame = distanceMatrix.loc[i, 0:i];
            
            for j in range(currentFrame.columns):

            

            


