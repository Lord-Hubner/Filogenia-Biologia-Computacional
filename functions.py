import numpy as np
import pandas as pd
import numpy as np
import re

class NeighborJoiningTree:
    _distanceMatrix:pd.DataFrame
    _nOTUs:int
    _remainingOTUs:int
    _tree:str
    _addedNodesCounter=0
    _unitedLabelsDict={}

    def __init__(self, distanceMatrix, nOTUs) -> None:
        self._distanceMatrix = distanceMatrix
        self._nOTUs = nOTUs
        self._remainingOTUs = nOTUs
        self._tree = ""

    def __calculateTotalDistances(self) -> list:
        itemsTotalDistances = np.empty(self._remainingOTUs, dtype=float)

        for i in range(len(itemsTotalDistances)):
            currentRow = self._distanceMatrix.iloc[i, :]
            itemsTotalDistances[i] = sum(currentRow)/(self._remainingOTUs-2)

        return itemsTotalDistances

    def __getSmallestDistancePair(self, itemsTotalDistances:int) -> list:
        #(math.factorial(self._nOTUs)/math.factorial(self._nOTUs-2)*math.factorial(2)) #Combinação de nOTUs elementos dois a dois
        distanceMatrix:pd.DataFrame = self._distanceMatrix
        itemsSmallestValues = np.empty([self._remainingOTUs,self._remainingOTUs], dtype=float)
        currentSmallestDistance = [0, 0]

        for i in range(self._remainingOTUs):
            if(i == 0): continue;

            currentFrame = distanceMatrix.iloc[i, 0:i];
            print(currentFrame)
            
            for j in range(currentFrame.count()):
                currentValue = (distanceMatrix.iloc[i,j] - itemsTotalDistances[i] - itemsTotalDistances[j])
                itemsSmallestValues[i,j] = currentValue
                if currentValue < currentSmallestDistance[0]:
                    currentSmallestDistance = [currentValue, [i,j]]

    
        return currentSmallestDistance
    
    def __calculateBranchLengths(self, smallestDistance:list, itemsTotalDistances:list) -> list:
        #Fórmula: Va = D_AB/2 + (ua - ub)/2
        aTotalDistance = itemsTotalDistances[smallestDistance[1][0]]
        bTotalDistance = itemsTotalDistances[smallestDistance[1][1]]

        aBranchLength =  round(self._distanceMatrix.iloc[smallestDistance[1][0], smallestDistance[1][1]] + (aTotalDistance - bTotalDistance)/2, 4)
        bBranchLength = round(self._distanceMatrix.iloc[smallestDistance[1][0], smallestDistance[1][1]] + (bTotalDistance - aTotalDistance)/2, 4)

        return [[self._distanceMatrix.index[smallestDistance[1][0]], aBranchLength], [self._distanceMatrix.index[smallestDistance[1][1]], bBranchLength]]
    
    def __addBranchesToTree(self, branches:list) -> list:

        labelsToFind = np.empty([2], dtype=list)

        for i in range(len(branches)):
            if branches[i][0][0] == 'U':
                unitedLabels = self._unitedLabelsDict[branches[i][0]]
                labelsToFind[i] = unitedLabels.split('+')

        nodeAPosition = self._tree.find(branches[0][0]) if not labelsToFind[0] else self._tree.find(labelsToFind[0][-1])
        nodeBPosition = self._tree.find(branches[1][0]) if not labelsToFind[1] else self._tree.find(labelsToFind[1][-1])
        if (nodeAPosition != -1):
            self.__addNodeHasPreviousBranch(branches, nodeAPosition, 1) #Nó A repetiu

        elif (nodeBPosition != -1):
            self.__addNodeHasPreviousBranch(branches, nodeBPosition, 0) #Nó B repetiu

        else:
            self.__addNodeNoPreviousBranch(branches)

        return [branches[0][0], branches[1][0]]

    def __addNodeNoPreviousBranch(self, branches:list) -> None:
        if (self._tree == ""):
            self._tree = f"""({branches[0][0]}:{branches[0][1]},{branches[1][0]}:{branches[1][1]})"""
            return

        self._tree = f"""{self._tree},({branches[0][0]}:{branches[0][1]},{branches[1][0]}:{branches[1][1]})"""

    def __addNodeHasPreviousBranch(self, branches:list, nodePosition, newBranch:int) -> None:
        positionToAdd, subTreeBeginning = self.__findPositionToAddClade(nodePosition)
        self._tree = f"""{self._tree[:subTreeBeginning]}({self._tree[subTreeBeginning:positionToAdd]}):{branches[~newBranch][1]},{branches[newBranch][0]}:{branches[newBranch][1]})"""

    def __findPositionToAddClade(self, nodePosition):
        unclosedParenthesis = 0
        for i in range(len(self._tree[:nodePosition+1])): #Descobrir equilíbrio de parantesis na posição do nodo repetido
            if self._tree[i] == '(':
                unclosedParenthesis += 1
            elif self._tree[i] == ')':
                unclosedParenthesis -= 1

        while unclosedParenthesis != 0:
            i += 1
            if self._tree[i] == '(':
                unclosedParenthesis += 1
            elif self._tree[i] == ')':
                unclosedParenthesis -= 1

        for k in range(nodePosition, 0, -1): #Encontrar parêntesis inicial da sub-árvore
            if self._tree[k] == '(':
                break

        return i, k #posição onde termina a sub-árvore a ser juntada em um clado + parêntesis inicial da sub-árvore
    
    def __updateDistanceMatrix(self, unitedLabels:list):
        print(self._distanceMatrix.index)
        labelsToCalculateDistances = [i for i in self._distanceMatrix.index if i not in unitedLabels]
        newMatrix = self._distanceMatrix.drop(columns=[unitedLabels[0], unitedLabels[1]], index=[unitedLabels[0], unitedLabels[1]]) #ainda sem o nó novo
        print(newMatrix)
        newMatrix = self.__calculateDistancesToNewNode(unitedLabels, newMatrix, labelsToCalculateDistances)

    def __calculateDistancesToNewNode(self, unitedLabels:list, newMatrix:pd.DataFrame, labelsToCalculateDistances:list):
        labelsDistances = list()
        distanceAToB = self._distanceMatrix.loc[unitedLabels[0], unitedLabels[1]]
        for label in labelsToCalculateDistances:
            distanceAtoNew = self._distanceMatrix.loc[unitedLabels[0], label]
            distanceBToNew = self._distanceMatrix.loc[unitedLabels[1], label] 
            labelsDistances.append([label, (distanceAtoNew+distanceBToNew-distanceAToB)/2])

        self.__addNewNode(unitedLabels, labelsDistances, newMatrix)


    def __addNewNode(self, unitedLabels:list, labelsDistances:list, newMatrix: pd.DataFrame):
        self._addedNodesCounter += 1
        newNodeLabel = f"U{self._addedNodesCounter}"
        self._unitedLabelsDict[newNodeLabel] = f"{unitedLabels[0]}+{unitedLabels[1]}"

        newDistances = [i[1] for i in labelsDistances]; newDistances.append(0.0);
        labels = [i[0] for i in labelsDistances]; labels.append(newNodeLabel);

        newRow = pd.DataFrame({i[0] : i[1] for i in labelsDistances}, index=[newNodeLabel])
        newColumn = pd.DataFrame({newNodeLabel: newDistances}, index=labels)

        print(newColumn)
        print(newRow)

        self._distanceMatrix = pd.concat(objs=[newMatrix, newRow], axis=0)
        print(self._distanceMatrix)
        self._distanceMatrix = pd.concat(objs=[self._distanceMatrix, newColumn], axis=1)
        
        self._remainingOTUs -= 1
        print(self._distanceMatrix)

    def buildPhylogeneticTree(tree):

        while tree._remainingOTUs > 2:

            itemsTotalDistances = tree.__calculateTotalDistances()

            smallestDistance = tree.__getSmallestDistancePair(itemsTotalDistances)
            print(smallestDistance)

            branches = tree.__calculateBranchLengths(smallestDistance, itemsTotalDistances)
            unitedLabels = tree.__addBranchesToTree(branches)

            tree.__updateDistanceMatrix(unitedLabels)

        # tree._tree = f"""({tree._tree}):{round(tree._distanceMatrix.iloc[0,1], 4)}"""
        match = re.search("\)[^:]", tree._tree)
        finalDistance = round(tree._distanceMatrix.iloc[0,1], 4)
        tree._tree = f"""({tree._tree[:match.regs[0][0]+1]}:{finalDistance}{tree._tree[match.regs[0][0]+1:]}:{finalDistance})"""
        return tree._tree

            


