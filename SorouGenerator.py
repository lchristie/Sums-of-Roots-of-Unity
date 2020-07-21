import support
import pickle
import itertools
import ReferenceTypes as rt
import TypeGen
import pdb
import time
from multiprocessing import Pool as ThreadPool 

def genBaseSorou (aMinVanType):
    working = aMinVanType[0]
    baseSorou = []
    for index in range(working[0]):
        baseSorou.append([support.rotate(working[1], working[0], index)])
    return baseSorou

#### https://stackoverflow.com/questions/6284396/permutations-with-unique-values
class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

#### https://stackoverflow.com/questions/6284396/permutations-with-unique-values
def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

#### https://stackoverflow.com/questions/6284396/permutations-with-unique-values
def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1

def findSubTypePermutations (n, p, myPrevPerms):
    if n > p: return False

    if (n,p) in myPrevPerms.keys():
        return myPrevPerms[(n,p)]

    working = [0] * p

    if n == 0:
        mPerm = tuple(working)
        return set([mPerm])

    for index in range(n):
        working[index] = index + 1

    unfiltered = list(perm_unique(working[1:]))
    output = set()
    for item in unfiltered:
        output.add(tuple([1]+ list(item)))

    myPrevPerms[(n,p)] = output

    return output

def makeBitMap (p):
    output = []

    for index in range(2 ** p):
        bitsString = bin(index)[2:]
        bits = [0] * (p - len(bitsString))
        for char in bitsString:
            bits.append(int(char))
        output.append(bits)

    return output

def selectSubSorous (aBits, aNewSorou):
    output = []

    for index in range(len(aBits)):
        output += aNewSorou[index][aBits[index]]

    return output

def shouldAdd (aSorou, aSorouSet):
    for item in aSorouSet:
        if support.isEquiv(aSorou, item):
            return False

    return True

def toStr(n, base):
    # http://interactivepython.org/courselib/static/pythonds/
    #        Recursion/pythondsConvertinganIntegertoaStringinAnyBase.html
    convertString = "0123456789ABCDEFGHIJK"
    if n < base:
        return convertString[n]
    else:
        return toStr(n // base, base) + convertString[n % base]

def toBaseN (number, base, length):
    inBase = toStr(number, base)

    toAdd = length - len(inBase)

    strToAdd = '0' * toAdd

    return strToAdd + inBase

def makeMap (m, n):
    if m > 19:
        print('Error, f_0 too long')

    output = []

    temp = [0] * m

    for index in range(n ** m):
        indexBaseM = list(toBaseN(index, n, m))
        indexArray = [int(item) for item in indexBaseM]
        output.append(indexArray)

    return output

def partitionF0 (aF0, n):
    #### Partitions aF0 into n sub sorous
    #### Returns a list of all such paritions

    output = []

    mappings = makeMap(len(aF0), n)

    for mapping in mappings:
        temp = [[]] * n
        for index in range(len(aF0)):
            temp[mapping[index]] = temp[mapping[index]] + [aF0[index]]
        shouldAdd = True
        for item in temp:
            if item == []:
                shouldAdd = False
        if shouldAdd:
            output.append(temp)

    return output

def canBeSubSorou (aSorou, anotherSorou):
    """
    Returns true if anotherSorou can be rotated so that it is a subSorou of aSorou, and false otherwise
    """

    output = False

    aSorouRelOrder = support.findRelativeOrder(aSorou)
    anotherSorouRelORder = support.findRelativeOrder(anotherSorou)
    crossRelativeOrder = support.lcm(aSorouRelOrder, anotherSorouRelORder)

    for index in range(crossRelativeOrder):
        rotatedSorou = support.rotate(anotherSorou, crossRelativeOrder, index)
        if support.isSubSorou(aSorou, rotatedSorou):
            output = True
            break

    return output

def matchTypeToSorous (aF0Parition, aSorouListList):
    #### Output is a three dimensional array
    ##      Index 1 - The sorou in aF0Parition
    ##      Index 2 - The subtype index
    ##      Index 3 - The specific sorou of subtypes[I2]

    output = []
    flagList = []

    for sorou in aF0Parition:
        t2 = []
        for aSorouList in aSorouListList:
            t1 = []
            for aSorou in aSorouList:
                if canBeSubSorou(aSorou, sorou):
                    t1.append(True)
                else:
                    t1.append(False)
            t2.append(t1)
        output.append(t2)
    
    return output

def simplifyArray (a3DBoolArray):
    rowN = len(a3DBoolArray)
    colN = len(a3DBoolArray[0])

    output = []

    for rowIndex in range(rowN):
        output.append([])
        for colIndex in range(colN):
            if (True in a3DBoolArray[rowIndex][colIndex]):
                output[rowIndex].append(True)
            else:
                output[rowIndex].append(False)

    return output

def findRotationIndicies (aSubSorou, aF0):
    output = set()

    for root in aSubSorou:
        a = root[0] * aF0[0][0]
        b = (root[0] * aF0[0][1] - aF0[0][0] * root[1]) % a
        output.add((a, b))

    return output

def checkF0OddParity (aF0):
    if len(aF0) == 1:
        return True
    elif aF0[0][1] % 2 == 0:
        return False 
    else:
        return True

def selectFromNewSorou (aNewSorou):
    output = [[]]
    for item in aNewSorou:
        temp = []
        for subItem in item:
            for i in output:
                temp.append(i+subItem)
        output = temp

    return output

def genAllSorousOfMinVanType (aMinVanType, aTypeSorouMap, aPermutationsMap):
    aTypeStr = support.typeToLatexString(aMinVanType)

    if aTypeStr in aTypeSorouMap.keys() and aTypeSorouMap[aTypeStr] != []:
        return aTypeSorouMap[aTypeStr]

    outputKeys = set()
    output = []

    subTypes = aMinVanType[0][2:]

    subSorous = []
    for subType in subTypes:
        subSorou = genAllSorousOfType(subType, aMinVanType[0][1], aTypeSorouMap, aPermutationsMap)
        subSorous.append(subSorou)

    support.tPrint("Generating SubType Permutations")
    subTypePermunations = findSubTypePermutations(len(subTypes), aMinVanType[0][0], aPermutationsMap)
    support.tPrint("{} SubType Permutations Generated".format(len(subTypePermunations)))

    subSorouPermutations = []

    for permutation in subTypePermunations:
        temp = []
        for index in permutation:
            if index == 0:
                temp.append([False])
            else:
                temp.append(subSorous[index-1])

        newSubSorouPermutations = [p for p in itertools.product(*temp)]
        for item in newSubSorouPermutations:
            if item not in subSorouPermutations:
                subSorouPermutations.append(item)

    support.tPrint("{} subSorouPermutations Generated".format(len(subSorouPermutations)))

    for i, permutation in enumerate(subSorouPermutations):
        support.tPrint(["New subSorou permutation", '\t', i+1, " of ", len(subSorouPermutations)])
        newSorou = genBaseSorou(aMinVanType)
        for index in range(aMinVanType[0][0]):
            if permutation[index]:
                rotationIndicies = findRotationIndicies(permutation[index], newSorou[index][0])

                for rotationIndex in rotationIndicies:
                    workingSubSorou = support.rotate(permutation[index], rotationIndex[0], rotationIndex[1])
                    if support.isSubSorou(workingSubSorou, newSorou[index][0]):
                        temp = support.subtractSorous(newSorou[index][0], workingSubSorou)
                        newSorou[index].append(temp)

        working = selectFromNewSorou(newSorou)
        support.tPrint(["All potential subsorours generated for this permutation"])
        for item in working:
            if len(item) == support.findWeightofType(aMinVanType):
                sortedItem = support.ssort(item)
                itemKey = support.sorouToLatexString(sortedItem)
                if itemKey not in outputKeys:
                    print(item)
                    outputKeys.add(itemKey)
                    output.append(item)

                # if shouldAdd(item, output):
                #     output.append(item)
                #     print(item)

    if output == []:
        print("ERROR: NO SOROUS GENERATED FOR {}".format(support.typeToLatexString(aMinVanType)))

    return output

def buildTypeSorouArray (aTypeMatchingArray):
    n = len(aTypeMatchingArray)
    toCheck = set()

    Is = itertools.permutations(range(n))
    Js = itertools.permutations(range(n))

    for Iperm in Is:
        for Jperm in Js:
            zipper = zip(Iperm, Jperm)
            toCheck.add(tuple(zipper))

    output = set()

    for item in toCheck:
        shouldAdd = True
        for indexPair in item:
            if not aTypeMatchingArray[indexPair[0]][indexPair[1]]:
                shouldAdd = False

        if shouldAdd:
            output.add(item)

    return output

def makeSubtractedSorouList (aSorou, aSorouList):
    """
    For each sorou in a sorou list, find all ways in which aSorou can be subtracted from it
    Returns the list of all possible results from the subractions
    """
    output = []
    firstRelOrder = support.findRelativeOrder(aSorou)

    for sorou in aSorouList:
        totalRelOrder = support.lcm(firstRelOrder, support.findRelativeOrder(sorou))

        for index in range(totalRelOrder):
            rotated = support.rotate(sorou, totalRelOrder, index)
            temp = support.subtractSorous(rotated, aSorou)

            if len(temp) == len(rotated) - len(aSorou) and shouldAdd(temp, output):
                output.append(temp)

    return output

def selectSorous (aSorouListList, aF0):
    """
    Takes [[f_1, f_2, ...]_j, ...] and returns the list [g_1, g_2, ...] 
    where each g_i is a sum sorous selected one from each list in aSorouListList + aF_0
    """

    output = []

    combinations = [p for p in itertools.product(*aSorouListList)]

    for combo in combinations:
        sorou = []
        for item in combo:
            sorou += item

        sorou += aF0
        output.append(sorou)

    return output

def genAllSorousOfType (aType, aF0, aTypeSorouMap, aPermutationsMap):
    if len(aType) == 1:
        return genAllSorousOfMinVanType(aType, aTypeSorouMap, aPermutationsMap)
    else:
        allSubTypeSorous = []
        for aMinVanType in aType:
            allSubTypeSorous.append(genAllSorousOfMinVanType([aMinVanType], aTypeSorouMap, aPermutationsMap))

        partitions = partitionF0(aF0, len(aType))

        output = []

        for partition in partitions:
            fullMatchingArray = matchTypeToSorous(partition, allSubTypeSorous)
            typeMatchingArray = simplifyArray(fullMatchingArray)

            subSorouToTypesIndicies = buildTypeSorouArray(typeMatchingArray)

            for matchingIndicies in subSorouToTypesIndicies:
                
                unselectedSorouLists = []

                for indexPair in matchingIndicies:
                    G_i = partition[indexPair[0]]
                    S_i = aType[indexPair[1]]
                    S_iSorous = allSubTypeSorous[indexPair[1]] 

                    subtractedSorouList = makeSubtractedSorouList(G_i, S_iSorous)
                    unselectedSorouLists.append(subtractedSorouList)

                sorouList = selectSorous(unselectedSorouLists, aF0)

                for item in sorouList:
                    if shouldAdd(item, output):
                        output.append(support.ssort(item))
                
        return output

def findAllParitiesOfMinVanType (aMinVanType):

    sorous = genAllSorousOfMinVanType(aMinVanType)

    output = set()

    for sorou in sorous:
        par = support.findParity(sorou)
        output.add(par)

    return output

def functionForMap (aType):
    typeStr = support.typeToLatexString(aType)

    if typeStr not in myTypeSorouMap.keys() or myTypeSorouMap[typeStr] == []:
        print(support.findWeightofType(aType), typeStr)

        newSorous = genAllSorousOfMinVanType(aType, myTypeSorouMap, myPermutationMap)
        myTypeSorouMap[typeStr] = newSorous
        print()
        with open('typeMap.txt', 'wb') as tm:
            pickle.dump(myTypeSorouMap, tm)
        return newSorous

    else:
        print()
        return myTypeSorouMap[typeStr]


def main():
    global myTypeSorouMap
    global myPermutationMap

    with open("types.txt", "rb") as tls:
        mTypeList = pickle.load(tls)
    
    with open("typeMap.txt", "rb") as tm:
        myTypeSorouMap = pickle.load(tm)

    myPermutationMap = {}

    #### This is the single threaded implementation
    #### Change the iterable here if only specific type's sorous are needed.
    for i,item in enumerate(mTypeList):
        if support.findWeightofType(item) < 22:
            print(i)
            print(item)
            functionForMap(item)


    #### Multithreading the sorou generation is possible using the following:
    # pool = ThreadPool(4) # Set the number of threads here
    # results = pool.map(functionForMap, typesToMap)

    print('All Sorous Generated')

    return 0

if __name__ == '__main__':
    main()

