import support
import pickle
import pdb


def genNextTypes (aPrevTypes):
    """
    Returns a list of all minimal vanishing types with weight one more than the highest
    weight in aPrevTypes
    """

    output = []

    prevWeight = support.findWeightofType(aPrevTypes[-1])
    currentWeight = prevWeight + 1

    baseTypeWeights = support.makeWeightMap(aPrevTypes)
    allTypeWeights = support.makeFullWeightMap(baseTypeWeights)

    primesToCheck = support.findPrimesToCheck(currentWeight)

    for prime in primesToCheck:
        partitions = support.makepart(currentWeight, prime)

        for partition in partitions:

            f0Possibilities = support.genf0s(partition[-1], prime)

            for f0 in f0Possibilities:

                pTypeList = support.genPTypeList(partition, allTypeWeights, prime, len(f0))
                subTypeListCombinations = support.genSubTypeCombos(pTypeList, aPrevTypes)

                for aSubTypeCombo in subTypeListCombinations:
                    trueSubTypeCombo = support.removeR2(aSubTypeCombo)
                    aNewType = [[prime, f0] + trueSubTypeCombo]
                    testSorou = support.genFromType(aNewType)

                    # if len(f0) > 1: pdb.set_trace()

                    if support.isMinVan(testSorou) and support.findWeightofType(aNewType) == currentWeight:
                        output.append(aNewType)

    return output

def main():

    with open("types.txt", "rb") as tls:
        previousTypeList = pickle.load(tls) 

    print('Starting Now')

    newTypes = genNextTypes(previousTypeList)

    allTypes = previousTypeList + newTypes

    with open("types.txt", "wb") as tls:
        pickle.dump(allTypes, tls) 

    prevWeight = support.findWeightofType(previousTypeList[-1])
    currentWeight = prevWeight + 1

    print('{0} Types Generated of weight {1}'.format(len(newTypes), currentWeight))

    print('--------------')
    for item in newTypes:
        print(support.findWeightofType(item) ,support.typeToLatexString(item))
    print('--------------')

    return 0

if __name__ == '__main__':
    main()







