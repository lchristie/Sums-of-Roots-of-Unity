""" Class Sorou:

    [ [ o_1, p_1 ], ..., [ o_j, p_j ] ]
    \nu_{o_1}^{p_1} + ... + \nu_{o_j}^{p_j}
"""

""" Class CommonRoot:

    [ o, p_1, ... p_j ]
    \nu_o^{p_1} + ... + \nu_o^{p_j}
"""

""" Class SubSorous:

    [ prime, [ p_1, f_1 ], ..., [ p_j, f_j ] ]
    \sum_{i = 1}^{j} \nu_prime^{p_i} \times f_i
"""

""" Class Type:

    [[ p, f_0, T_1, ... T_j ], ..., [ q, g_0, S_1, ... S_j ]]
    ( R_p : f_0 : T_1, ..., T_j ) (+) ( R_q : g_0 : S_1, ... S_j )
"""

import math
import cmath
import fractions
import pickle
import pdb
import itertools
import time

from sympy import minimal_polynomial, sqrt, solve, QQ, I, simplify, expand
from sympy.functions import re, im
from sympy.abc import x, y

vanishingLimit = (0.1) ** (6)

q = { 1: [[1]] } # Used for partition generation

def decompose(n):
    try:
        return q[n]
    except:
        pass

    result = [[n]]

    for i in range(1, n):
        a = n-i
        R = decompose(i)
        for r in R:
            if r[0] <= a:
                result.append([a] + r)

    q[n] = result
    return result
    
def filterToLength (partitions, l):
    output = []
    for item in partitions:
        if len(item) == l:
            output += [item]
    return output

def filterToTwo (aPartition):
    output = []
    for item in aPartition:
        if 1 not in item and 2 in item:
            output += [item]
    return output

def filterP (aPartition, l):
    output = filterToLength(aPartition, l)
    return output

def makepart (n, l):
    output = filterP(decompose(n),l)
    return output

def reverseList (aList):
    output = []
    for index in range(len(aList)):
        output.append(aList[-index - 1])
    return output

def isPrime(number):
    if number > 1:
        if number == 2:
            return True
        if number % 2 == 0:
            return False
        for current in range(3, int(math.sqrt(number) + 1), 2):
            if number % current == 0: 
                return False
        return True
    return False

def getPrimes(number):
    while True:
        if isPrime(number):
            yield number
        number += 1

def findNthPrime(n):
    if n == 0:
        return 1
    
    primeGenerator = getPrimes(1)

    for index in range(n - 1):
        next(primeGenerator)

    return next(primeGenerator)

def findNtoMPrimes(n, m):

    output = []

    primeGenerator = getPrimes(1)

    for index in range(n - 1):
        next(primeGenerator)

    for index in range(m - n):
        output.append(next(primeGenerator))

    return output

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    value = a * b // gcd(a, b)
    return value

def lcmm(args):
    tempLcm = 1
    index = 0
    while index < len(args):
        tempLcm = lcm(tempLcm, args[index])
        index += 1
    return tempLcm

def findMaxPrime (aWeightLimit):
    index = 0
    while findNthPrime(index + 1) < aWeightLimit:
        index += 1
    return findNthPrime(index)

def findMaxPrimeIndex (aWeightLimit):
    ps = findMaxPrime(aWeightLimit)
    index = 0
    while findNthPrime(index) != ps:
        index += 1
    if findNthPrime(index) == ps:
        return index
    else:
        return False

def pFactor (n):
    pFactors = []

    prime = 1
    primes = getPrimes(1)

    while prime <= n:
        if n % prime == 0:
            pFactors.append(prime)
        prime = next(primes)

    return pFactors

def containsEquiPar (aParList):
    for par in aParList:
        if par[0] == par[1]:
            return True
    return False

def toRoot (z):
    if abs(z - 1) > vanishingLimit:
        print("This number is not a root of Unity")
        return False
    else:
        phi = cmath.phase(z) / (2 * cmath.pi)
        ratio = fractions.from_float(phi)
        return simplifyRoot(ratio.numerator, ratio.denominator)

def fromRoot (aRoot):
    z = cmath.exp(2*cmath.pi*aRoot[1]*1j / aRoot[0])
    return z

def simplifyRoot (aOmega, aPower):
    tempOmega = aOmega
    tempPower = aPower % aOmega

    newOmega = tempOmega // gcd(tempOmega, tempPower)
    newPower = tempPower // gcd(tempOmega, tempPower)

    return [newOmega, newPower]

def multiplyRoots (aOmega1, aPower1, aOmega2, aPower2):
    tempOmega = aOmega1 * aOmega2
    tempPower = aPower1 * aOmega2 + aOmega1 * aPower2
    return simplifyRoot(tempOmega, tempPower)

def findRootOrder (aRoot):
    flag = False
    index = 1
    while not flag:
        if index * aRoot[1] % aRoot[0] == 0:
            flag = True
            return index
        index += 1
        if index > 1000000000:
            print("Order too high")
            return False
    print("Encountered some error")
    return False

def findOrder (aSorou):
    orders = []
    for item in aSorou:
        orders.append(findRootOrder(item))
    temp = set(orders)
    uniqueOrders = []
    for item in temp:
        uniqueOrders.append(item)
    return lcmm(uniqueOrders)

def findRelativeOrder (aSorou):
    rootRatios = []
    for item in aSorou:
        rootRatios.append(item)
        for other in aSorou:
            newRoot = [item[0] * other[0], item[0] * other[1] - item[1] * other[0]]
            rootRatios.append(newRoot)
    return findOrder(rootRatios)

def rotate (aSorou, aOmega, aPower):
    output = []
    for item in aSorou:
        output.append(multiplyRoots(item[0],item[1],aOmega,aPower))
    return output

def rotateFor1 (aSorou):
    # Rotates a sorou so that it includes a 1
    firstRoot = aSorou[0]
    return rotate(aSorou, firstRoot[0], (firstRoot[0] - firstRoot[1]) % firstRoot[0])

def orderRoots (root1, root2):
    if root1[0] > root2[0]:
        return 1
    elif root2[0] > root1[0]:
        return 2
    elif root1[1] > root2[1]:
        return 1
    elif root2[1] > root1[1]:
        return 2
    else:
        return 0

def ssort (aSorou):
    tempSorou = []
    for root in aSorou:
        tempSorou.append(simplifyRoot(root[0],root[1]))
    return sorted(tempSorou, key=lambda tup: (tup[0],tup[1]) )

def isEqual (Sorou1, Sorou2):
    if len(Sorou1) != len(Sorou2):
        return False

    sSorou1 = ssort(Sorou1)
    sSorou2 = ssort(Sorou2)

    for i in range(len(sSorou1)):
        root1 = simplifyRoot(sSorou1[i][0], sSorou1[i][1])
        root2 = simplifyRoot(sSorou2[i][0], sSorou2[i][1])
        if root1[0] != root2[0]:
            return False
        if root1[1] != root2[1]:
            return False

    return True

def isEquiv (Sorou1, Sorou2):
    relOrder1 = findRelativeOrder(Sorou1)
    relOrder2 = findRelativeOrder(Sorou2)
    relOrder = lcm(relOrder1, relOrder2)

    tempSorou = Sorou2
    for i in range(relOrder):
        if isEqual(Sorou1, tempSorou):
            return True
        else:
            tempSorou = rotate(tempSorou, relOrder, 1)

    return False

def genRP (aPrime):
    RP = []
    for i in range(aPrime):
        RP.append([aPrime,i])
    return RP

def genAllSubSorous (aSorou):
    output = []
    for root in aSorou:
        newItems = []
        for item in output:
            newItems.append(item + [root])
        output.append([root])
        if newItems: output += newItems
    return output

def findWeightofSorou (aSorou):
    weight = len(aSorou)
    return weight

def findWeightofType (aTypeListSum):
    weight = 0
    for typelist in aTypeListSum:
        weightList = [typelist[0], findWeightofSorou(typelist[1])]
        # print("weightList: ", weightList)
        for subType in typelist[2:]:
            if len(subType) == 1 & len(subType[0]) == 2:
                weightList.append(subType[0][0])
            else:
                weightList.append(findWeightofType(subType))
        # print("weightList with SubTypes: ", weightList)
        weightPartition = []
        for item in weightList[2:]:
            weightPartition.append(item - weightList[1])
        while len(weightPartition) < typelist[0]:
            weightPartition.append(weightList[1])
        # print("weightPartition: ", weightPartition)
        weight += sum(weightPartition)
    return weight

def findWeightPart (aMinVanType):
    '''
    This takes a minimal vanishing type: [[p, f0, T1, T2, ...]]
    and returns it's weight partition : [a,b,c,...]
    '''
    output = []
    wf0 = len(aMinVanType[0][1])
    p = aMinVanType[0][0]

    for item in aMinVanType[0][2:]:
        w = findWeightofType(item)
        output.append(w - wf0)

    for index in range(p - len(aMinVanType[0][2:])):
        output.append(wf0)

    return output

def findHeightofSorou (aSorou):
    if len(aSorou) == 0:
        return 0

    working = ssort(aSorou)

    listsOfHeights = [1]
    
    for index in range(1,len(working)):
        if working[index] == working[index-1]:
            listsOfHeights[-1] += 1
        else:
            listsOfHeights.append(1)

    return max(listsOfHeights)

def toCommonRoot (aSorou):
    """
    Takes a sorou [[o,p],...,[o,p]] and converts it to sum of powers of a common root

    Stored as a Common Root list [relorder, l_1, l_2, ..., l_n] = sum_{i=1}^n \nu_{relorder}^{l_i}
    """
    relOrder = findRelativeOrder(aSorou)
    output = [relOrder]
    for item in aSorou:
        factor = relOrder // item[0]
        output.append(item[1]*factor)

    return output

def fromCommonRoot (aCommonRoot):
    """
    Takes a Common Root list [relorder, l_1, l_2, ..., l_n] and returns a Sorou List [[o,p],...,[o,p]]
    """
    sorou = []
    relOrder = aCommonRoot[0]
    powers = aCommonRoot[1:]
    for power in powers:
        sorou.append([relOrder,power])
    return ssort(sorou)

def splitRoot (aRoot, aTopPrime, aOtherPrimes):
    """
    This takes a root w = v_relorder^a: [relorder, a]
    And returns two roots with product w, one a power of v_topPrime and the other a power of v_otherPrimes
    Returns FALSE if error
    """
    temp = simplifyRoot(aRoot[0], aRoot[1])

    if temp[0] == aTopPrime:
        return (temp, [1,0])

    root1 = [ aTopPrime, 0 ]
    root2 = [ aOtherPrimes, 0 ]

    for a in range(aTopPrime):
        for b in range(aOtherPrimes):
            test = multiplyRoots(root1[0],root1[1],root2[0],root2[1])
            if test[0] == temp[0] and test[1] == temp[1]:
                return (simplifyRoot(root1[0],root1[1]),simplifyRoot(root2[0],root2[1]))
            else:
                root2[1] += 1
        root1[1] += 1
    print("ERROR: Cannot split root")
    return False

def formatted (aCommonRoot):
    """
    Helper function for toSubSorous
    Takes a list of common roots and returns it in subsorou format.
    Only use for prime order common roots (returns False otherwise)
    """
    output = []

    if aCommonRoot[0] not in primes:
        return False

    rootMap = {x : [] for x in range(aCommonRoot[0])}

    for power in aCommonRoot[1:]:
        rootMap[power].append([1,0])

    for item in rootMap.keys():
        output.append([item, rootMap[item]])

    return output

def rotateSSorousW (ssorou):
    """
    Rotates the subsorou form [p, [o,f], ..., [o,f]] so that

        w(f_0) < w(f_j) for all j

    """
    output = [ssorou[0]]
    working = ssorou[1:]

    minIndex = 0
    minWeight = findWeightofSorou(working[0][1])
    for index in range(len(working)):
        newWeight = findWeightofSorou(working[index][1])
        if (newWeight < minWeight and newWeight > 0) or minWeight == 0:
            minIndex = index
            minWeight = newWeight

    rotateOrder = ssorou[0] - working[minIndex][0]

    for index in range(len(working)):
        originalIndex = (index + rotateOrder) % ssorou[0]
        newSubSorou = [index, working[(index + rotateOrder) % ssorou[0]][1]]
        output.append(newSubSorou)

    return output

def rotateSSorous1 (ssorou):
    """
    rotates the subsorou form [p, [o,f], ..., [o,f]] so that

        1 -< f_0

    """
    output = [ssorou[0]]
    working = ssorou[1:]

    firstRoot = ssorou[1][1][0]

    rotateOmega = firstRoot[0]
    rotatePower = (firstRoot[0] - firstRoot[1]) % rotateOmega

    for subSorou in working:
        o = subSorou[0]
        newSubSorou = rotate(subSorou[1], rotateOmega, rotatePower)
        output.append([o, newSubSorou])

    return output

def rotateSubSorous (ssorou):
    temp = rotateSSorousW(ssorou)
    return rotateSSorous1(temp)

def toSubSorous (aSorou):
    """
    Takes a sorou [[o,p],...,[o,p]] and converts it to the sum of products of roots of order p_s a sub sorous

    Stored as [p, [o,f], ..., [o,f]]
    """
    output = []
    working = toCommonRoot(aSorou)

    relOrder = working[0]
    powers = working[1:]

    thesePrimes = pFactor(relOrder)

    if len(thesePrimes) == 1:
        output.append(thesePrimes[0])
        temp = formatted(working)
        if temp: output += temp
        return output

    topPrime = thesePrimes[-1]

    output.append(topPrime)

    rootMap = {x: [] for x in range(topPrime)}

    rotators = []
    splits = []
    SubSorous = []

    # Splits out top prime roots
    for power in powers:
        root = simplifyRoot(relOrder, power)
        first, second = splitRoot(root, topPrime, relOrder // topPrime)
        rotators.append(first)
        splits.append(second)
    for index in range(len(rotators)):
        root = rotators[index]
        rootMap[root[1]] += [splits[index]]
    for index in range(topPrime):
        output.append([index, rootMap[index]])

    temp = rotateSubSorous(output)

    output2 = [output[0]]
    for item in temp[1:]:
        if item[1]: output2.append(item)

    return output2

def printSubSorous (aSubSorou):
    print(aSubSorou[0])
    for item in aSubSorou[1:]:
        print(item)

def fromSubSorous (ssorou):
    topPrime = ssorou[0]
    output = []

    for item in ssorou[1:]:
        for item2 in item[1]:
            newroot = multiplyRoots( topPrime, item[0], item2[0], item2[1] )
            output.append(newroot)

    return output

def evalSorouN (aSorou):
    """
    Takes a Sorou list [[o,p],...,[o,p]] and returns the value as
    a complex number z
    """
    nums = []
    for root in aSorou:
        nums.append(cmath.exp(2*cmath.pi*root[1]*1j / root[0]))
    return sum(nums)

def isVanishingN (Sorou):
    value = evalSorouN(Sorou)
    if abs(value) < vanishingLimit:
        return True
    else:
        return False

def findPrimitiveRoot (n):
    if (n == 1):
        return solve(x - 1)[0]
    else:        
        unorderedRoots = solve(x**n - 1)
        if unorderedRoots[0] != 1:
            currentMax = unorderedRoots[0]
        else:
            currentMax = unorderedRoots[1]
        for item in unorderedRoots:
            if re(item) >= re(currentMax) and im(item) > 0:
                currentMax = item
        return currentMax

def raiseAlgebraicToPower(aPrimitiveRoot, aPower):
    output = 1

    for index in range(aPower):
        output = simplify(expand(output * aPrimitiveRoot))

    return output

def isVanishing (aSorou):
    algebraicRoots = []

    for root in aSorou:
        primitiveRoot = findPrimitiveRoot(root[0])
        algebraicRoots.append(simplify(expand(primitiveRoot ** root[1])))

    return sum(algebraicRoots) == 0

def isMinVan (aSorou):
    """
    Based on Prop 2.3, check if the sorou is minimal vanishing.
    """
    if not aSorou:
        return False

    working = toSubSorous(aSorou)
    SubSorous = working[1:]

    if not isVanishingN(aSorou):
        return False

    # First criteria
    if isVanishingN(working[1][1]):
        return False

    # Second & Third Criteria
    values = []
    for item in SubSorous:
        valsToAdd = set()
        allSubSorous = genAllSubSorous(item[1])[:-1]
        for subitem in allSubSorous:
            valsToAdd.add(evalSorouN(subitem))
            if isVanishingN(subitem):
                if len(subitem) != len(item[1]):
                    return False
        values.append(valsToAdd)

    # Third Criteria
    if len(values[0]) == 0:
        return True

    intersection = values[0]
    for vals in values[1:]:
        intersection = intersection & vals

    if len(intersection) == 0:
        return True

    return False

def orderTypes (typeSum1, typeSum2):
    """
    Compares the types type1 and type2 and returns either 1 or 2 depending on lower order type.
    If the types are the same returns 0

        typeSum1 = [[[p, f0] + [t_1, ..., t_j]], ..., [...]]
        typeSum2 = [[[q, g0] + [s_1, ..., s_i]], ..., [...]]

        (with t_a and s_a listed in ascending order)

    Order is defined by type1 > type2 if and only if (in order or presidence):

        I.  Weight Type1 > Weight Type2
        I.  Number of minimal vanishing sorous of type1 > type2
        II. For minvan types listed in decending order:
            1.  p >= q
            2.  w(f0) > w(g0)
            3.  For roots e^{\phi_a i} in f0, e^{\theta_a i} in g0:
                A. \phi_a > \theta_a
            4.  j > i
            5.  For subTypes t_a, s_a:
                A.  t_a > s_a

    """
    if findWeightofType(typeSum1) > findWeightofType(typeSum2):
        return 1
    elif findWeightofType(typeSum2) >  findWeightofType(typeSum1):
        return 2

    if len(typeSum1) > len(typeSum2):
        return 1
    elif len(typeSum2) > len(typeSum1):
        return 2

    for index in range(len(typeSum1)):
        type1 = typeSum1[index]
        type2 = typeSum2[index]

        if type1 == type2:
            continue

        if type1[0] > type2[0]:
            return 1
        elif type2[0] > type1[0]:
            return 2

        f0 = type1[1] 
        g0 = type2[1]
        wf0 = findWeightofSorou(f0)
        wg0 = findWeightofSorou(g0)

        if wf0 > wg0:
            return 1
        elif wg0 > wf0:
            return 2

        for sorouIndex in range(wf0):
            phi1 = f0[sorouIndex][1] / f0[sorouIndex][0]
            phi2 = g0[sorouIndex][1] / g0[sorouIndex][0]
            if phi1 > phi2:
                return 1
            elif phi2 > phi1:
                return 2

        if len(type1) == 2 and len(type2) == 2:
            continue
        elif len(type1) > len(type2):
            return 1
        elif len(type2) > len(type1):
            return 2

        for index in range(2,len(type1)):
            subTypeOrder = orderTypes(type1[index], type2[index])
            if subTypeOrder == 1:
                return 1
            elif subTypeOrder == 2:
                return 2

    return 0

def subtractSorous (Sorou1, Sorou2):
    '''
    returns sorou1 - sorou2
    '''
    output = []
    sumand = ssort(Sorou1)
    negand = ssort(Sorou2)

    index1 = 0
    index2 = 0

    tempSumand = []
    tempNegand = []

    while index1 < len(sumand) and index2 < len(negand):
        root1 = sumand[index1]
        root2 = negand[index2]

        tester = orderRoots(root1, root2)
        if tester == 0:
            index1 += 1
            index2 += 1
        elif tester == 1:
            tempNegand.append(root2)
            index2 += 1
        else:
            tempSumand.append(root1)
            index1 += 1

    if index1 < len(sumand):
        tempSumand += sumand[index1:]
    if index2 < len(negand):
        tempNegand += negand[index2:]

    output += tempSumand
    output += rotate(tempNegand, 2, 1)

    return output

def subtractF0 (f0, aSubSorou):
    working = aSubSorou.copy()

    for index in range(len(aSubSorou)):
        tester = subtractSorous(f0, working)
        if len(tester) == len(aSubSorou) - len(f0):
            return tester

        temp = working[1:] + [working[0]]
        rotateOmega = temp[0][0]
        rotatePower = (rotateOmega - temp[0][1]) % rotateOmega
        working = rotate(temp, rotateOmega, rotatePower)

    return False

def isSubSorou (aSorou1, aSorou2):
    # Checks if aSorou2 is contained within aSorou1
    working = subtractSorous(aSorou1, aSorou2)
    if len(working) == len(aSorou1) - len(aSorou2):
        return True
    else:
        return False

def makeWeightMap (aTypeList):
    theTypes = {}
    for T in aTypeList:
        w = findWeightofType(T)
        if w in theTypes.keys():
            theTypes[w] += [T]
        else:
            theTypes[w] = [T]
    maxWeight = sorted(theTypes.keys())[-1]
    for w in range(maxWeight):
        if w not in theTypes.keys():
            theTypes[w] = []
    return theTypes

def makeFullWeightMap (aWeightMap):
    weightLimit = sorted(aWeightMap.keys())[-1]

    output = {0:[], 1:[]}

    for index in range(2,weightLimit):
        newTypes = []
        # indexPartitions = filterToTwo(decompose(index))
        indexPartitions = decompose(index)
        for partition in indexPartitions:
            weightList = [aWeightMap[innerIndex] for innerIndex in partition]
            combos = [p for p in itertools.product(*weightList)]
            for item in combos:
                newType = []
                for aType in item:
                    newType += aType
                if checkTypesIsOrdered(item) and newType not in newTypes:
                    newTypes.append(newType)
        output[index] = newTypes

    return output

def orderByWeight (aTypeList):
    output = []
    theTypes = makeWeightMap(aTypeList)
    for index in sorted(theTypes.keys()):
        output += theTypes[index]
    return output

def checkTypesIsOrdered (aTypeList):
    for index in range(1,len(aTypeList)):
        if orderTypes(aTypeList[index-1], aTypeList[index]) == 2:
            return False
    return True

def genFromTypeForF0 (SorouTypeSum, f0):
    """
    Generates a sorou for each minimial vanishing type and rotates them so that it contains f_0
    """
    sorouList = []
    for innerType in SorouTypeSum:
        minVanType = [innerType]
        temp = genFromType(minVanType)
        sorouList.append(rotateFor1(temp))

    if len(sorouList) < len(f0):
        return False

    output = []

    for index in range(len(f0)):
        rotateOmega = f0[index][0]
        rotatePower = f0[index][1]
        output += rotate(sorouList[index],rotateOmega, rotatePower)

    for subSorou in sorouList[len(f0):]:
        output += subSorou

    return output

def genFromType (SorouTypeSum):
    """
    Takes a type [[ p, f_0, T_1, T_2, .., T_j ], ... ]
    and returns a representative sorou with this type
    """
    sorou = []

    for SorouType in SorouTypeSum:
        if len(SorouType) == 2:
            sorou += genRP(SorouType[0])
        else:
            subTypes = SorouType[2:]
            if len(subTypes) >= SorouType[0]: return False
            subSorous = [SorouType[1]]
            differences = []
            for subType in subTypes:
                if len(subType) == 1:
                    temp = genFromType(subType)
                else:
                    temp = genFromTypeForF0(subType, SorouType[1])
                if not temp: return False
                differences.append(ssort(temp))

            for item in differences:
                subSorou = subtractF0(SorouType[1], item)
                if not subSorou: return False
                subSorous.append(subSorou)

            tempLen = SorouType[0] - len(subTypes) - 1
            if tempLen:
                for index in range(tempLen):
                    subSorous.append(SorouType[1])

            sSorou = [SorouType[0]]
            for index in range(SorouType[0]):
                sSorou.append([index, subSorous[index]])

            sorou += fromSubSorous(sSorou)

    return ssort(sorou)

def findParity (aSorou):
    """
    Note that a root (-1)^n \nu_o^p in minimal terms can be rewritten +\nu_{o'}^{p}
    With o' % 2 == 0 if and only if n % 2 = 1
    """
    odds = 0
    evens = 0
    for root in aSorou:
        newRoot = simplifyRoot(root[0], root[1])
        if newRoot[0] % 2 == 0:
            evens += 1
        else:
            odds += 1
    return (max(odds, evens), min(odds,evens))

def findTopPrimeIndex (aWeightLimit):
    index = 0 
    while findNthPrime(index) < aWeightLimit:
        index += 1
    return index

def genP1MinVanTypeList (aWeightLimit):
    """
    generates the first types based on [primes] defined above
    All the ones based on R_p and then all the ones based on R_{p+1}
    returns a list of the types
    """

    n = findTopPrimeIndex(aWeightLimit)

    primeTypes = [ [[p, [[1,0]]]] for p in findNtoMPrimes(2,n) ]
    allTypes = [primeTypes[0]]

    for item in primeTypes[1:]:
        topPrime = item[0][0]
        toAdd = [[item]]

        for index in range(topPrime-1):
            innerToAdd = []
            previous = toAdd[index]

            for prevType in previous:
                for extraType in allTypes:
                    if len(prevType[0]) == 2:
                        newType = [prevType[0] + [extraType]]
                        w = findWeightofType(newType)
                        if w < aWeightLimit:
                            innerToAdd.append(newType)
                    elif orderTypes(extraType, prevType[0][-1]) in [0,2]:
                        newType = [prevType[0] + [extraType]]
                        w = findWeightofType(newType)
                        if w < aWeightLimit:
                            innerToAdd.append(newType)
            toAdd.append(innerToAdd)

        for innerList in toAdd:
            for subInnerList in innerList:
                allTypes.append(subInnerList)

    return allTypes

def genVanTypes (aWeightLimit):
    output = []
    minVanTypes = genP1MinVanTypeList(aWeightLimit)

    vanillaTypeList = reverseList([[[2,[[1,0]]]]] + minVanTypes)

    for i in range(len(vanillaTypeList)):
        firstType = vanillaTypeList[i]
        toAdd = [[firstType]]
        for index in range(aWeightLimit - findWeightofType(firstType)):
            innerToAdd = []
            try:
                previous = toAdd[index]
            except:
                continue

            for prevType in previous:
                for extraType in vanillaTypeList[i:]:
                    if findWeightofType(prevType + extraType) < aWeightLimit:
                        innerToAdd.append(prevType + extraType)
            if len(innerToAdd) > 0: toAdd.append(innerToAdd)
        for inner in toAdd:
            output += inner

    return output

def checkForAMinVanType (aTypeList, minVanTypes):
    for aType in aTypeList:
        if aType in minVanTypes:
            return True
    return False

def reduceByPs (someTypes, aTopPrime):
    output = []
    for item in someTypes:
        if item[0][0] < aTopPrime:
            output.append(item)
    return output    

def genPTypeList (aPartition, aTypeMap, aTopPrime, aWF0):
    output = []
    for index in aPartition:
        someTypes = aTypeMap[index+aWF0]
        filtered = reduceByPs(someTypes, aTopPrime)
        output.append(filtered)
    return output

def genSubTypeCombos (aPTypeList, minVanTypes):
    output = []

    endIndex = len(aPTypeList)          # Note that this is the top prime
    listLens = []
    for item in aPTypeList:
        listLens.append(len(item))

    indexLists = []

    for i in range(endIndex):
        newIndexLists = []
        if i == 0:
            for index in range(listLens[i]):
                newIndexLists.append([index])
        else:    
            for item in indexLists:
                for index in range(listLens[i]):
                    newIndexLists.append( item + [index] )

        indexLists = newIndexLists

    for l in indexLists:
        typeCombo = []
        for i in range(endIndex):
            typeCombo.append(aPTypeList[i][l[i]])
        if checkForAMinVanType(typeCombo, minVanTypes) and checkTypesIsOrdered(typeCombo):
            output.append(typeCombo)

    return list(output)

def isR2Sum (aType):
    if aType == [[2,[[1,0]]]]:
        return True
    if aType == [[2,[[1,0]]],[2,[[1,0]]]]:
        return True
    if aType == [[2,[[1,0]]],[2,[[1,0]]],[2,[[1,0]]]]:
        return True
    return False

def removeR2 (aSubTypeCombo):
    output = []
    for item in aSubTypeCombo:
        if not isR2Sum(item):
            output.append(item)
    return output

def prodPrevPrimes (aPrime):
    index = 0
    output = 1
    while findNthPrime(index) < aPrime:
        output *= findNthPrime(index)
        index += 1
    return output

def genf0s (aWeight, aTopPrime):
    layerOutput = [[[1,0]]]
    relOrder = prodPrevPrimes(aTopPrime)
    for index in range(aWeight - 1):
        nextLayerOutput = []
        for item in layerOutput:
            if index == 0:
                for i in range(1,relOrder // 2):
                    nextLayerOutput.append(item + [[relOrder, i]])
            else:
                for i in range(1, item[-1][-1]):
                    nextLayerOutput.append(item + [[relOrder, i]])
        layerOutput = nextLayerOutput
    return layerOutput

def findPrimesToCheck(aWeight):
    # Note that we only consider primes stricly greater than 5
    # Because the maximum weight type with top prime 5 is
    # (R_5 : 4R_3) which is already hard coded into the reference types   
    output = []
    primeGenerator = getPrimes(6)

    p = next(primeGenerator)
    while p <= aWeight:
        output.append(p)
        p = next(primeGenerator)

    return output

######## Functions for formatting output:

def printNestedList (aList, counter = 0):
    for item in aList:
        if type(item) != list:
            print(' ' * counter, item)
        elif type(item[0]) != list:
            print(' ' * counter, item)
        elif type(item[0][0]) != list:
            print(' ' * counter, item)
        else:
            print(' ' * counter, '[')
            printNestedList(item, counter + 1)
            print(' ' * counter, ']')

def rootToLatexString (root):
    """
    Takes a root [a,b] and returns a string "\nu_{a}^{b}"
    """
    temp = simplifyRoot(root[0],root[1])
    if temp == [1,0]:
        return str(1)
    else:
        return "\\nu_{" + str(temp[0]) + "}^{" + str(temp[1]) + "}"

def sorouToLatexString (Sorou):
    """
    Takes a sorou: 
    and returns a latex formatted string "\nu_o^p + ..."
    """
    output = ""
    for index in range(len(Sorou)):
        output += rootToLatexString(Sorou[index])
        if index != len(Sorou) - 1:
            output += "+"
    return output

def sorouToTikzString (aSorou):
    print('\\begin\{tikzpicture\}')
    print('\\draw[<->, gray] (-4,0) -- (4,0);')
    print('\\draw[<->, gray] (0,-4) -- (0,4);')
    print('\\draw[dashed, thin, gray] (0,0) circle (3);')

    for aRoot in aSorou:
        phi = (360 / aRoot[0]) * aRoot[1]
        psi = round(phi)
        colour = 'black'
        print('\\draw[fill = {}] (0,0) -- ({}:3) circle (0.7mm);'.format(colour, psi))

    print('\\end\{tikzpicture\}')

def npMinVanTypeToLatexString (typelist):
    """
    Takes a type list [ p, f_0 ] + [T_1, T_2, .., T_p ]
    and creates a string "(R_p : f_0 : T_1, ... T_j)"
    """
    output = "R_"

    output += '{' + str(typelist[0]) + '}'
    f0str = sorouToLatexString(typelist[1])
    if f0str != "1": output += " : " + f0str

    if len(typelist) == 2:
        return output
    else:
        output += " : "

    index = 2
    counter = 1
    while index < len(typelist):
        curType = typelist[index]
        nextType = typelist[index + 1] if index + 1 < len(typelist) else 0

        if curType == nextType:
            counter += 1
            index += 1
        elif curType[0] != [2,[[1,0]]]:
            if counter > 1:
                output += str(counter) + typeToLatexString(curType)
            else:
                output += typeToLatexString(curType)

            if index != len(typelist) - 1:
                output += "; " # comma for latex, semi-colon for csvs            
            index += 1
            counter = 1
        else:
            if index == len(typelist) - 1:
                output = output[:-2]  
            index += 1
            counter = 1

    return output

def minVanTypeToLatexString (typelist):
    # Add brackets to npTypeString
    if len(typelist) == 2:
        return npMinVanTypeToLatexString(typelist)
    else:
        return "( " + npMinVanTypeToLatexString(typelist) + " )"

def typeToLatexString (typeListSum):
    output = ""
    if len(typeListSum) > 1:
        output = "("

    for innerType in typeListSum[:-1]:
        output += minVanTypeToLatexString(innerType) + ' \oplus '
    output += minVanTypeToLatexString(typeListSum[-1])

    if len(typeListSum) > 1:
        output += ")"

    return output

def tPrint (anArgs):

    output = time.asctime() + "; "

    for item in anArgs:
        output = output + str(item)

    print(output)

######## Build a full type list:

def main():
    
    return 0

if __name__ == '__main__':
    main()







