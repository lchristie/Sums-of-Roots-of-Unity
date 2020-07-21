import support
import pickle

def filterTL (aTypeList):
    output = []
    for item in aTypeList:
        if support.findWeightofType(item) > 21:
            pass
        else:
            output.append(item)
    return output

def hasEquisigned (aPars):
    for item in aPars:
        if item[0] == item[1]:
            return True
    return False

def writeCSV (aTypeList, aTypeMap, aDocName):
    f = open(aDocName, 'w')
    firstLine = 'Weight,\tTop Prime,\tRelative Order,\tWeight Partition,\tType,\tHeight,\tParities,\tHasEquisigned\n'
    f.write(firstLine)

    for item in aTypeList:
        typeKey = support.typeToLatexString(item)
        try:
            tSorous = aTypeMap[typeKey]
        except:
            print("Error for type: " + typeKey)
            continue
        w = support.findWeightofType(item)

        ts = []

        for aSorou in tSorous:
            if support.isMinVan(aSorou): ts.append(aSorou)

        if ts == []: print("ERROR: Type {} has no minimal vanishing sorou in this typemap".format(typeKey))

        h = set()
        pars = set()
        for aSorou in ts:
            h.add(support.findHeightofSorou(aSorou))
            newPar = support.findParity(aSorou)
            pars.add(newPar)
        wp = support.findWeightPart(item)

        tp = item[0][0]
        relOrder = support.findRelativeOrder(ts[0])
        partition = support.reverseList(wp)

        partitionStr = '('
        for a in partition:
            partitionStr += str(a) + ';'
        partitionStr = partitionStr[:-1]
        partitionStr += ')'

        parityStr = str()
        for aPar in pars:
            # parityStr += str(aPar) + ';'
            parityStr += '(' + str(aPar[0]) + ';' + str(aPar[1]) + ');'
        parityStr = parityStr[:-1]
		
        heightStr = str()
        for height in h:
            heightStr += str(height) + ';'
        heightStr = heightStr[:-1]

        line = '{0},\t{1},\t{2},\t{3},\t{4},\t{5},\t{6},\t{7}\n'.format(w, 
                                                                   tp, 
                                                                   relOrder, 
                                                                   partitionStr, 
                                                                   support.typeToLatexString(item), 
                                                                   heightStr, 
                                                                   parityStr,
                                                                   hasEquisigned(pars))

        f.write(line) 
        print(line)  

    f.close()

def main ():

    with open("types.txt", "rb") as tl:
        typeList = pickle.load(tl)

    with open("typeMap.txt", "rb") as tm:
        typeMap = pickle.load(tm)

    tl = filterTL(typeList)

    print("Writing CSV of type list containing {} types".format(len(tl)))
    writeCSV(tl, typeMap, 'types.csv')
    print("Finsihed writing CSV")

    return 0

if __name__ == '__main__':
    main()
