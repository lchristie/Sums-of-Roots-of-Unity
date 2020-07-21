# Sums-of-Roots-of-Unity
This is the code relating to the paper linked. It contains the most recent version of the code as well as examples and outputs. It is split over _ files as explained below. The project it written in Python and depends on the following packages:

Support.py
  This file contains many useful functions for manipulating sorou that are used by the other files. It also contains some tools for formatting the output of such sorou.
  
Types.txt
  This is a pickle dump of a python list containing the currently generated types. 
  
typeMap.txt
  This is a pickle dump of a python dictionary. The keys are the type_strings generated by the support.typeToLatexString() function. The values are lists of all possible sorous of the key's type, up not rotational symmetry. 
  
ReferenceTypes.py
  This file contains a hard coded version of the types of sorou up to weight 16. It is needed in order to begin the type generation process. Running this file resets the "types.txt" file to this hard coded list. 
  
WriteTypeList.py
  Running this file takes the data in the "types.txt" file and the "typeMap.txt" file and writes a csv table containing the following information on each type: weight, top prime, relative order, weight partition, height, and possible parities.

TypeGen.py
  Running this file takes the list of types in "types.txt" and generates all possible types of minimal vanishing sorou of the next weight. It then appends these to the list in "types.txt".

SorouGenerator.py
  Running this file takes the data in the "types.txt" file and the "typeMap.txt" file and generates every possible (up to rotational equivalence) sorou for each type not already in the typeMap.
