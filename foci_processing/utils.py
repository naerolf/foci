def extractPointsFromFile(file):
    pointArr = []
    lineNumber = 0
    #currentRoiId = 0 # python 3 unified int and long int
    with open(file, 'r', encoding='utf8') as lunID:
        for line in lunID:
            if lineNumber == 0: # first line
                frameNum = line.split('\t')[0]
            else:
                aux1 = ''+line
                aux2 = aux1.split('\t')
                ptX = parseScientificNotation(aux2[0])
                ptY = parseScientificNotation(aux2[1])
                lum = parseScientificNotation(aux2[2])
                ID = parseScientificNotation(aux2[3])

                pointArr.append((ptX,ptY,lum,int(ID)))
            lineNumber += 1
        return pointArr,frameNum
    
def parseScientificNotation(string_s):
    if string_s.find('E') >= 0:
        float_line = float(string_s.split('E')[0].replace(',','.'))
    else:
        float_line = float(string_s)
    return float_line