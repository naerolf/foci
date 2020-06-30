#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-processing script for SR-Tesseler localization and convex hull files in 2D.

This script is a supplemental file for the scientific research article 
"dSTORM microscopy evidences in HeLa cells clustered and scattered gammaH2AX 
nanofoci sensitive to ATM, DNA-PK, and ATR kinase inhibitor"
by Liddle P., Jara-Wilde J., Lafon-Hughes L., Castro I., Härtel S. and Folle G.A.
Currently under revision for publication.

__author__ = "Iván Castro"
__email__ = "ivan.castro@ug.uchile.cl"
__credits__ = "Iván Castro, Jorge Jara-Wilde, Pablo Liddle"
__license__ = "MIT"
__version__ = "0.9-2020.05.17"
__status__ = "Prototype"
"""

import tkinter as tk
from tkinter import filedialog
import csv
import cv2
import numpy as np
from scipy.spatial import ConvexHull

root = tk.Tk()
root.withdraw()


"""Input parameters

distThreshold: Distance threshold to expand the neighborhood of a foci convex hull polygon.
               Set to a value greater that or equal to 0.0

splitChr:      Separator character for the input file. Common characters are
               semicolon ';' , comma ',', or white space ' '.

filetype:      Set 'to txt' or 'csv' as file type.

inputenc:      Specifies the input file encoding. Set by default to 'utf8'.

outputenc:     Specifies the input file encoding. Set by default to 'utf8'.

"""
distThreshold = 0.0

splitChr = ' '

filetype = 'txt'

inputenc = 'utf8'

outputenc = 'utf8'


inputLocationsFull = filedialog.askopenfilename(title = 'Select TXT file with all the locations',
                                                filetypes = (("Text File", "*.txt"),))
if filetype == 'csv':
    inputLocationsROIs = filedialog.askopenfilename(title = 'Select file with the points to define ROIs of exclusion',
                                                    filetypes = (("CSV File", "*.csv"),))
elif filetype == 'txt':
    inputLocationsROIs = filedialog.askopenfilename(title = 'Select file with the points to define ROIs of exclusion',
                                                    filetypes = (("Text File", "*.txt"),))
extPos = inputLocationsROIs.find('.csv')

distThresholdSuffix = '_d' + "".join(str(distThreshold).split()) if distThreshold >= 0 else ''

if filetype == 'csv':
    outputLocations = inputLocationsROIs.split('.csv')[0] + distThresholdSuffix + '.txt'
elif filetype == 'txt':
    outputLocations = inputLocationsROIs.split('.txt')[0] + distThresholdSuffix + '.txt'


def extractFrameNumFromFile(file):
    """Extracts the frame number from a localization filename.
    """
    lineNumber = 0
    with open(file, 'r', encoding = inputenc) as lunID:
        for line in lunID:
            if lineNumber == 0: # first line
                frameNum = line.split('\t')[0]
                break
    return frameNum


def extractRoisFromFile(file):
    """Extracts the frame number from a localization filename.
    """
    roiList = []
    lunID = open(file, 'r', encoding = inputenc)
    while(True):
        header_line = lunID.readline()
        if header_line == '':
            break
        else:
            #print('header_line is ' + header_line)
            roi_id = int(header_line.split(splitChr)[0])
            roi_len = int(header_line.split(splitChr)[1])
            currentRoi = []
            for i in range(roi_len):
                roiPointLine = lunID.readline()
                #print('roiPointLine is ' + roiPointLine)
                point_x = float(roiPointLine.split(splitChr)[0])
                point_y = float(roiPointLine.split(splitChr)[1])
                current_point = (point_x,point_y)
                currentRoi.append(current_point)
            roiList.append(np.asarray(currentRoi))
    return roiList


def generateCoordinatesArray(complementRoiFilteredPointList):
    """Constructs an array of x,y coordinates from the input point list.
    """
    roiCoordinatesArray = []
    for roi in complementRoiFilteredPointList:
        roiCoordinates = []
        for points in roi:
            coordinates = points[:-2]
            roiCoordinates.append(coordinates)
        roiCoordinates = np.asarray(roiCoordinates)
        roiCoordinatesArray.append(roiCoordinates)
    return roiCoordinatesArray


def generateConvexHullArray(roiCoordinatesArray):
    """Computes and returns the convex hull of the input point coordinates array.
    It also returns: the indices and the coordinates of the points belonging to 
    the convex hull, and the convex hull area.
    """
    convexHullArray = []
    survivalIndex = []
    index = 0
    convexHullVerticesArray = []
    holeAreaArray = []
    for roi in roiCoordinatesArray:
        #try:
        if len(roi) > 2:  # checks whether a ROI is correctly "calculated"
            survivalIndex.append(index)
            ch = ConvexHull(roi)
            convexHullArray.append(ch)
            vertices_list = []
            for vertex in ch.vertices:
                vertices_list.append(roi[vertex])
            convexHullVerticesArray.append(vertices_list)
            holeAreaArray.append(ch.volume)
        else:
            print('not enough points to estimate ch:', len(roi))
            print('index problem:', index)
        index += 1
    return convexHullArray, np.asarray(survivalIndex), convexHullVerticesArray, holeAreaArray


def parseScientificNotation(string_s):
    """Taks a numeric string using exponential notation with 'E' to represent 
    a number with <base>E<exponent> format. converts it to a decimal notation 
    (not equivalent), by replacing 'E' with a decimal separator '.'.
    """
    if string_s.find('E') >= 0:
        float_line = float(string_s.split('E')[0].replace(',', '.'))
    else:
        float_line = float(string_s)
    return float_line


def extractPointsFromFile(file):
    """ Reads the content of an input file with tabulator separations 
    """
    pointArr = []
    lineNumber = 0
    fileSep = '\t'
    #currentRoiId = 0 # python 3 unified int and long int
    with open(file, 'r', encoding = inputenc) as lunID:
        for line in lunID:
            if lineNumber == 0: # first line
                frameNum = line.split(fileSep)[0]
            else:
                aux1 = '' + line
                aux2 = aux1.split(fileSep)
                ptX = parseScientificNotation(aux2[0])
                ptY = parseScientificNotation(aux2[1])
                lum = parseScientificNotation(aux2[2])
                ID = parseScientificNotation(aux2[3])
                pointArr.append((ptX, ptY, lum, int(ID)))
            lineNumber += 1
        return pointArr, frameNum


def distanceFromContourFilter(ROI, pList, distance):
    """Using an input region (ROI), filters the points from the input list (pList),
    according to the threshold distance input parameter for proximity. Only the
    points that lie within the proximity threshold are kept.
    """
    filteredPointList = []
    complementFilteredPointList = []
    factor = 1e0
    fRoi = np.rint(ROI*factor).astype(int)
    for point in pList:
        #print(point[0],point[1])
        fPoint = (int(point[0] * factor), int(point[1] * factor))
        ppt = cv2.pointPolygonTest(fRoi, fPoint, True)
        if -ppt >= distance * factor:
            filteredPointList.append(point)
        else:
            complementFilteredPointList.append(point)
    return filteredPointList,complementFilteredPointList 


def diff_list(first, second):
    """Computes the difference between two input sets,
    with the elements in the first set that are not in the second.
    """
    second = set(second)
    return [item for item in first if item not in second]


pointList, frameNum = extractPointsFromFile(inputLocationsFull) 
roiList = extractRoisFromFile(inputLocationsROIs) 
roiFilteredPointList = []
complementRoiFilteredPointList = []
convexHullExpandedList = []

for i in range(len(roiList)):
    print('Processing ROI ', i + 1, ' of ', len(roiList))
    filteredPointList, complementFilteredPointList = distanceFromContourFilter(ROI = roiList[i], pList = pointList, distance = distThreshold)
    roiFilteredPointList.append(filteredPointList)
    complementRoiFilteredPointList.append(complementFilteredPointList)

result = set(complementRoiFilteredPointList[0]) # intersection of points
for s in complementRoiFilteredPointList[1:]:
    result = result.union(s)
result_c = sorted(result, key = lambda tup: tup[3]) # verify

result = diff_list(pointList, result_c)

roiCoordinatesArray = generateCoordinatesArray(complementRoiFilteredPointList)
convexHullArray, survivalIndex, convexHullVerticesArray, holeAreaArray = generateConvexHullArray(roiCoordinatesArray)
totalHoleArea = np.sum(holeAreaArray)

countInPts  = len(result)
countOutPts = len(pointList) - countInPts
with open(outputLocations, 'w', encoding = outputenc) as out:
    csv_out = csv.writer(out, delimiter = '\t', lineterminator = '\n')
    header = (int(frameNum), countInPts)
    csv_out.writerow(header)
    for row in result:
        csv_out.writerow(row)

outputLocations_c = outputLocations.split('.txt')[0] + '_c.txt'

with open(outputLocations_c, 'w', encoding = outputenc) as out:
    csv_out=csv.writer(out, delimiter = '\t', lineterminator = '\n')
    header = (int(frameNum),countOutPts)
    csv_out.writerow(header)
    for row in result_c:
        csv_out.writerow(row)

print('Finished! Number of points in/out/total:', countInPts, countOutPts, countInPts + countOutPts)
print('Output file:', outputLocations)
print('Complement output file:', outputLocations_c)

outputLocations_c_areas = outputLocations.split('.txt')[0] + '_Areas.csv' 

with open(outputLocations_c_areas, 'w', encoding = outputenc) as out:
    csv_out = csv.writer(out, delimiter = ';', lineterminator = '\n')
    header = ('object index', 'area')
    csv_out.writerow(header)
    index = 0
    for area in holeAreaArray:
        csv_out.writerow((index + 1, area))
        index += 1
    csv_out.writerow(('total', totalHoleArea))

print('Area output file:', outputLocations_c_areas)

outputLocations_c_convex_hull_vertices = outputLocations.split('.txt')[0] + '_CHVertices.csv' 

with open(outputLocations_c_convex_hull_vertices, 'w',encoding = outputenc) as out:
    csv_out=csv.writer(out, delimiter = ';', lineterminator = '\n')
    roi_index = 0
    total_roi = len(convexHullVerticesArray)
    header = ('total roi number', total_roi)
    csv_out.writerow(header)
    header = ('roi number', roi_index)
    csv_out.writerow(header)
    for convexHullVertices in convexHullVerticesArray:
        header = ('total vertices', len(convexHullVertices))
        csv_out.writerow(header)
        for vertex in convexHullVertices:
            csv_out.writerow((vertex[0], vertex[1]))
        roi_index += 1
        if roi_index < total_roi:
            header = ('roi number', roi_index)
            csv_out.writerow(header)

print('Convex hull output file:', outputLocations_c_convex_hull_vertices)
