# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:04:21 2018
This function cleans the MLC sequence of a VMAT sequence

@author: byron wilson, byron.wilson@bccancer.bc.ca
"""
import pydicom as dcm
import numpy as np

def loadContours(RS):
    contours = []
    for i, structure in enumerate(RS.StructureSetROISequence):
        if 'PTV' in structure.ROIName:
            print(structure.ROIName)
#        if 'PTV' in structure.ROIName:
            for contour in RS.ROIContourSequence[i].ContourSequence:
                contours.extend(contour.ContourData)
    contours = np.array(contours).reshape([len(contours)//3,3])
    
    return contours

def rotationMatrix(theta,direction):
    c, s = np.cos(theta), np.sin(theta)
    if direction == 'x':
        R = [[1,0,0],
             [0,c,-s],
             [0,s,c]]
    elif direction == 'y':
        R = [[c,0,s],
             [0,1,0],
             [-s,0,c]]
    else: # direction assumed to be 'z'
        R = [[c,-s,0], 
             [s,c,0],
             [0,0,1]] 
    return np.array(R)

def rotatePoints(points,direction,angle):
    R = rotationMatrix(angle,direction)
    return np.array([R.dot(point) for point in points])


def transformPoints(contours, iso, gantryAngle, collimatorAngle, couchAngle):
    '''Transform the points to the beam axis coordinate system 
    '''
    ##Subtract the isocentre
    contours = contours - np.matmul(iso,np.ones((1,contours.shape[0]))).T
    
    ##Rotate phantom along y axis (ventral) by couchAngle
    contours = rotatePoints(contours,'y',couchAngle)
    ##Rotate by Gantry Angle along z axis (superior) so that the long axis of
    ##the beam is aligned with the y axis
    contours = rotatePoints(contours, 'z',gantryAngle)
    ##Rotate the phantom along the y (which is now aligned with long axis of beam)
    contours = rotatePoints(contours,'y', collimatorAngle)
    
    ##Correct for the divergence of the beam
    contours = correctDivergence(contours,1000)
    
    
    return contours

def correctDivergence(contours, SAD):
    '''Fix for the divergence of the beam'''
    #The MLC 
    scalingFactor = SAD/(SAD+(contours[:,1]))
    #multiply x and z by the scaling factor. 
    contours[:,0] = contours[:,0]*scalingFactor
    contours[:,2] = contours[:,2]*scalingFactor
    return contours


def calculateMLCEdges():
    ''''''
    MLCEdge = [-110]
    iEdge = -110 #Field edge 110 mm away from top
    for i in range(14):
        iEdge += 5
        MLCEdge.append(iEdge)
    for i in range(32):
        iEdge += 2.5
        MLCEdge.append(iEdge)
        
    for i in range(14):
        iEdge += 5
        MLCEdge.append(iEdge)
        
    return MLCEdge


if __name__ == '__main__':
    #Load the files
  
#    RP = dcm.dcmread("RP.dcm")   
#    RS = dcm.dcmread("RS.dcm")
    RP = dcm.dcmread("RP.CMAT1.dcm")   
    RS = dcm.dcmread("RS.CMAT1.dcm")
    contours = loadContours(RS)
    # for each beam check wether the open fluences overlap with the PTV 
    for beam in RP.BeamSequence:
        #For each control point
        for i,CP in enumerate(beam.ControlPointSequence):
            CP = beam.ControlPointSequence[i]
            #grab information from ControlPoint (CP) structure pertinent to calculation
            if i == 0:
                MLC = CP.BeamLimitingDevicePositionSequence[2]['300a', '011c']
                iso = np.array([CP.IsocenterPosition]).T
                couchAngle = CP.PatientSupportAngle*np.pi/180
                collimatorAngle = CP.BeamLimitingDeviceAngle *np.pi/180
                      
            else:
                MLC = CP.BeamLimitingDevicePositionSequence[0]['300a', '011c']
            
            MLCA = MLC[:60]
            MLCB = MLC[60:]
            
            gantryAngle = -CP.GantryAngle*np.pi/180
            
            #Transform the points to the beam coordinate system
            iContours = transformPoints(contours,iso,gantryAngle, collimatorAngle,couchAngle)
#            print(iContours[0,:],gantryAngle)
            MLCEdges = calculateMLCEdges()
            for j in range(60): #For each leaf pair
                
                #Find Contour points which are in the band of the MLC Leaf Pair
                MLCBandContours = iContours[np.logical_and(iContours.T[2]>MLCEdges[j],iContours.T[2]<MLCEdges[j+1])]
               
                contoursBetweenLeafPairs = MLCBandContours[np.logical_and(MLCBandContours.T[0]>MLCA[j],MLCBandContours.T[0]<MLCB[j])]
                if not(len(contoursBetweenLeafPairs)): #If no points between leafs
                    #Close the leaves (but leave 1mm gap)
                    MLCReplaceValue = (MLCA[j] + MLCB[j])/2
                    MLC.value[j] = MLCReplaceValue - 0.5
                    MLC.value[j+60] = MLCReplaceValue + 0.5 #1mm leaf gap
                                  
            
            if i == 0:
                CP.BeamLimitingDevicePositionSequence[2]['300a', '011c'] = MLC
            else:
                CP.BeamLimitingDevicePositionSequence[0]['300a', '011c'] = MLC
            beam.ControlPointSequence[i] = CP
    
    RP.save_as('FixedMLC1.dcm')

