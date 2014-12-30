#!/usr/bin/env python

import numpy as np
import SimpleITK as sitk
import scipy.ndimage

def ComputeTubeIntensityProfileStats(intensityProfiles):
    
    profileMedian = np.percentile( intensityProfiles, 50, 0 )
    profileLowerQuartile = np.percentile( intensityProfiles, 25, 0 )
    profileUpperQuartile = np.percentile( intensityProfiles, 75, 0 )        
    
    return (profileMedian, profileLowerQuartile, profileUpperQuartile)
    
    
def ComputeTubeIntensityProfiles(tubeFlyThroughImageFile, 
                                 tubeFlyThroughMaskFile):
    """
    Given a tube's fly through image and mask files, this function computes
    the intensity profiles along a cross-section perpendicular to the tube's 
    tangent at each point on the tube's center line
    
    Usage:
    
    (intensityProfiles, 
     distToTubeCenter) = ComputeTubeIntensityProfiles(tubeFlyThroughImageFile, 
                                                      tubeFlyThroughMaskFile)
    Input:

        tubeFlyThroughImageFile:  
        
            Path to an image file containing the fly through 
            image data of a tube.
        
        tubeFlyThroughMaskFile:   
        
            Path to an image file where in the non-zero pixels
            denote the pixels inside the tube in the tube's fly
            through image file.
    
    Output:
    
        intensityProfiles: 
        
            A 2D matrix where in the i'th row denotes the 
            cross-sectional intensity profile at the i'th 
            point on the tube's center line
                           
        distToTubeCenter:  
        
            A 1D array where in the ith element contains the
            distance to the tube's center of the i'th column
            in the intensityProfiles matrix        
    """

    # load tube image     
    imTubeSitk = sitk.ReadImage( tubeFlyThroughImageFile )
    imTube = sitk.GetArrayFromImage( imTubeSitk )

    imageSize = imTubeSitk.GetSize()
    imageSpacing = imTubeSitk.GetSpacing()        
    
    # load tube mask
    imTubeMaskSitk = sitk.ReadImage( tubeFlyThroughMaskFile )
    imTubeMask = sitk.GetArrayFromImage( imTubeMaskSitk )

    # Extract mid-y aka coronal cross section    
    if np.size( imageSize ) == 3:
        imCrossSecTube = imTube[:,np.uint(0.5*imageSize[1]),:]
        imCrossSecTubeMask = imTubeMask[:,np.uint(0.5*imageSize[1]),:]
    elif np.size( imageSize ) == 2:
        imCrossSecTube = imTube
        imCrossSecTubeMask = imTubeMask        

    # compute max tube radius
    maxRadius = 0
    
    for i in range(imCrossSecTubeMask.shape[0]):
        
        nzInd = imCrossSecTubeMask[i,:].nonzero()[0]
        if np.size(nzInd) == 0:
            continue;
            
        curRad = nzInd[-1] - nzInd[0]
        
        if curRad > maxRadius:
            maxRadius = curRad

    # extract intensity profile of each row in the cross section image       
    tubeIntensityProfiles = np.zeros((imCrossSecTubeMask.shape[0], 
                                  maxRadius))
                                  
    interpolationSplineOrder = 3 # cubic spline interpolation
    
    distToTubeCenter = (np.arange(maxRadius) - maxRadius/2) * imageSpacing[0]
    
    for i in range(imCrossSecTubeMask.shape[0]):
        
        nzInd = imCrossSecTubeMask[i,:].nonzero()[0]
        if np.size(nzInd) == 0:
            continue;
    
        curTubeIntensities = imCrossSecTube[i, nzInd[0]:nzInd[-1]]
        
        curZoom = np.float( maxRadius ) / np.size(curTubeIntensities) 
    
        curProfile = scipy.ndimage.zoom(curTubeIntensities, 
                                        zoom=curZoom, 
                                        order=interpolationSplineOrder)
        
        tubeIntensityProfiles[i, :] = curProfile    
        

    return (tubeIntensityProfiles, distToTubeCenter)
        
def main():    
    
    import os
    import argparse
    import matplotlib.pyplot as plt
    
    description = """
    Given a tube's fly through image and mask files, computes the intensity 
    profiles along a cross-section perpendicular to the tube's tangent at 
    each point on the tube's center line, and generates a plot depicting the
    median, lower and upper quartiles of the extracted intesnity profiles.
    """

    
    inputParser = argparse.ArgumentParser(description=description)
    
    inputParser.add_argument("tubeFlyThroughImageFile",
                              help = "Path to an image file containing the \
                                      fly through image data of a tube.")
                             

    inputParser.add_argument("tubeFlyThroughMaskFile",
                              help = "Path to an image file where in the \
                                      non-zero values denote pixels \
                                      inside the tube in the tube's fly \
                                      through image file.")
                                      
    inputParser.add_argument("outputTubeIntensityProfilePlotFile", 
                              nargs = "?",
                              help = "Path to a file in which to save \
                                      tube intensity profile plot. If \
                                      this argument is not provided then \
                                      the plot will be displayed to the \
                                      screen")                                      

    inputParser.add_argument("-t", "--title", 
                             default = "Tube intensity profile",
                             metavar = "<string>",
                             help = "Title of the intensity profile plot") 

    inputParser.add_argument("-d", "--display", 
                             default = False,
                             action = 'store_true',
                             help = "Display the intensity profile plot.") 
                             
    args = inputParser.parse_args()

    # compute tube intensity profiles    
    (tubeIntensityProfiles, 
     distToTubeCenter) = ComputeTubeIntensityProfiles(
                                     args.tubeFlyThroughImageFile,
                                     args.tubeFlyThroughMaskFile )
                                     
    # compute profile stats                                         
    profileMedian = np.percentile( tubeIntensityProfiles, 50, 0 )
    profileLowerQuartile = np.percentile( tubeIntensityProfiles, 25, 0 )
    profileUpperQuartile = np.percentile( tubeIntensityProfiles, 75, 0 )   

    # display
    plt.plot(distToTubeCenter, profileMedian, linewidth=2.0)
    
    plt.fill_between(distToTubeCenter, 
                     profileLowerQuartile, profileUpperQuartile, 
                     color='g', alpha=0.2)
                     
    plt.plot([], [], color='g', alpha=0.2, linewidth=10)
    
    plt.xlabel('Distance from tube center (mm)')
    plt.ylabel('Intensity (HU)')    
    plt.legend(['Median', 'IQR'])
    plt.title(args.title)
    plt.grid(True)    
    
    if args.outputTubeIntensityProfilePlotFile is not None:
        plt.savefig(args.outputTubeIntensityProfilePlotFile)
        
        if args.display:
            plt.show()            
    else:
        plt.show()        
    
if __name__ == "__main__":
    main()