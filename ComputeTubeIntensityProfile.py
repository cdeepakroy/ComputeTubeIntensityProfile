#!/usr/bin/env python

import numpy as np
import SimpleITK as sitk

#
# ComputeTubeIntensityProfile
#

class ComputeTubeIntensityProfile:
    def __init__( self, parent ):
        parent.title = "ComputeTubeIntensityProfile"
        parent.categories = ["Informatics"]
        parent.dependencies = []
        parent.contributors = ["Deepak Roy Chittajallu (Kitware)"]
        parent.helpText = """
        Given a tube's fly through image and mask data, computes the intensity 
        profiles along a cross-section perpendicular to the tube's tangent at 
        each point on the tube's center line, and generates a plot depicting the
        median, lower and upper quartiles of the extracted intesnity profiles.
        """
        self.parent = parent
        
#
# qLComputeTubeIntensityProfileWidget
#

class ComputeTubeIntensityProfileWidget:
    def __init__( self, parent = None ):
        
        #from __main__ import vtk, qt, ctk, slicer
        self.slicer = __import__("slicer", fromlist="__main__")
        self.qt = __import__("qt", fromlist="__main__")
        self.ctk = __import__("ctk", fromlist="__main__")
        
        if not parent:
            self.parent = self.slicer.qMRMLWidget()
            self.parent.setLayout(self.qt.QVBoxLayout())
            self.parent.setMRMLScene(self.slicer.mrmlScene)
        else:
            self.parent = parent

        self.layout = self.parent.layout()
        if not parent:
            self.setup()
            self.parent.show()
            
        self.isDeveloperMode = self.qt.QSettings().value('QtTesting/Enabled')        
        
    def setup(self):
        # Instantiate and connect widgets ...

        #
        # Reload and Test area
        #  
        if self.isDeveloperMode:
            
            reloadCollapsibleButton = self.ctk.ctkCollapsibleButton()
            reloadCollapsibleButton.text = "Reload && Test"
            reloadCollapsibleButton.collapsed = True
            self.layout.addWidget(reloadCollapsibleButton)
            reloadFormLayout = self.qt.QFormLayout(reloadCollapsibleButton)
            
            self.reloadButton = self.qt.QPushButton("Reload")
            self.reloadButton.toolTip = "Reload this module."
            reloadFormLayout.addWidget(self.reloadButton)
            self.reloadButton.connect('clicked()', self.onReload)

        #
        # I/O
        #
        ioCollapsibleButton = self.ctk.ctkCollapsibleButton()
        ioCollapsibleButton.text = "Input/Output"
        ioCollapsibleButton.collapsed = False
        self.layout.addWidget(ioCollapsibleButton)
        
        # Layout within the parameters collapsible button
        ioFormLayout = self.qt.QFormLayout(ioCollapsibleButton)
        
        # input tube fly through image selector
        self.inTubeImageSelector = self.slicer.qMRMLNodeComboBox()
        self.inTubeImageSelector.nodeTypes = (("vtkMRMLScalarVolumeNode"), "")
        self.inTubeImageSelector.addAttribute("vtkMRMLScalarVolumeNode", 
                                              "LabelMap", 0 )
        self.inTubeImageSelector.selectNodeUponCreation = True
        self.inTubeImageSelector.addEnabled = False
        self.inTubeImageSelector.removeEnabled = False
        self.inTubeImageSelector.noneEnabled = False
        self.inTubeImageSelector.showHidden = False
        self.inTubeImageSelector.showChildNodeTypes = False
        self.inTubeImageSelector.setMRMLScene( self.slicer.mrmlScene )
        self.inTubeImageSelector.setToolTip( "Pick the tube fly through \
                                              image." )
        ioFormLayout.addRow("Input Tube Fly-through Image: ", 
                            self.inTubeImageSelector)
        

        # input tube fly through mask image selector
        self.inTubeMaskSelector = self.slicer.qMRMLNodeComboBox()
        self.inTubeMaskSelector.nodeTypes = (("vtkMRMLScalarVolumeNode"), "")
        self.inTubeMaskSelector.addAttribute("vtkMRMLScalarVolumeNode", 
                                             "LabelMap", 1 )
        self.inTubeMaskSelector.selectNodeUponCreation = True
        self.inTubeMaskSelector.addEnabled = False
        self.inTubeMaskSelector.removeEnabled = False
        self.inTubeMaskSelector.noneEnabled = False
        self.inTubeMaskSelector.showHidden = False
        self.inTubeMaskSelector.showChildNodeTypes = False
        self.inTubeMaskSelector.setMRMLScene(self.slicer.mrmlScene)
        self.inTubeMaskSelector.setToolTip("Pick the tube fly through \
                                            mask image." )
        ioFormLayout.addRow("Input Tube Fly-through Mask: ", 
                            self.inTubeMaskSelector)
        
        #
        # output tube intensity profile array selector
        #
        self.outTubeProfileSelector = self.slicer.qMRMLNodeComboBox()
        self.outTubeProfileSelector.nodeTypes = (("vtkMRMLDoubleArrayNode"), "")
        self.outTubeProfileSelector.addEnabled = True
        self.outTubeProfileSelector.removeEnabled = True
        self.outTubeProfileSelector.noneEnabled = False
        self.outTubeProfileSelector.showHidden = False
        self.outTubeProfileSelector.showChildNodeTypes = False
        self.outTubeProfileSelector.setMRMLScene(self.slicer.mrmlScene)
        self.outTubeProfileSelector.setToolTip( "Pick the output tube \
                                                 intensity profile array." )
        ioFormLayout.addRow("Output Tube Intensity Profile: ", 
                            self.outTubeProfileSelector)

        #
        # Apply Button
        #
        self.applyButton = self.qt.QPushButton("Apply")
        self.applyButton.toolTip = "Run the algorithm."
        self.applyButton.enabled = False
        ioFormLayout.addRow(self.applyButton)
        self.onSelect()
    
        # connections
        self.applyButton.connect('clicked(bool)', self.onApplyButton)
        self.inTubeImageSelector.connect("currentNodeChanged(vtkMRMLNode*)", 
                                         self.onSelect )
        self.inTubeMaskSelector.connect("currentNodeChanged(vtkMRMLNode*)", 
                                        self.onSelect )
        self.outTubeProfileSelector.connect("currentNodeChanged(vtkMRMLNode*)", 
                                            self.onSelect)
    
        # Add vertical spacer
        self.layout.addStretch(1)

    def onSelect(self):
        self.applyButton.enabled = self.inTubeImageSelector.currentNode() and \
                                   self.inTubeMaskSelector.currentNode() and \
                                   self.outTubeProfileSelector.currentNode()

    def onApplyButton(self):       
        logic = ComputeTubeIntensityProfileLogic()
        logic.run(self.inTubeImageSelector.currentNode(), 
                  self.inTubeMaskSelector.currentNode(),
                  self.outTubeProfileSelector.currentNode())

    def onReload(self,moduleName="ComputeTubeIntensityProfile"):
        """Generic reload method for any scripted module.
        ModuleWizard will subsitute correct default moduleName.
        """
        globals()[moduleName] = self.slicer.util.reloadScriptedModule(moduleName)

#
# ComputeTubeIntensityProfileLogic
#

class ComputeTubeIntensityProfileLogic:
    """This class should implement all the actual 
    computation done by your module.  The interface 
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget
    """
    def __init__(self):
        self.chartNodeID = None
        self.slicer = __import__("slicer", fromlist="__main__")
        self.vtk = __import__("vtk", fromlist="__main__")
        pass

    def run(self, inTubeImage, inTubeMask, outTubeProfile):
        """
        Run the actual algorithm
        """

        def vtk_image_to_sitk_image(imVTK):
            
            dims = imVTK.GetDimensions()
            
            imArr = self.vtk.util.numpy_support.vtk_to_numpy(imVTK.GetPointData().GetScalars())
            imArr = imArr.reshape(dims[2], dims[1], dims[0])
            
            return sitk.GetImageFromArray(imArr)
            
        # convert images to sitk format
        imTubeSitk = vtk_image_to_sitk_image(inTubeImage.GetImageData())
        imTubeMaskSitk = vtk_image_to_sitk_image(inTubeMask.GetImageData())                
                
        # compute tube intensity profiles    
        (tubeIntensityProfiles, 
         nmzdDistToTubeCenter) = ComputeTubeIntensityProfiles(imTubeSitk,
                                                              imTubeMaskSitk)
                                         
        # compute profile stats  
        profileMedian = np.median(tubeIntensityProfiles, 0)
        
        # update output array
        outTubeProfileArr = outTubeProfile.GetArray()
        outTubeProfileArr.SetNumberOfTuples( np.size(profileMedian) )     
        
        for i in range(len(profileMedian)):
            outTubeProfileArr.SetComponent(i, 0, nmzdDistToTubeCenter[i])
            outTubeProfileArr.SetComponent(i, 1, profileMedian[i])
            outTubeProfileArr.SetComponent(i, 2, 0)
        
        # update chart
        ln = self.slicer.util.getNode(pattern='vtkMRMLLayoutNode*')
        ln.SetViewArrangement(24)

        cn = None
        if self.chartNodeID:
            cn = self.slicer.mrmlScene.GetNodeByID(self.chartNodeID)
            
        if not cn:
            cn = self.slicer.mrmlScene.AddNode(self.slicer.vtkMRMLChartNode())
            self.chartNodeID = cn.GetID()
            cn.SetProperty('default', 'title', 
                           'Median Tube Intensity profile')
            cn.SetProperty('default','xAxisLabel', 
                           'Distance to tube center / max tube radius')
            cn.SetProperty('default', 'yAxisLabel', 
                           'Intensity (HU)')
                           
        cn.AddArray(inTubeImage.GetName(), outTubeProfile.GetID())                           
        
        cvn = self.slicer.util.getNode(pattern="vtkMRMLChartViewNode*")
        cvn.SetChartNodeID(cn.GetID())
        cvn.Modified()
        
def ComputeTubeIntensityProfileStats(intensityProfiles):
    
    profileMedian = np.percentile( intensityProfiles, 50, 0 )
    profileLowerQuartile = np.percentile( intensityProfiles, 25, 0 )
    profileUpperQuartile = np.percentile( intensityProfiles, 75, 0 )        
    
    return (profileMedian, profileLowerQuartile, profileUpperQuartile)
    
    
def ComputeTubeIntensityProfiles(imTubeSitk, 
                                 imTubeMaskSitk):
    """
    Given a tube's fly through image and mask as sitk image objects, this 
    function computes the intensity profiles along a cross-section 
    perpendicular to the tube's tangent at each point on the tube's 
    center line
    
    Usage:
    
    (intensityProfiles, 
     distToTubeCenter) = ComputeTubeIntensityProfiles(imTubeSitk, 
                                                      imTubeMaskSitk)
    Input:

        imTubeSitk:  
        
            Tube fly through image as an sitk object
        
        tubeFlyThroughMaskFile:   
        
            Tube fly through mask as an sitk object where in non-zero pixels
            denote the pixels inside the tube.
    
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

    # Grab input
    imTube = sitk.GetArrayFromImage( imTubeSitk )
    imTubeMask = sitk.GetArrayFromImage( imTubeMaskSitk )

    imageSize = imTubeSitk.GetSize()

    # Extract mid-y aka coronal cross section    
    if np.size( imageSize ) == 3:
        imCrossSecTube = imTube[:,np.uint(0.5*imageSize[1]),:]
        imCrossSecTubeMask = imTubeMask[:,np.uint(0.5*imageSize[1]),:]
    elif np.size( imageSize ) == 2:
        imCrossSecTube = imTube
        imCrossSecTubeMask = imTubeMask        

    # compute max tube radius
    maxTubeDiameter = 0
    
    for i in range(imCrossSecTubeMask.shape[0]):
        
        nzInd = imCrossSecTubeMask[i,:].nonzero()[0]
        if np.size(nzInd) == 0:
            continue;
            
        curDiameter = nzInd[-1] - nzInd[0] + 1
        
        if curDiameter > maxTubeDiameter:
            maxTubeDiameter = curDiameter

    # extract intensity profile of each row in the cross section image      
    tubeIntensityProfiles = np.zeros((imCrossSecTubeMask.shape[0], 
                                      maxTubeDiameter))
                                  
    for i in range(imCrossSecTubeMask.shape[0]):
        
        # extract the part of the current profile inside the tube
        nzInd = imCrossSecTubeMask[i,:].nonzero()[0]
        if np.size(nzInd) == 0:
            continue;
        
        curTubeIntensities = imCrossSecTube[i, nzInd[0]:nzInd[-1]]
        minVal = curTubeIntensities.min()
        
        # pad values
        curTubeIntensities = np.append(curTubeIntensities, minVal)
        curTubeIntensities = np.insert(curTubeIntensities, 0, minVal)
        
        # resize to maxTubeDiameter
        rzFactor = np.size(curTubeIntensities) / np.float(maxTubeDiameter)
        
        resampler = sitk.ResampleImageFilter()
        resampler.SetTransform(sitk.Transform())
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetSize([int(maxTubeDiameter), 2])
        resampler.SetOutputSpacing( [rzFactor, rzFactor])  
        resampler.SetDefaultPixelValue(minVal.astype('double'))
        curProfile = resampler.Execute(sitk.GetImageFromArray(np.tile(curTubeIntensities,(2,1))) )
        curProfile = sitk.GetArrayFromImage(curProfile)[0,:]
        
        tubeIntensityProfiles[i, :] = curProfile 
        
    nmzdDistToTubeCenter = (np.arange(maxTubeDiameter) - (maxTubeDiameter * 0.5)) / (maxTubeDiameter * 0.5)

    return (tubeIntensityProfiles, nmzdDistToTubeCenter)
        
def main():    
    
    import argparse
    import matplotlib.pyplot as plt
    
    description = """
    Given a tube's fly through image and mask files, computes the intensity 
    profiles along a cross-section perpendicular to the tube's tangent at 
    each point on the tube's center line, and generates a plot depicting the
    median, lower and upper quartiles of the extracted intesnity profiles, 
    and displays and/or saves the plot to a file.
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
                              type = argparse.FileType('w'),
                              help = "Path to a file in which to save \
                                      tube intensity profile plot. If \
                                      this argument is not provided then \
                                      the plot will be displayed onto the \
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

    # read tube fly through image file
    imTubeSitk = sitk.ReadImage( args.tubeFlyThroughImageFile )

    # read tube mask file    
    imTubeMaskSitk = sitk.ReadImage( args.tubeFlyThroughMaskFile )
    
    # compute tube intensity profiles    
    (tubeIntensityProfiles, 
     nmzdDistToTubeCenter) = ComputeTubeIntensityProfiles(imTubeSitk,
                                                          imTubeMaskSitk )
                                     
    # compute profile stats                                         
    (profileMedian, 
     profileLowerQuartile, 
     profileUpperQuartile) = ComputeTubeIntensityProfileStats( 
                                                 tubeIntensityProfiles )                                    

    # display
    plt.plot(nmzdDistToTubeCenter, profileMedian, linewidth=2.0)
    
    plt.fill_between(nmzdDistToTubeCenter, 
                     profileLowerQuartile, profileUpperQuartile, 
                     color='g', alpha=0.2)
                     
    plt.plot([], [], color='g', alpha=0.2, linewidth=10)
    
    plt.xlabel('Distance to tube center / max tube radius')
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