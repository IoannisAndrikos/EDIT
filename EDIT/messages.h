#ifndef MESSAGES_H
#define MESSAGES_H

#pragma once

#define _SCL_SECURE_NO_WARNINGS


namespace warningMessages {

	static const string firtFrameSegmentationFailed = "Failed to segment the first frame. Plese decrease the setting Level Size or select another frame as initial.";

	static const string frameSegmentationFailed = "Failed to segment the frame: ";

	static const string cannotReadTheDicomFile = "Failed to read the selected file. Please check the file exctension is .dcm";

	static const string cannotGetAllNecessaryDicomTags = "The selected Dicom does not contain all the necessary dicom tags for the process. Plese check if the following dicom Tags are available:\n 1) PhysicalDeltaX \n 2) PhysicalDeltaY \n 3) RegionLocationMinX0 \n 4) RegionLocationMaxX1 \n 5) RegionLocationMinY0 \n 6) RegionLocationMaxY1";

	static const string cannotFindThickness = "Cannot find thickness. Please check all configurations and repeat the process";
	
	//----------------------------------------About 3D---------------------------------

	static const string problemDuringTriangulation = "A problem has occurred during the triangulation process. Please check the segmentation for any potential sharp points.";

	static const string cannotProduce3D = "Cannot produce 3D model for only one frame. Please perform segmentation process for more than one frames.";
}

	
#endif 