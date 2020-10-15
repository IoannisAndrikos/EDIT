#ifndef PROCESS_3D_H
#define PROCESS_3D_H

#pragma once

#define _SCL_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <iostream>
#include <direct.h>
#include <string>
#include <locale>
#include <codecvt>
#include <experimental/filesystem>
#include <opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>

// ----------------------------VTK Triangulation----------------------
#include "vtk-9.0/vtkCardinalSpline.h"
#include "vtk-9.0/vtkPoints.h"
#include "vtk-9.0/vtkDoubleArray.h"
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0/vtkPolyData.h>
#include <vtk-9.0/vtkIntersectionPolyDataFilter.h>
#include <vtk-9.0/vtkLine.h>
#include <vtk-9.0/vtkOBBTree.h>
#include <vtk-9.0/vtkModifiedBSPTree.h>
#include <vtk-9.0/vtkCellArray.h>
#include <vtk-9.0/vtkTriangle.h>
#include <vtk-9.0/vtkSTLWriter.h>
#include <vtk-9.0/vtkSmoothPolyDataFilter.h>
#include <vtk-9.0/vtkPolyDataNormals.h>
#include <vtk-9.0/vtkBooleanOperationPolyDataFilter.h>
#include <vtk-9.0/vtkFillHolesFilter.h>
#include <vtk-9.0/vtkPolyDataMapper.h>
#include <vtk-9.0/vtksubdivisionFilter.h>
#include <vtk-9.0/vtkAdaptiveSubdivisionFilter.h>
#include <vtk-9.0/vtkLinearSubdivisionFilter.h>
#include <vtk-9.0/vtkButterflySubdivisionFilter.h>
#include <vtk-9.0/vtkLoopSubdivisionFilter.h>
#include <vtk-9.0/vtkCleanPolyData.h>
#include <vtk-9.0/vtkSelectEnclosedPoints.h>
#include <vtk-9.0/vtkSTLReader.h>
//---------------------------------------------------------------------


using namespace std;
using namespace cv;


class process_3D {
public:
	//functions
	process_3D(); //constructor
	~process_3D(); //destructor

	enum STLType { BLADDER, SKIN, THICKNESS, OXY, DeOXY};

	string triangulation(vector<vector<Point2f>> point_cloud, STLType type); //string path
	string surface_smoothing(vtkSmartPointer<vtkPolyData> surface, STLType type); //string path
	string findPixelsArePlacedIntoGeometries(vector<vector<Point3f>> pixels3D, STLType type);
	string findPixelsArePlacedIntoGeometries2(vector<vector<vector<Point3f>>> sharderPixels, vector<vector<vector<Point3f>>> interpolatedPixels, STLType type);

	void setMainOutputDirectory(string studyDir) {
		this->studyDir = studyDir;

		this->outputObjectsDir = this->studyDir + separator() +  "stl_objects";
		_mkdir(outputObjectsDir.c_str());
	}

	void openLogger(bool open) {
		if (open) {
			this->loggertxt = this->studyDir + separator() + "process_3D_logger.txt";
			this->logFile.open(this->loggertxt);
		}
	}

	void closeLogger() {
		if (this->logFile.is_open()) {
			this->logFile.close();
		}
	}

	//variables
	double xspace;
	double yspace;
	double distanceBetweenFrames;
	ofstream logFile;
	Point2f imageCenter;
	bool fillHoles = true;
	string bladderGeometry;
	string thicknessGeometry;


private:
	//variable
	string studyDir;
	string outputObjectsDir;
	string loggertxt;

	//functions
	string getSTLName(STLType type);

	void LoggerMessage(string message) {
		if (this->logFile.is_open()) this->logFile << " --> " << message << endl;
	}

	inline char separator()
	{
    #ifdef _WIN32
		return '\\';
    #else
		return '/';
    #endif
	}

};

#endif 