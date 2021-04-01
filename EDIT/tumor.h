#ifndef TUMOR_H
#define TUMOR_H

#pragma once

#define _SCL_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <direct.h>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <locale>
#include <codecvt>
#include <experimental/filesystem>
#include "DicomReader.h"
#include "CImg.h"
#include "spline.h"

#include "ALGLIB_src/interpolation.h"

#include <opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/opencv.hpp>


// ----------------------------VTK Triangulation----------------------
#include "vtk-9.0/vtkSmartPointer.h"
#include "vtk-9.0/vtkCardinalSpline.h"
#include "vtk-9.0/vtkPoints.h"
#include "vtk-9.0/vtkPolyData.h"
#include "vtk-9.0/vtkDoubleArray.h"
//---------------------------------------------------------------------

using namespace std;
using namespace cv;
using namespace cimg_library;
using namespace tk;

class tumor {
public:
	tumor(); //constructor
	~tumor(); //destructor

	void sortClockwise(vector<vector<vector<Point2f>>> points);

	vector<vector<Point2f>> getTumorPoints() {
		return this->tumorPoints;
	}


private:
	//methods
	bool IsClockwise(vector<Point2f> points);
	vector<Point2f> smoothContour(vector<Point2f> contour, int num_spline, bool closedContour = false);
	vector<Point2f> smoothCurve(vector<vector<vector<double>>> centerline, int num_spline = 0);

	//variables
	vector<vector<Point2f>> tumorPoints;
	Point2f getCenterOfGravity(vector<Point2f> points);
	Point2f getCenterOfMass(vector<Point2f> points);

	

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