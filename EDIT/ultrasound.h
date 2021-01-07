#ifndef ULTRASOUND_H
#define ULTRASOUND_H

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



class ultrasound {
public:
	ultrasound(); //constructor
	~ultrasound(); //destructor

	string exportImages(string dicomPath);

	string processing(int initialFrame, int lastFrame, cv::Point clickPoint); //cropping and filtering //vector<vector<Point2f>>

	string extactTumor2D(Point clickPoint, vector<vector<Point2f>> thickness);
	vector<vector<Point2f>> getTumorBorders();

	string fixArtifact(Point clickPoint, vector<vector<Point2f>> points);
	

	void finalizeAllBladderContours(vector<vector<Point2f>> points);
	void extractSkinPoints(vector<vector<Point2f>> bladderPoints);
	void writePointsAndImages();


	void setMainOutputDirectory(string mainOutputDirectory) {
		this->mainOutputDirectory = mainOutputDirectory;
		creatDirectories();
	}

	string getOutputImagesDir() {
		return outputImagesDir;
	}

	string getStudyDir() {
		return studyDir;
	}

	vector<Mat> getImages() {
		return images;
	}

	Mat getIthImage(int i) {
		return images[i];
	}

	vector<double> getTags() {
		return tags;
	}

	string getDicomPath() {
		return dicomPath;
	}

	string getOutputPath() {
		return outputPath;
	}

	string getFilename() {
		return filename;
	}

	void openLogger(bool open) {
		if (open) {
			this->loggertxt = this->studyDir + separator() + "ultrasound_logger.txt";
			this->logFile.open(this->loggertxt);
		}
	}

	void closeLogger() {
		if (this->logFile.is_open()) {
			this->logFile.close();
		}
	}

	int getInitialFrame() {
		return initialFrame;
	}

	int getLastFrame() {
		return lastFrame;
	}

	void setInitialFrame(int initialFrame) {
		this->initialFrame = initialFrame;
	}

	void setLastFrame(int lastFrame) {
		this->lastFrame = lastFrame;
	}

	vector<vector<Point2f>> getlumenPoints() {
		return lumenPoints;
	}

	vector<Mat> getTumorImages() {
		return this->tumorImages;
	}

	vector<vector<Point3f>> getlumen3DPoints() {
		return lumen3DPoints;
	}

	vector<vector<Point2f>> getSkinPoints() {
		return skinPoints;
	}

	vector<double> getLumenArea() {
		return lumenArea;
	}

	void clearUltrasoundProcess() {
		vector<vector<Point2f>>().swap(this->lumenPoints);
		vector<vector<Point2f>>().swap(this->skinPoints);
		LoggerMessage("The Segmentation process was repeated");
	}


	//variables
	int repeats = 20;
	double smoothing = 3;
	double lamda1 = 1.0;
	double lamda2 = 1.0;
	int levelsetSize = 40;
	boolean applyEqualizeHist = false;
	boolean enableLogging = true;
	ofstream logFile;
	Point2f imageCenter;
	
private:
	//functions
	enum ResultOfProcess { SUCCESS, FAILURE };
	enum interpolationMethod { CUBIC, AKIMA, CATMULLROM, MONOTONE, LINEAR };
	CImg<unsigned char> *cvImgToCImg(Mat &cvImg);
	Mat CImgtoCvImg(CImg<unsigned char> img);
	void creatDirectories();
	vector<Point2f> smoothContour(vector<Point2f> contour, int num_spline, bool closedContour = false);
	vector<Point2f> smoothCurve(vector<vector<vector<double>>> centerline, int num_spline = 0);
	vector<Point2f> interpolateConvexPoints(vector<Point2f> p, interpolationMethod method = interpolationMethod::AKIMA);
	CImg<unsigned char> circle_levelset(int height, int width, const array<int, 2>& center, double radius, double scalerow = 1.0);
	vector<Point2f> sortBasedEuclideanDistance(vector<Point2f> points);
	ResultOfProcess centerAndPointsOfContour(Mat processed, vector<Point2f> *points, Point2f *center, cv::Point *highestWhitePixel = &cv::Point(0,0));
	ResultOfProcess sortUsingPolarCoordinates(vector<Point2f> *p, int iter, Point2f *center, Mat image, int skinDistance);
	ResultOfProcess sortClockwise(vector<Point2f>* p, Point2f* center, int iter);
	int findLongestVector(vector<vector<cv::Point>> vec);
	Point2f getCenterOfGravity(vector<Point2f> points);
	Point2f getCenterOfGravity(vector<Point> points);
	Point2f getCenterOfMass(vector<Point2f> points);
	bool  IsClockwise(vector<Point2f> points);


	void LoggerMessage(string message) {
		if (this->logFile.is_open()) this->logFile << " --> " << message << endl;
	}

	//variables
	const string success = "success";
	string mainOutputDirectory;
	string dicomPath;
	string studyDir;
	string outputImagesDir;
	string outputSegmentedImagesDir;
	string outputPointsDir;
	string outputPath;
	string filename;
	vector<Mat> images;
	vector<double> tags;
	int initialFrame;
	int lastFrame;
	vector<vector<Point2f>> lumenPoints;
	
	vector<Mat> tumorImages;


	vector<vector<Point3f>> lumen3DPoints;
	vector<vector<Point2f>> skinPoints;
	vector<double> lumenArea;
	CImg<unsigned char> levelsetImage;

	string loggertxt;

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