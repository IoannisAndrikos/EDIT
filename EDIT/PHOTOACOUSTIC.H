#ifndef PHOTOACOUSTIC_H
#define PHOTOACOUSTIC_H

#pragma once

#define _SCL_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <direct.h>
#include <iostream>
#include <numeric>
#include <string>
#include <locale>
#include <codecvt>
#include <experimental/filesystem>
#include "DicomReader.h"


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


class photoAcoustic {
public:
	photoAcoustic(); //constructor
	~photoAcoustic(); //destructor

	enum Point3DType { OXY, DeOXY };

	void exportOXYImages(string dicomPath);
	void exportDeOXYImages(string dicomPath);
	void thicknessExtraction(int frame = 0); //frame = 0 means that the process will be performed for all photoAcoustic sequenence
	void finalizeAllThicknessContours(vector<vector<Point2f>> points);
	void extractOXYandDeOXYPoints(vector<vector<Point2f>> bladderContours, vector<vector<Point2f>> thicknessContours, Point3DType type);
	void extractOXYandDeOXYPoints2(vector<vector<Point2f>> bladderContours, vector<vector<Point2f>> thicknessContours, Point3DType type);
	void writeThicknessPoints();

	void setMainOutputDirectory(string mainOutputDirectory) {
		this->mainOutputDirectory = mainOutputDirectory;
		creatDirectories();
	}

	string getOutputOXYImagesDir() {
		return outputOXYImagesDir;
	}
	
	string getOutputDeOXYImagesDir() {
		return outputDeOXYImagesDir;
	}

	string getStudyDir() {
		return studyDir;
	}

	vector<Mat> getOXYImages() {
		return OXYimages;
	}
	
	vector<Mat> getdeOXYImages() {
		return deOXYimages;
	}

	Mat getIthImage(int i) {
		return OXYimages[i];
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
			this->loggertxt = this->studyDir + separator() + "photoAcoustic_logger.txt";
			this->logFile.open(this->loggertxt);
		}
	}

	void closeLogger() {
		if (this->logFile.is_open()) {
			this->logFile.close();
		}
	}

	vector<vector<Point2f>> getThicknessPoints() {
		return thicknessPoints;
	}

	vector<Point2f> getContourForFix() {
		return this->contourForFix;
	}


	vector<double> getMeanThickness() {
		return meanThickness;
	}

	vector<vector<Point2f>> getFinalThicknessPoints() {
		return finalThicknessPoints;
	}

	vector<vector<Point3f>> getOXYPoints() {
		return OXYPoints;
	}


	vector<vector<Point3f>> getDeOXYPoints() {
		return DeOXYPoints;
	}

	vector<vector<Point3f>> getNotSharderPoints() {
		return notSharderPoints;
	}

	vector<vector<vector<Point3f>>> getSharderPoints() {
		return sharderPoints;
	}
	vector<vector<vector<Point3f>>> getInterpolatedPoints() {
		return interpolatedPoints;
	}

	//------------------------------------------THESE SHOULD BE GIVEN BY ULTRASOUND PART---------------------------
	void setInitialFrame(int initialFrame){
		this->initialFrame = initialFrame;
	}

	void setLastFrame(int lastFrame) {
		this->lastFrame = lastFrame;
	}

	void setlumenPoints(vector<vector<Point2f>> lumenPoints) {
	//	vector<vector<Point2f>>().swap(lumenPoints);
		this->lumenPoints = lumenPoints;
	}

	void setContourForFix(vector<Point2f> contourForFix) {
	//	vector<Point2f>().swap(contourForFix);
		this->contourForFix = contourForFix;
	}

	//-------------------------------------------------------------------------------------------------------------

	void clearThicknessProcess() {
		vector<vector<Point2f>>().swap(this->thicknessPoints);
		vector<vector<Point2f>>().swap(this->lumenPoints);
		vector<vector<Point2f>>().swap(this->finalThicknessPoints);
		vector<Point2f>().swap(this->contourForFix);
		vector<double>().swap(meanThickness);
		LoggerMessage("The Segmentation process was repeated");
	}


	//variables
	bool enableLogging = true;
	ofstream logFile;
	double maxThickness = 0.7;
	double minThickness = 0.3;
	double minDegree = 180.0;
	double maxDegree = 360.0;
	int pixelDistance = 40;
	double distanceBetweenFrames = 0.203;

	Point2f imageCenter;
	double xspace;
	double yspace;
	int fixedMeanPixelThickness = 10;


private:
	//functions
	enum totalSequenceOrCorrecton { TOTAL, CORRECTION};
	void creatDirectories();
	void process(int frame, totalSequenceOrCorrecton type);
	vector<Point2f> smoothContour(vector<Point2f> contour, int num_spline = 0);
	vector<Point2f> smoothCurve(vector<vector<vector<double>>> centerline, int num_spline = 0);
	Point2f findCenterOfContour(vector<Point2f> contour);
	vector<Point3f> findPixelsBetweenThicknessAndBladder(Mat image, vector<Point2f> bladderContours, vector<Point2f> thicknessContours, int iter);
	vector<Point2f> sortUsingPolarCoordinates(vector<Point2f> p, int num_spline);


	void LoggerMessage(string message) {
		if (this->logFile.is_open()) this->logFile << " --> " << message << endl;
	}

	//variables
	string mainOutputDirectory;
	string dicomPath;
	string studyDir;
	string outputOXYImagesDir;
	string outputDeOXYImagesDir;
	string outputSegmentedImagesDir;
	string outputPointsDir;
	string outputPath;
	string filename;
	vector<Mat> OXYimages;
	vector<Mat> deOXYimages;
	vector<vector<Point3f>> OXYPoints, DeOXYPoints;
	vector<Point3f> alreadProcessedFramePoints3D;
	vector<vector<Point3f>> notSharderPoints;
	vector<vector<vector<Point3f>>> sharderPoints;
	vector<vector<vector<Point3f>>> interpolatedPoints;
	vector<double> tags;
	int initialFrame;
	int lastFrame;
	vector<vector<Point2f>> thicknessPoints;
	vector<vector<Point2f>> finalThicknessPoints;
	vector<double> meanThickness;
	vector<vector<Point2f>> lumenPoints;
	vector<Point2f> contourForFix;
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
