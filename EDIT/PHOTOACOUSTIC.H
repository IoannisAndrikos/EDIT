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
#include "spline.h"

#include "ALGLIB_src/interpolation.h"
#include "BSpline.h"


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

	enum Point3DType { OXY, DeOXY, GNR };
	enum tumorSize { BIG, SMALL };

	string exportOXYImages(string dicomPath);
	string exportDeOXYImages(string dicomPath);
	string exportGNRImages(string dicomPath);
	string thicknessExtraction(int frame = 0); //frame = 0 means that the process will be performed for all photoAcoustic sequenence
	void finalizeAllThicknessContours(vector<vector<Point2f>> points);
	void extractOXYandDeOXYPoints(vector<vector<Point2f>> bladderContours, vector<vector<Point2f>> thicknessContours, Point3DType type);
	void writeThicknessPoints();

	void extractTumorPoints(vector<Mat> tumorImages);


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


	string getOutputGNRImagesDir() {
		return outputGNRImagesDir;
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

	vector<Mat> getGNRImages() {
		return GNRimages;
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
		cout << "KKPP" << contourForFix.size() << endl;

		return this->contourForFix;
	}


	vector<double> getMeanThickness() {
		return this->meanThickness;
	}

	double getMeanThicknessOfContourForFix() {
		return this->meanThicknessOfContourForFix;
	}

	vector<vector<Point2f>> getFinalThicknessPoints() {
		return this->finalThicknessPoints;
	}

	vector<vector<Point3f>> getOXYPoints() {
		return this->OXYPoints;
	}


	vector<vector<Point3f>> getDeOXYPoints() {
		return this->DeOXYPoints;
	}

	vector<vector<Point3f>> getNotSharderPoints() {
		return this->notSharderPoints;
	}

	vector<vector<vector<Point3f>>> getSharderPoints() {
		return this->sharderPoints;
	}
	vector<vector<vector<Point3f>>> getInterpolatedPoints() {
		return this->interpolatedPoints;
	}


	vector<vector<Point3f>> getNotSharderPoints_Tumor() {
		return this->notSharderPoints_Tumor;
	}

	vector<vector<vector<Point3f>>> getSharderPoints_Tumor() {
		return this->sharderPoints_Tumor;
	}
	vector<vector<vector<Point3f>>> getInterpolatedPoints_Tumor() {
		return this->interpolatedPoints_Tumor;
	}


	//------------------------------------------THESE SHOULD BE GIVEN BY ULTRASOUND PART---------------------------
	void setInitialFrame(int initialFrame){
		this->initialFrame = initialFrame;
	}

	void setLastFrame(int lastFrame) {
		this->lastFrame = lastFrame;
	}

	void setlumenPoints(vector<vector<Point2f>> lumenPoints) {
		vector<vector<Point2f>>().swap(this->lumenPoints);
		this->lumenPoints = lumenPoints;
	}

	void setContourForFix(vector<Point2f> contourForFix) {
		vector<Point2f>().swap(this->contourForFix);
		this->contourForFix = contourForFix;
	}

	
	void setThicknessPoints(vector<vector<Point2f>> thicknessPoints) {
		vector<vector<Point2f>>().swap(this->thicknessPoints);
		vector<vector<Point2f>>().swap(this->finalThicknessPoints);

		this->thicknessPoints = thicknessPoints;
		this->finalThicknessPoints = thicknessPoints;
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
	tumorSize tumor = tumorSize::SMALL;
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
	enum interpolationMethod {CUBIC, AKIMA, CATMULLROM, MONOTONE, LINEAR};
	enum totalSequenceOrCorrecton { TOTAL, CORRECTION};
	void creatDirectories();
	string process(int frame, totalSequenceOrCorrecton type);
	vector<Point2f> smoothContour(vector<Point2f> contour, int num_spline = 0, bool closedContour = false);
	vector<Point2f> smoothCurve(vector<vector<vector<double>>> centerline, int num_spline = 0);
	Point2f findCenterOfContour(vector<Point2f> contour);
	vector<Point3f> findPixelsBetweenThicknessAndBladder(Mat image, vector<Point2f> bladderContours, vector<Point2f> thicknessContours, int iter);
	vector<Point2f> sortClockwise(vector<Point2f> p, int num_spline);
	Point2f getCenterOfGravity(vector<Point2f> points);
	Point2f getCenterOfMass(vector<Point2f> points);
	vector<Point2f> interpolateConvexPoints(vector<Point2f> p, interpolationMethod method = interpolationMethod::AKIMA);
	vector<Point2f> findOppositePoints(vector<Point2f> p);
	bool IsClockwise(vector<Point2f> points);

	void LoggerMessage(string message) {
		if (this->logFile.is_open()) this->logFile << " --> " << message << endl;
	}

	//variables
	const string success = "success";
	string mainOutputDirectory;
	string dicomPath;
	string studyDir;
	string outputOXYImagesDir;
	string outputDeOXYImagesDir;
	string outputGNRImagesDir;
	string outputSegmentedImagesDir;
	string outputPointsDir;
	string outputPath;
	string filename;
	vector<Mat> OXYimages;
	vector<Mat> deOXYimages;
	vector<Mat> GNRimages;
	vector<vector<Point3f>> OXYPoints, DeOXYPoints;
	vector<Point3f> alreadProcessedFramePoints3D;
	vector<vector<Point3f>> notSharderPoints;
	vector<vector<vector<Point3f>>> sharderPoints;
	vector<vector<vector<Point3f>>> interpolatedPoints;


	vector<Point3f> alreadProcessedFramePoints3D_Tumor;
	vector<vector<Point3f>> notSharderPoints_Tumor;
	vector<vector<vector<Point3f>>> sharderPoints_Tumor;
	vector<vector<vector<Point3f>>> interpolatedPoints_Tumor;


	vector<double> tags;
	int initialFrame;
	int lastFrame;
	vector<vector<Point2f>> thicknessPoints;
	vector<vector<Point2f>> finalThicknessPoints;
	vector<double> meanThickness;
	double meanThicknessOfContourForFix;
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
