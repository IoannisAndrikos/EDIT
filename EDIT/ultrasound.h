#ifndef ULTRASOUND_H
#define ULTRASOUND_H

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
#include "CImg.h"
#include "spline.h"



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
	ultrasound(string dicomPath, string outputPath);
	~ultrasound(); //destructor

	void exportImages(string dicomPath);

	void processing(int initialFrame, int lastFrame, Point clickPoint); //cropping and filtering //vector<vector<Point2f>>
	void finalizePoints(vector<vector<Point2f>> p);
	void extractSkinPoints();
	void writePointsAndImages();


	void setMainOutputDirectory(string mainOutputDirectory) {
		this->mainOutputDirectory = mainOutputDirectory;
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

	void openLogger() {
		this->loggertxt = this->studyDir + "/ultrasound_logger.txt";
		this->logFile.open(this->loggertxt);
	}

	void closeLogger() {
		this->logFile.close();
	}

	vector<vector<Point2f>> getlumenPoints() {
		return lumenPoints;
	}

	vector<vector<Point2f>> getSkinPoints() {
		return skinPoints;
	}

	void clearLumenAndSkinPoints() {
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
	
private:
	//functions
	enum ResultOfProcess { SUCCESS, FAILURE };
	CImg<unsigned char> *cvImgToCImg(Mat &cvImg);
	Mat CImgtoCvImg(CImg<unsigned char> img);
	void creatDirectories();
	vector<Point2f> smoothCenterline(vector<vector<vector<double>>> centerline, int num_spline);
	CImg<unsigned char> circle_levelset(int height, int width, const array<int, 2>& center, double radius, double scalerow = 1.0);
	vector<vector<vector<double>>> sortBasedEuclideanDistance(vector<Point2f> points);
	ResultOfProcess centerAndPointsOfContour(Mat processed, vector<Point2f> *points, Point2f *center, Point *highestWhitePixel = &Point(0,0));
	void sortUsingPolarCoordinates(vector<Point2f> *p, int iter, Point2f *center, Mat image, int skinDistance);
	int findLongestVector(vector<vector<Point>> vec);


	void LoggerMessage(string message) {
		if (this->logFile.is_open()) this->logFile << " --> " << message << endl;
	}

	//variables
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
	vector<vector<Point2f>> skinPoints;
	CImg<unsigned char> levelsetImage;

	string loggertxt;


};

#endif 