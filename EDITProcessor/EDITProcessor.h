#pragma once
#ifndef EDITPROCESSOR_H
#define EDITPROCESSOR_H

#include < vcclr.h >
#include "CVPoint.h"
#include <msclr\marshal_cppstd.h>
//
#include "../EDIT/ultrasound.h"
#include "../EDIT/process_3D.h"
//
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0/vtkPolyData.h>

#include <opencv2/core/core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/opencv.hpp>


#include <fstream>
#include <istream>
#include <iostream>
#include <string>

using namespace cv;
using namespace System;
using namespace System::Windows;
using namespace System::Drawing;
using namespace System::Collections::Generic;
using namespace EDITCore;

namespace EDITProcessor {
	public ref class Processor
	{
	private:
		//----OBJECTS----
		ultrasound *ultr;
		process_3D *proc;
		//----OBJECTS-----
		 
		//-----variables-----
		bool isLoggerEnabled;


		//--------------------------------some functions for internal use---------------------------------
		List<List<EDITCore::CVPoint^>^> ^vectorPointsTOListPoints(vector<vector<Point2f>> lumenPoints) {

			List<List<EDITCore::CVPoint^>^> ^BladderPoints = gcnew List<List<EDITCore::CVPoint^>^>();
			for (int i = 0; i < lumenPoints.size(); i++) {
				List<EDITCore::CVPoint^> ^contour = gcnew List<EDITCore::CVPoint^>();
				for (int j = 0; j < lumenPoints[i].size(); j++) {
					contour->Add(gcnew EDITCore::CVPoint(lumenPoints[i][j].x, lumenPoints[i][j].y));
				}
				BladderPoints->Add(contour);
			}

			return BladderPoints;
		}

		vector<vector<Point2f>> listPointsToVectorPoints(List<List<EDITCore::CVPoint^>^> ^listPoints) {

			vector<vector<Point2f>> vectorPoints;

			for (int i = 0; i < listPoints->Count; i++) {
				vector<Point2f> contour;
				List<EDITCore::CVPoint^> ^listcontour = gcnew List<EDITCore::CVPoint^>();
				listcontour = listPoints[i];
				for (int j = 0; j < listcontour->Count; j++) {
	
					contour.push_back(Point2f(listcontour[j]->GetX(), listcontour[j]->GetY()));
				}
				vectorPoints.push_back(contour);
			}

			return vectorPoints;
		}
		//------------------------------------------------------------------------------------------------

	public:
		Processor() {
			ultr = new ultrasound();
			proc = new process_3D();
		}

		~Processor() {
		}

		

		 void setExaminationsDirectory(System::String ^examDir) {
			string exam_dir = msclr::interop::marshal_as<std::string>(examDir);
			ultr->setMainOutputDirectory(exam_dir);
		}


		 void setSegmentationConfigurations(int repeats, int smoothing, double lamda1, double lamda2, int levelsetSize, bool applyEqualizeHist) {
			 ultr->repeats = repeats;
			 ultr->smoothing = smoothing;
			 ultr->lamda1 = lamda1;
			 ultr->lamda2 = lamda2;
			 ultr->levelsetSize = levelsetSize;
			 ultr->applyEqualizeHist = applyEqualizeHist;
		 }

		 //export images
		 System::String ^exportImages(System::String ^dicomFile, bool isLoggerEnabled) {
			 this->isLoggerEnabled = isLoggerEnabled;
			string Dicom_file = msclr::interop::marshal_as<std::string>(dicomFile);
			ultr->exportImages(Dicom_file); // 1
			if (this->isLoggerEnabled) ultr->openLogger(); //2
			string outputimagesDir = ultr->getOutputImagesDir();
			return msclr::interop::marshal_as<System::String ^>(outputimagesDir);
		}


		 List<double> ^getPixelSpacing() {
			 vector<double> tags = ultr->getTags();
			 List<double>^ pixelSpacing = gcnew List<double>();
			 pixelSpacing->Add(tags[0]*10);
			 pixelSpacing->Add(tags[1]*10);
			 return pixelSpacing;
		 }

		 //extract bladder
		 List<List<EDITCore::CVPoint^>^> ^extractBladder(int startingFrame, int endingFrame, EDITCore::CVPoint ^userPoint) {
			ultr->processing(startingFrame, endingFrame, cv::Point(round(userPoint->GetX()), round(userPoint->GetY())));
			vector<vector<Point2f>> lumenPoints = ultr->getlumenPoints();
			return vectorPointsTOListPoints(lumenPoints);
		 }

		 //extract STL
		 void extractBladderSTL(List<List<EDITCore::CVPoint^>^> ^bladderPoints) {

			 vector<double> Tags = ultr->getTags();
			 proc->xspace = Tags[0] * 10;
			 proc->yspace = Tags[1] * 10;
			 proc->distanceBetweenFrames = 0.203;
			 proc->imageCenter.x = (Tags[3] - Tags[2]) / 2; //center_x = (Xmax - Xmin)/2
			 proc->imageCenter.y = (Tags[5] - Tags[4]) / 2; //center_y = (Ymax - Ymin)/2 
			 string stydyDir = ultr->getStudyDir();
			 proc->setStudyDir(stydyDir);
			 if (this->isLoggerEnabled) proc->openLogger();

			 ultr->finalizePoints(listPointsToVectorPoints(bladderPoints));

			 proc->triangulation(ultr->getlumenPoints(), process_3D::STLType::BLADDER);
			 //proc->triangulation(listPointsToVectorPoints(bladderPoints));
		 }

		 void extractSkin() {
			 ultr->extractSkinPoints();
			 vector<vector<Point2f>> skinPoints = ultr->getSkinPoints();
			 proc->triangulation(skinPoints, process_3D::STLType::SKIN);
		 }


		 void writePointsAndImages() {
			 ultr->writePointsAndImages();
		 }

		 void repeatSegmentation() {
			 ultr->clearLumenAndSkinPoints();
		 }


		 List<double>^ getLumenMetrics() {
			 vector<double> lumenVectorArea = ultr->getLumenArea();
			 List<double>^ lumenListArea;
			 for each (double area in lumenVectorArea)
			 {
				 lumenListArea->Add(area);
			 }
			 vector<double>().swap(lumenVectorArea);

			 return lumenListArea;
		 }

	};
}
#endif