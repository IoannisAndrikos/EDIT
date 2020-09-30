#pragma once
#ifndef EDITPROCESSOR_H
#define EDITPROCESSOR_H

#include < vcclr.h >
#include "CVPoint.h"
#include <msclr\marshal_cppstd.h>
//
#include "../EDIT/ultrasound.h"
#include "../EDIT/process_3D.h"
#include "../EDIT/photoAcoustic.h"
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
		photoAcoustic *photo;
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

		List<EDITCore::CVPoint^>^ vectorPointsTOListPoints(vector<Point2f> lumenPoints) {

			List<EDITCore::CVPoint^>^ BladderPoints = gcnew List<EDITCore::CVPoint^>();
			for (int i = 0; i < lumenPoints.size(); i++) {
				BladderPoints->Add(gcnew EDITCore::CVPoint(lumenPoints[i].x, lumenPoints[i].y));
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

		vector<Point2f> listPointsToVectorPoints(List<EDITCore::CVPoint^>^ listPoints) {

			vector<Point2f> vectorPoints;

			for (int i = 0; i < listPoints->Count; i++) {
				vectorPoints.push_back(Point2f(listPoints[i]->GetX(), listPoints[i]->GetY()));
			}
			return vectorPoints;
		}

		void freeMemory(List<List<EDITCore::CVPoint^>^>^ listPoints) {
			listPoints->Clear();
		}

		void freeMemory(List<EDITCore::CVPoint^>^ listPoints) {
			listPoints->Clear();
		}
		//------------------------------------------------------------------------------------------------

	public:
		Processor() {
			ultr = new ultrasound();
			photo = new photoAcoustic();
			proc = new process_3D();
		}

		~Processor() {
		}

		//------------------------------- U L T R A S O U N D - P A R T ----------------------------------

		 void setExaminationsDirectory(System::String ^examDir) {
			string exam_dir = msclr::interop::marshal_as<std::string>(examDir);
			ultr->setMainOutputDirectory(exam_dir);
			photo->setMainOutputDirectory(exam_dir);
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
		 System::String ^extractBladderSTL(List<List<EDITCore::CVPoint^>^> ^bladderPoints) {

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
			 freeMemory(bladderPoints);

			 string STLPath = proc->triangulation(ultr->getlumenPoints(), process_3D::STLType::BLADDER);
			 return msclr::interop::marshal_as<System::String^>(STLPath);
		 }

		 System::String ^extractSkinSTL() {
			 ultr->extractSkinPoints();
			 vector<vector<Point2f>> skinPoints = ultr->getSkinPoints();
			 string STLPath = proc->triangulation(skinPoints, process_3D::STLType::SKIN);
			 return msclr::interop::marshal_as<System::String^>(STLPath);
		 }


		 void writePointsAndImages() {
			 ultr->writePointsAndImages();
		 }

		 void repeatSegmentation() {
			 ultr->clearUltrasoundProcess();
			 photo->clearThicknessProcess();
		 }

		 // ------------------------------ P H O T O A C O U S T I C - P A R T ----------------------------------

		  //export photaccoustic images
		 System::String^ exportPhotoAcousticImages(System::String^ dicomFile, bool isLoggerEnabled) {
			 this->isLoggerEnabled = isLoggerEnabled;
			 string Dicom_file = msclr::interop::marshal_as<std::string>(dicomFile);
			 photo->exportImages(Dicom_file); // 1
			 if (this->isLoggerEnabled) photo->openLogger(); //2
			 string outputimagesDir = photo->getOutputImagesDir();
			 return msclr::interop::marshal_as<System::String^>(outputimagesDir);
		 }

		 void setPhotoAcousticSegmentationConfigurations(double minThickness, double maxThickness) {
			 photo->minThickness = minThickness;
			 photo->maxThickness = maxThickness;
		 }


		 List<List<EDITCore::CVPoint^>^> ^extractThickness(List<List<EDITCore::CVPoint^>^>^ bladderPoints) {
			// ultr->finalizePoints(listPointsToVectorPoints(bladderPoints));
			
			 photo->setInitialFrame(ultr->getInitialFrame());
			 photo->setLastFrame(ultr->getLastFrame());
			 photo->setlumenPoints(listPointsToVectorPoints(bladderPoints));
			 freeMemory(bladderPoints);
			 photo->thicknessExtraction();
			 vector<vector<Point2f>> thicknessPoints = photo->getThicknessPoints();
			 return vectorPointsTOListPoints(thicknessPoints);
		 }

		 List<EDITCore::CVPoint^>^ extractThicknessForUniqueFrame(int frame, List<EDITCore::CVPoint^>^ bladderPoints) {
			 photo->setContourForFix(listPointsToVectorPoints(bladderPoints));
			 freeMemory(bladderPoints);
			 photo->thicknessExtraction(frame);
			 vector<Point2f> thicknessPoints = photo->getContourForFix();
			 return vectorPointsTOListPoints(thicknessPoints);
		 }

		 //extract STL
		 System::String ^extractThicknessSTL(List<List<EDITCore::CVPoint^>^>^ thicknessPoints) {

			 vector<double> Tags = photo->getTags();
			 proc->xspace = Tags[0] * 10;
			 proc->yspace = Tags[1] * 10;
			 proc->distanceBetweenFrames = 0.203;
			 proc->imageCenter.x = (Tags[3] - Tags[2]) / 2; //center_x = (Xmax - Xmin)/2
			 proc->imageCenter.y = (Tags[5] - Tags[4]) / 2; //center_y = (Ymax - Ymin)/2 
			 string stydyDir = photo->getStudyDir();
			 proc->setStudyDir(stydyDir);
			 if (this->isLoggerEnabled) proc->openLogger();

			 photo->finalizeAllThicknessContours(listPointsToVectorPoints(thicknessPoints));
			 freeMemory(thicknessPoints);

			 string STLPath = proc->triangulation(photo->getFinalThicknessPoints(), process_3D::STLType::THICKNESS);
			 return msclr::interop::marshal_as<System::String^>(STLPath);
		 }


		 List<double>^ getMeanThickness() {
			 vector<double> meanThicknessVec = photo->getMeanThickness();
			 List<double>^ meanThickness = gcnew List<double>();
			 for each (double mT in meanThicknessVec)
			 {
				 meanThickness->Add(mT);
			 }
			 vector<double>().swap(meanThicknessVec);

			 return meanThickness;
		 }


		void writeThicknessPoints(List<List<EDITCore::CVPoint^>^>^ thicknessPoints) {
			photo->finalizeAllThicknessContours(listPointsToVectorPoints(thicknessPoints));
			freeMemory(thicknessPoints);
			photo->writeThicknessPoints();
		 }

	};
}
#endif