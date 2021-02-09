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
		EDITResponse ^response;
		//----OBJECTS-----
		 
		//-----variables-----
		bool isLoggerEnabled = true;

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

		void releaseMemory(List<List<EDITCore::CVPoint^>^>^ listPoints) {
			listPoints->Clear();
		}

		void releaseMemory(List<EDITCore::CVPoint^>^ listPoints) {
			listPoints->Clear();
		}

		void releaseMemory(vector<vector<Point2f>> points) {
			vector<vector<Point2f>>().swap(points);
		}
		
		void releaseMemory(vector<Point2f> points) {
			vector<Point2f>().swap(points);
		}
		
		void releaseMemory(vector<double> values) {
			vector<double>().swap(values);
		}
		//------------------------------------------------------------------------------------------------

	public:
		Processor() {
			ultr = new ultrasound();
			photo = new photoAcoustic();
			proc = new process_3D();
			response = EDITResponse::Instance;
		}

		~Processor() {
		}


		void setLoggingOnOff(bool OnOff) {
			if (OnOff) {
				ultr->openLogger(OnOff);
				photo->openLogger(OnOff);
				proc->openLogger(OnOff);
			}
			else if (!OnOff) {
				ultr->closeLogger();
				photo->closeLogger();
				proc->closeLogger();
			}
		}

		//------------------------------- U L T R A S O U N D - P A R T ----------------------------------

		 void setExaminationsDirectory(System::String ^examDir) {
			string exam_dir = msclr::interop::marshal_as<std::string>(examDir);
			ultr->setMainOutputDirectory(exam_dir);
			photo->setMainOutputDirectory(exam_dir);
			proc->setMainOutputDirectory(exam_dir);

			setLoggingOnOff(false); //enable logging by default
		}

		
		 void setSegmentationConfigurations(int repeats, int smoothing, double lamda1, double lamda2, int levelsetSize, bool applyEqualizeHist) {
			 ultr->repeats = repeats;
			 ultr->smoothing = smoothing;
			 ultr->lamda1 = lamda1;
			 ultr->lamda2 = lamda2;
			 ultr->levelsetSize = levelsetSize;
			 ultr->applyEqualizeHist = applyEqualizeHist;
		 }

		 void setStudySettings(double distanceBetweenFrames, double xspace, double yspace) {
			 proc->distanceBetweenFrames = distanceBetweenFrames;
			 photo->distanceBetweenFrames = distanceBetweenFrames;

			 proc->xspace = xspace;
			 proc->yspace = xspace;

			 photo->xspace = yspace;
			 photo->yspace = yspace;
		 }


		 void setDicomTags() {
			 try {
				 vector<double> Tags = ultr->getTags();
				 proc->xspace = Tags[0] * 10;
				 proc->yspace = Tags[1] * 10;
				 proc->imageCenter.x = (Tags[3] - Tags[2]) / 2; //center_x = (Xmax - Xmin)/2
				 proc->imageCenter.y = (Tags[5] - Tags[4]) / 2; //center_y = (Ymax - Ymin)/2 


				 photo->xspace = Tags[0] * 10;
				 photo->yspace = Tags[1] * 10;
				 photo->imageCenter.x = (Tags[3] - Tags[2]) / 2; //center_x = (Xmax - Xmin)/2
				 photo->imageCenter.y = (Tags[5] - Tags[4]) / 2; //center_y = (Ymax - Ymin)/2 

				 response->setSuccessOrFailure("success");
			 }
			 catch (exception e) {
				 response->setSuccessOrFailure("Cannot load Dicom Tags");
			 } 
		 }


		 //export images
		 void exportUltrasoundImages(System::String ^dicomFile) {
			string Dicom_file = msclr::interop::marshal_as<std::string>(dicomFile);
			string errorMessage = ultr->exportImages(Dicom_file); // 1
			response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			if (response->isSuccessful()) {
				setDicomTags();
				string outputimagesDir = ultr->getOutputImagesDir();
				vector<double> tags = ultr->getTags();
				//Pixel Spacing
				response->addNumericData(tags[0] * 10);
				response->addNumericData(tags[1] * 10);
				//image size
				response->addNumericData(tags[3] - tags[2]);
				response->addNumericData(tags[5] - tags[4]);
				releaseMemory(tags);
				response->setData(msclr::interop::marshal_as<System::String^>(outputimagesDir));
			}
		}


		 //extract bladder
		 void extractBladder(int startingFrame, int endingFrame, EDITCore::CVPoint ^userPoint, bool fixArtifact) {

			 vector<double> Tags = ultr->getTags();
			 ultr->imageCenter.x = (Tags[3] - Tags[2]) / 2; //center_x = (Xmax - Xmin)/2
			 ultr->imageCenter.y = (Tags[5] - Tags[4]) / 2;

			string errorMessage = ultr->processing(startingFrame, endingFrame, cv::Point(round(userPoint->GetX()), round(userPoint->GetY())));
			response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			if (response->isSuccessful()) {
				if (fixArtifact) {
					errorMessage = ultr->fixArtifact(cv::Point(round(userPoint->GetX()), round(userPoint->GetY())), ultr->getlumenPoints());
					response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
					if (response->isSuccessful()) {
						response->setData(vectorPointsTOListPoints(ultr->getlumenPoints()));
					}
				}
				else {
					response->setData(vectorPointsTOListPoints(ultr->getlumenPoints()));
				}
			}
		 }

		 //extract bladder
		 void extractBladderForUniqueFrame(int frame, List<EDITCore::CVPoint^>^ bladderPoints, bool fixArtifact) {
			string errorMessage = ultr->recalculate(frame, listPointsToVectorPoints(bladderPoints));
			releaseMemory(bladderPoints);
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 if (fixArtifact) {
					 errorMessage = ultr->fixArtifact(frame, ultr->getContourForFix());
					 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
					 if (response->isSuccessful()) {
						 response->setData(vectorPointsTOListPoints(ultr->getContourForFix()));
					 }
				 }
				 else {
					 response->setData(vectorPointsTOListPoints(ultr->getContourForFix()));
				 }
				
			 }
		 }

		/* void fixArtifact(EDITCore::CVPoint^ userPoint, List<List<EDITCore::CVPoint^>^>^ bladderPoints) {
			 string errorMessage = ultr->fixArtifact(cv::Point(round(userPoint->GetX()), round(userPoint->GetY())), listPointsToVectorPoints(bladderPoints));

			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 response->setData(vectorPointsTOListPoints(ultr->getlumenPoints()));
			 }
		 }*/


		 //extract STL
		 void extractBladderSTL(List<List<EDITCore::CVPoint^>^>^ bladderPoints, bool fillHoles) {
			// setDicomTags();
			 proc->fillHoles = fillHoles;
			
			 ultr->finalizeAllBladderContours(listPointsToVectorPoints(bladderPoints));
			 releaseMemory(bladderPoints);
			 string errorMessage = proc->triangulation(ultr->getlumenPoints(), process_3D::STLType::BLADDER);
			
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 string STLPath = proc->getBladderGeometry();
				 response->setData(msclr::interop::marshal_as<System::String^>(STLPath));
			 }
		 }

		 void extractSkinSTL(List<List<EDITCore::CVPoint^>^>^ bladderPoints, bool fillHoles) {
			// setDicomTags();
			 proc->fillHoles = fillHoles;
			 ultr->extractSkinPoints(listPointsToVectorPoints(bladderPoints));
			 releaseMemory(bladderPoints);
			
			 string errorMessage = proc->triangulation(ultr->getSkinPoints(), process_3D::STLType::SKIN);

			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 string STLPath = proc->getSkinGeometry();
				 response->setData(msclr::interop::marshal_as<System::String^>(STLPath));
			 }
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
		 void exportOXYImages(System::String^ dicomFile) {
			 string Dicom_file = msclr::interop::marshal_as<std::string>(dicomFile);
			 string errorMessage = photo->exportOXYImages(Dicom_file); // 
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 vector<double> tags = photo->getTags();
				 //Pixel Spacing
				 response->addNumericData(tags[0] * 10);
				 response->addNumericData(tags[1] * 10);
				 //image size
				 response->addNumericData(tags[3] - tags[2]);
				 response->addNumericData(tags[5] - tags[4]);
				 releaseMemory(tags);
				 string outputimagesDir = photo->getOutputOXYImagesDir();
				 response->setData(msclr::interop::marshal_as<System::String^>(outputimagesDir));
			 }
		 }

		 //export photaccoustic images
		 void exportDeOXYImages(System::String^ dicomFile) {
			 string Dicom_file = msclr::interop::marshal_as<std::string>(dicomFile);
			 string errorMessage = photo->exportDeOXYImages(Dicom_file); // 1
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 vector<double> tags = photo->getTags();
				 //Pixel Spacing
				 response->addNumericData(tags[0] * 10);
				 response->addNumericData(tags[1] * 10);
				 //image size
				 response->addNumericData(tags[3] - tags[2]);
				 response->addNumericData(tags[5] - tags[4]);
				 releaseMemory(tags);
				 string outputimagesDir = photo->getOutputDeOXYImagesDir();
				 response->setData(msclr::interop::marshal_as<System::String^>(outputimagesDir));
			 }
		 }

		 void setPhotoAcousticSegmentationConfigurations(double minThickness, double maxThickness, bool bigTumor) {
			 photo->minThickness = minThickness;
			 photo->maxThickness = maxThickness;
			 (bigTumor) ? photo->tumor = photoAcoustic::tumorSize::BIG : photo->tumor = photoAcoustic::tumorSize::SMALL;
		 }


		 void extractThickness(List<List<EDITCore::CVPoint^>^>^ bladderPoints) {
			 //setDicomTags();
			 if (response->isSuccessful()) { //ckeck dico Tags
				 photo->setInitialFrame(ultr->getInitialFrame());
				 photo->setLastFrame(ultr->getLastFrame());
				 photo->setlumenPoints(listPointsToVectorPoints(bladderPoints));
				 releaseMemory(bladderPoints);
				 string errorMessage = photo->thicknessExtraction();
				 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
				 if (response->isSuccessful()) {
					 vector<vector<Point2f>> thicknessPoints = photo->getThicknessPoints();

					 vector<double> meanThicknessVec = photo->getMeanThickness();
					 List<double>^ meanThickness = gcnew List<double>();
					 for each (double mT in meanThicknessVec)
					 {
						 meanThickness->Add(mT);
					 }
					 releaseMemory(meanThicknessVec);
					 response->setData(meanThickness); //mean thickness
					 response->setData(vectorPointsTOListPoints(thicknessPoints));
				 }
			 }
		 }

		 void extractThicknessForUniqueFrame(int frame, List<EDITCore::CVPoint^>^ bladderPoints) {
			 photo->setContourForFix(listPointsToVectorPoints(bladderPoints));
			 releaseMemory(bladderPoints);
			 string errorMessage = photo->thicknessExtraction(frame);
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 vector<Point2f> thicknessPoints = photo->getContourForFix();
				 double contourMeanThickness = photo->getMeanThicknessOfContourForFix();
				 response->addNumericData(contourMeanThickness);
				 response->setData(vectorPointsTOListPoints(thicknessPoints));
			 }
		 }


		 //extract STL
		 void extractThicknessSTL(List<List<EDITCore::CVPoint^>^>^ thicknessPoints, bool fillHoles) {
			 vector<double> Tags = photo->getTags();
			 proc->fillHoles = fillHoles;
			 //setDicomTags();
			 photo->finalizeAllThicknessContours(listPointsToVectorPoints(thicknessPoints));
			 releaseMemory(thicknessPoints);
			 string errorMessage = proc->triangulation(photo->getFinalThicknessPoints(), process_3D::STLType::THICKNESS);
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 string STLPath = proc->getThicknessGeometry();
				 response->setData(msclr::interop::marshal_as<System::String^>(STLPath));
			 }
		 }


		 void extractOXYandDeOXYPoints(List<List<EDITCore::CVPoint^>^>^ bladderPoints, List<List<EDITCore::CVPoint^>^>^ thicknessPoints, System::String^ BladderFilePath, System::String^ thicknessFilePath) {
			
			
			proc->saveGeometryPath(msclr::interop::marshal_as<std::string>(BladderFilePath), process_3D::STLType::BLADDER);
			proc->saveGeometryPath(msclr::interop::marshal_as<std::string>(thicknessFilePath), process_3D::STLType::THICKNESS);
		
			photo->extractOXYandDeOXYPoints(listPointsToVectorPoints(bladderPoints), listPointsToVectorPoints(thicknessPoints), photoAcoustic::Point3DType::OXY);
			string errorMessage1 = proc->findPixelsArePlacedIntoGeometries(photo->getSharderPoints(), photo->getInterpolatedPoints(), process_3D::STLType::OXY);
			photo->extractOXYandDeOXYPoints(listPointsToVectorPoints(bladderPoints), listPointsToVectorPoints(thicknessPoints), photoAcoustic::Point3DType::DeOXY);
			string errorMessage2 = proc->findPixelsArePlacedIntoGeometries(photo->getSharderPoints(), photo->getInterpolatedPoints(), process_3D::STLType::DeOXY);

			response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage2));
			if (response->isSuccessful()) {
				List<System::String^>^ paths = gcnew List<System::String^>();
				paths->Add(msclr::interop::marshal_as<System::String^>(proc->getOXYGeometry()));
				paths->Add(msclr::interop::marshal_as<System::String^>(proc->getDeOXYGeometry()));
				response->setData(paths);
			}
		 }



		 //extract bladder
		 void extractTumor2D(EDITCore::CVPoint^ userPoint, List<List<EDITCore::CVPoint^>^>^ bladderPoints, List<List<EDITCore::CVPoint^>^>^ thicknessPoints, System::String^ BladderFilePath, System::String^ thicknessFilePath) {
 
			 string errorMessage = ultr->extactTumor2D(cv::Point(round(userPoint->GetX()), round(userPoint->GetY())), listPointsToVectorPoints(bladderPoints), listPointsToVectorPoints(thicknessPoints));
			 releaseMemory(thicknessPoints);
			 releaseMemory(bladderPoints);

			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 response->setData(vectorPointsTOListPoints(ultr->getTumorBorders()));
			 }
		 }


		 void extractTumor3D(){
			 photo->extractTumorPoints(ultr->getTumorImages());
			 string errorMessage = proc->findPixelsArePlacedIntoGeometries(photo->getSharderPoints_Tumor(), photo->getInterpolatedPoints_Tumor(), process_3D::STLType::Tumor);
			 response->setSuccessOrFailure(msclr::interop::marshal_as<System::String^>(errorMessage));
			 if (response->isSuccessful()) {
				 List<System::String^>^ paths = gcnew List<System::String^>();
				 paths->Add(msclr::interop::marshal_as<System::String^>(proc->getTumorGeometry()));
				 response->setData(paths);
			 }
		 }


		void writeThicknessPoints(List<List<EDITCore::CVPoint^>^>^ thicknessPoints) {
			photo->finalizeAllThicknessContours(listPointsToVectorPoints(thicknessPoints));
			releaseMemory(thicknessPoints);
			photo->writeThicknessPoints();
		 }

		//----------------------------------------------FILL VARIABLES WHEN LOAD DATA FROM UI---------------------------------
		void fill2DVariablesWhenLoadDataFromUI(int startingFrame, int endingFrame){ //List<List<EDITCore::CVPoint^>^>^ bladderPoints
			setDicomTags();
			if (response->isSuccessful()) {
				ultr->setInitialFrame(startingFrame);
				ultr->setLastFrame(endingFrame);
				photo->setInitialFrame(startingFrame);
				photo->setLastFrame(endingFrame);
				//photo->setlumenPoints(listPointsToVectorPoints(bladderPoints));
			}
		}

	};
}
#endif