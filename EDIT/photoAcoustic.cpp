#pragma once
#include "photoAcoustic.h"

struct sortclass {
	bool operator() (cv::Point2f pt1, cv::Point2f pt2) {
		return (pt1.y < pt2.y); 
	}
} sortPolar;



photoAcoustic::photoAcoustic() {
}


photoAcoustic::~photoAcoustic() {
}


void photoAcoustic::creatDirectories() {
	//Create all neccessary directories
	this->studyDir = this->mainOutputDirectory; //this->filename;
	_mkdir(studyDir.c_str());

	this->outputOXYImagesDir = studyDir + separator() + "OXY_images";
	_mkdir(outputOXYImagesDir.c_str());

	this->outputDeOXYImagesDir = studyDir + separator() + "deOXY_images";
	_mkdir(outputDeOXYImagesDir.c_str());


	this->outputSegmentedImagesDir = studyDir + separator() + "photoacoustic_segmented_images";
	_mkdir(outputSegmentedImagesDir.c_str());

	this->outputPointsDir = studyDir + separator() + "photoacoustic_points";
	_mkdir(outputPointsDir.c_str());

	/*this->loggertxt = this->studyDir + "/logger.txt";
	if(this->enableLogging) this->logFile.open(this->loggertxt);*/
}


void photoAcoustic::exportOXYImages(string dicomPath) {

	this->dicomPath = dicomPath;
	string nameAndExt = experimental::filesystem::path(dicomPath).filename().generic_string();

	size_t lastindex = nameAndExt.find_last_of(".");
	this->filename = nameAndExt.substr(0, lastindex);
	//create dir for all outputs

	wstring dcm_str(dicomPath.length(), L' ');
	copy(dicomPath.begin(), dicomPath.end(), dcm_str.begin());

	//clear folder of images
	for (const auto& entry : std::experimental::filesystem::directory_iterator(this->outputOXYImagesDir))
		std::experimental::filesystem::remove_all(entry.path());
	//clear memory of images
	vector<Mat>().swap(OXYimages);


	CDicomReader* reader = new CDicomReader();
	this->tags = reader->GetDicomInfo(dcm_str.c_str()); //read dicom tags
	this->OXYimages = reader->dcmimage_Mat(dcm_str.c_str(), this->outputOXYImagesDir, tags[2], tags[3], tags[4], tags[5], CDicomReader::ImageChannel::RED);

	reader->~CDicomReader();
}


void photoAcoustic::exportDeOXYImages(string dicomPath) {

	this->dicomPath = dicomPath;
	string nameAndExt = experimental::filesystem::path(dicomPath).filename().generic_string();

	size_t lastindex = nameAndExt.find_last_of(".");
	this->filename = nameAndExt.substr(0, lastindex);
	//create dir for all outputs
	

	wstring dcm_str(dicomPath.length(), L' ');
	copy(dicomPath.begin(), dicomPath.end(), dcm_str.begin());

	//clear folder of images
	for (const auto& entry : std::experimental::filesystem::directory_iterator(this->outputDeOXYImagesDir))
		std::experimental::filesystem::remove_all(entry.path());
	//clear memory of images
	vector<Mat>().swap(deOXYimages);

	CDicomReader* reader = new CDicomReader();
	this->tags = reader->GetDicomInfo(dcm_str.c_str()); //read dicom tags
	this->deOXYimages = reader->dcmimage_Mat(dcm_str.c_str(), this->outputDeOXYImagesDir, tags[2], tags[3], tags[4], tags[5], CDicomReader::ImageChannel::BLUE);

	reader->~CDicomReader();
}


void photoAcoustic::thicknessExtraction(int frame) {
	if (frame == 0) {
		vector<vector<Point2f>>().swap(this->thicknessPoints);
		vector<vector<Point2f>>().swap(finalThicknessPoints);
		vector<double>().swap(meanThickness);
		vector<Point2f>().swap(contourForFix);
		for (int k = this->initialFrame; k <= this->lastFrame; k++) {
			process(k, totalSequenceOrCorrecton::TOTAL);
		}
	}
	else {
		process(frame, totalSequenceOrCorrecton::CORRECTION);
	}
}


void photoAcoustic::process(int frame, totalSequenceOrCorrecton type) {

	double xspace = tags[0] * 10;
	double yspace = tags[1] * 10;

	Mat OXYImage, edge_OXY, red_PA;
	OXYImage = this->OXYimages[frame];
	bilateralFilter(OXYImage, red_PA, 9, 75, 75);
	//threshold(red_PA, red_PA, 100, 255, THRESH_BINARY);
	Canny(red_PA, edge_OXY, 90, 270, 3);

	vector<Point2f> bladder;// = this->lumenPoints[frame - this->initialFrame]; //--------------------------->??????
	if (type == totalSequenceOrCorrecton::TOTAL) {
		convexHull(this->lumenPoints[frame - this->initialFrame], bladder, true);

		bladder.push_back(bladder[0]);
		bladder = smoothContour(bladder, 100);
	}
	else {
		convexHull(this->contourForFix, bladder, true);
		vector<Point2f>().swap(this->contourForFix);
		bladder.push_back(bladder[0]);
		bladder = smoothContour(bladder, 100);
	}

	Point2f center = findCenterOfContour(bladder);

	//shift bladder points to (0,0)
	transform(bladder.begin(), bladder.end(), bladder.begin(), std::bind2nd(std::minus<Point2f>(), center));

	//convert to polar
	Mat xpts(bladder.size(), 1, CV_32F, &bladder[0].x, 2 * sizeof(float));
	Mat ypts(bladder.size(), 1, CV_32F, &bladder[0].y, 2 * sizeof(float));
	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);
	vector<Point2f> polarXY_outer, outer_mean;

	vector<double> thickness(bladder.size());
	fill(thickness.begin(), thickness.end(), 1.0); //initialize all thickness values with 1.0

	for (int j = 0; j < bladder.size(); j++) {
		double degree = 180 * angle.at<float>(j) / CV_PI;
		if (degree < this->minDegree || degree > this->maxDegree) {
			thickness[j] = 0.0;
		}
		polarXY_outer.push_back(Point2f(magnitude.at<float>(j) + this->pixelDistance, angle.at<float>(j)));
		outer_mean.push_back(Point2f(magnitude.at<float>(j), angle.at<float>(j)));
	}

	Mat mag(polarXY_outer.size(), 1, CV_32F, &polarXY_outer[0].x, 2 * sizeof(float));
	Mat ang(polarXY_outer.size(), 1, CV_32F, &polarXY_outer[0].y, 2 * sizeof(float));
	Mat xnew, ynew;
	polarToCart(mag, ang, xnew, ynew);
	vector<Point> expantedBladder;


	for (int j = 0; j < polarXY_outer.size(); j++) {
		expantedBladder.push_back(Point(xnew.at<float>(j) + center.x, ynew.at<float>(j) + center.y));
		bladder[j].x += center.x;
		bladder[j].y += center.y;
	}

	//-----------------------ALGORITHM----------------------

	vector<Point2f> thicknessPoints(thickness.size());
	fill(thicknessPoints.begin(), thicknessPoints.end(), Point2f(0.0, 0.0));
	vector<double> distances(thickness.size());
	fill(distances.begin(), distances.end(), 0);

	int count = 0;
	double meanDistance = 0;

	for (int i = 0; i < thickness.size(); i++) { 

		if (thickness[i] > 0.0) {
			LineIterator it(edge_OXY, bladder[i], expantedBladder[i], 8, false);
			LineIterator itBil(red_PA, bladder[i], expantedBladder[i], 8, false);


			vector<Point2f> potencialPositions;
			vector<int> potencialPositionsBasedOnIntesityCount;
			vector<int> potencialPositionsBasedOnPixelCount;
			vector<Point2f> coordinates(it.count);


			int intensityCount = 0;
			int countPixel = 0;
			for (int j = 0; j < it.count; j++) 
			{
				Vec3b intensity = Vec3b(*it);
				Vec3b intensityBil = Vec3b(*itBil);
				if ((int)intensity.val[1] > 0) {

					double dist = sqrt(pow(countPixel * xspace, 2) + pow(countPixel * yspace, 2));
					if (dist > this->maxThickness) {
						//cout << "greated than maxThickness" << endl;
						break;
					}
					if (dist > this->minThickness) { //&& countPixel > 15
						potencialPositions.push_back(it.pos());
						potencialPositionsBasedOnIntesityCount.push_back(intensityCount);
						potencialPositionsBasedOnPixelCount.push_back(countPixel);
					}
				}

				if ((int)intensityBil.val[1] < 10) {
					intensityCount = 0;
				}
				intensityCount += (int)intensityBil.val[1];
				it++;
				itBil++;
				countPixel++;
			}
			if (!potencialPositionsBasedOnIntesityCount.empty()) {
				int IndexOfmax_1 = max_element(potencialPositionsBasedOnIntesityCount.begin(), potencialPositionsBasedOnIntesityCount.end()) - potencialPositionsBasedOnIntesityCount.begin();
				Point winner = potencialPositions[IndexOfmax_1];
				thicknessPoints[i] = winner;

				double distance = potencialPositionsBasedOnPixelCount[IndexOfmax_1];//sqrt(pow(winner.x - bladder[i].x, 2) + pow(winner.y - bladder[i].y, 2));
				distances[i] = distance;
				count++;
				//cout << "---->Distance: " << distance << endl;
				meanDistance += distance;

			}

			vector<Point2f>().swap(potencialPositions);
			vector<int>().swap(potencialPositionsBasedOnIntesityCount);
			vector<int>().swap(potencialPositionsBasedOnPixelCount);
			vector<Point2f>().swap(coordinates);

			//line(PA_initial, bladder_updated[i], expantedBladder[i], Scalar(255, 255, 255), 1);
		}


	}

	(count != 0) ? meanDistance /= count : meanDistance = 0;

	int countNew = 0;
	int newMeanDistance = 0;
	for (int j = 0; j < thicknessPoints.size(); j++) {
		if (distances[j] > 0.7 * meanDistance && distances[j] < 1.3 * meanDistance) {

			//cout << "point: " << j << ":  " << distances[j] << endl;
			newMeanDistance += meanDistance;
			countNew++;
		}
		else {
			thicknessPoints[j] = Point2f(0.0, 0.0);
		}
	}

	(countNew != 0) ? newMeanDistance /= countNew : newMeanDistance = 0;

	//---------------------------MEAN TICKNESS-------------------
	for (int i = 0; i < outer_mean.size(); i++) {
		outer_mean[i].x += round(newMeanDistance);
	}
	Mat magT(outer_mean.size(), 1, CV_32F, &outer_mean[0].x, 2 * sizeof(float));
	Mat angT(outer_mean.size(), 1, CV_32F, &outer_mean[0].y, 2 * sizeof(float));
	Mat xnewT, ynewT;

	vector<Point2f> meanThickness;
	polarToCart(magT, angT, xnewT, ynewT);
	for (int j = 0; j < polarXY_outer.size(); j++) {
		meanThickness.push_back(Point2f(xnewT.at<float>(j) + center.x, ynewT.at<float>(j) + center.y));
	}
	//-----------------------------------------------------------

	double meanThick = sqrt(pow(meanDistance * xspace, 2) + pow(meanDistance * yspace, 2));
	
	if (type == totalSequenceOrCorrecton::TOTAL) {
		this->meanThickness.push_back(meanThick);
	}
	else {
		this->meanThickness[frame - this->initialFrame] = meanThick;
	}
	

	for (int j = 0; j < thicknessPoints.size(); j++) {
		if (thicknessPoints[j] == Point2f(0.0, 0.0)) {
			thicknessPoints[j] = meanThickness[j];
		}
	}
	if (type == totalSequenceOrCorrecton::TOTAL) {
		this->thicknessPoints.push_back(smoothContour(thicknessPoints, 0));
	}
	else {
		this->contourForFix = smoothContour(thicknessPoints, 0);
		this->thicknessPoints[frame - this->initialFrame] = this->contourForFix;
	}

	//some memory release
	OXYImage.release();
	red_PA.release();
	edge_OXY.release();
	mag.release();
	ang.release();
	magT.release();
	angT.release();
	vector<Point2f>().swap(polarXY_outer);
	vector<Point2f>().swap(outer_mean);
	vector<Point2f>().swap(bladder);
	vector<Point2f>().swap(thicknessPoints);
	vector<Point2f>().swap(meanThickness);

}


void photoAcoustic::extractOXYandDeOXYPoints(vector<vector<Point2f>> bladderContours, vector<vector<Point2f>> thicknessContours, Point3DType type){

	if (type == Point3DType::OXY) {
		vector<vector<Point3f>>().swap(this->OXYPoints);
	}
	else if (type == Point3DType::DeOXY) {
		vector<vector<Point3f>>().swap(this->DeOXYPoints);
	}


	Point2f imageCenter;

	imageCenter.x = (this->tags[3] - this->tags[2]) / 2; //center_x = (Xmax - Xmin)/2
	imageCenter.y = (this->tags[5] - this->tags[4]) / 2; //center_y = (Ymax - Ymin)/2 
	double xspace = this->tags[0] * 10;
	double yspace = this->tags[1] * 10;
	

	for (int i = 0; i < bladderContours.size(); i++) {
		
		bladderContours[i] = smoothContour(bladderContours[i], 100);
		thicknessContours[i] = smoothContour(thicknessContours[i], 100);

		Mat image, black, result;

		if (type == Point3DType::OXY) {
			image = this->OXYimages[i + this->initialFrame];
		}
		else if (type == Point3DType::DeOXY) {
			image = this->deOXYimages[i + this->initialFrame];
		}
		
		threshold(image, image, 100, 255, THRESH_BINARY);

		black = Mat::zeros(image.size(), image.type());
		
		vector<Point>  rounded_thicknessContoursPoint;
		vector<Point>  rounded_bladderContoursPoint;
		for (int j = 0; j < bladderContours[i].size(); j++) {
			rounded_thicknessContoursPoint.push_back(Point(round(thicknessContours[i][j].x), round(thicknessContours[i][j].y)));
			rounded_bladderContoursPoint.push_back(Point(round(bladderContours[i][j].x), round(bladderContours[i][j].y)));
		}

		fillPoly(black, rounded_thicknessContoursPoint, Scalar(255, 255, 255));
		bitwise_and(black, image, result);

		fillPoly(result, rounded_bladderContoursPoint, Scalar(0, 0, 0));


		Mat Pixels;   // output, locations of non-zero pixels

		findNonZero(result, Pixels);

		vector<Point3f> FramePoints;


		for (int k = 0; k < Pixels.total(); k++) {
			FramePoints.push_back(Point3f((Pixels.at<Point>(k).x - imageCenter.x) * xspace, (Pixels.at<Point>(k).y - imageCenter.y) * yspace, this->distanceBetweenFrames * i));
		}

		if (type == Point3DType::OXY) {
			this->OXYPoints.push_back(FramePoints);
		}
		else if (type == Point3DType::DeOXY) {
			this->DeOXYPoints.push_back(FramePoints);
		}

		image.release();
		black.release();
		result.release();
		Pixels.release();;
		vector<Point3f>().swap(FramePoints);
		vector<Point>().swap(rounded_thicknessContoursPoint);
		vector<Point>().swap(rounded_bladderContoursPoint);
	}
}



void photoAcoustic::extractOXYandDeOXYPoints2(vector<vector<Point2f>> bladderContours, vector<vector<Point2f>> thicknessContours, Point3DType type) {

	vector<vector<Point3f>>().swap(notSharderPoints);
	vector<vector<vector<Point3f>>>().swap(sharderPoints);
	vector<vector<vector<Point3f>>>().swap(interpolatedPoints);

	auto comp = [](const Point3f& lhs, const Point3f& rhs) {return ((lhs.x == rhs.x) && (lhs.y == rhs.y)); };

	for (int i = 0; i < bladderContours.size() - 1; i++) {

		//i-th frame
		Mat imageFrameB;

		if (type == Point3DType::OXY) {
			imageFrameB = this->OXYimages[i + 1 + this->initialFrame];
		}
		else if (type == Point3DType::DeOXY) {
			imageFrameB = this->deOXYimages[i + 1 + this->initialFrame];
		}

		vector<Point3f> Frame1Points3D;
		if (i == 0) {
			Mat imageFrameA;
			if (type == Point3DType::OXY) {
				imageFrameA = this->OXYimages[i + this->initialFrame];
			}
			else if (type == Point3DType::DeOXY) {
				imageFrameA = this->deOXYimages[i + this->initialFrame];
			}
			bladderContours[i] = smoothContour(bladderContours[i], 100);
			thicknessContours[i] = smoothContour(thicknessContours[i], 100);
			Frame1Points3D = findPixelsBetweenThicknessAndBladder(imageFrameA, bladderContours[i], thicknessContours[i], i);
			imageFrameA.release();
		}
		else {
			Frame1Points3D = this->alreadProcessedFramePoints3D;
		}
		bladderContours[i + 1] = smoothContour(bladderContours[i + 1], 100);
		thicknessContours[i + 1] = smoothContour(thicknessContours[i + 1], 100);
		vector<Point3f> Frame2Points3D = findPixelsBetweenThicknessAndBladder(imageFrameB, bladderContours[i+1], thicknessContours[i+1], i+1);
		vector<Point3f>().swap(alreadProcessedFramePoints3D);
		this->alreadProcessedFramePoints3D = Frame2Points3D;

		vector<Point3f> sharedPoints;

		set_union(Frame1Points3D.begin(), Frame1Points3D.end(), Frame2Points3D.begin(), Frame2Points3D.end(), std::back_inserter(sharedPoints), comp);
		
		//------------------------------------------------------------------------------------------
		
		vector<vector<Point3f>> frameSharedPoints;
		vector<vector<Point3f>> frameInterpolationPointsPoints;
		double lower, upper;
		for (int k = 0; k < sharedPoints.size(); k++) {
			vector<Point3f> coupleOfSharedPoints;
			coupleOfSharedPoints.push_back(Point3f(sharedPoints[k].x, sharedPoints[k].y, this->distanceBetweenFrames * i));
			coupleOfSharedPoints.push_back(Point3f(sharedPoints[k].x, sharedPoints[k].y, this->distanceBetweenFrames * (i+1)));
			frameSharedPoints.push_back(coupleOfSharedPoints);
			vector<Point3f>().swap(coupleOfSharedPoints);

			vector<Point3f> interpolationPoints;

			lower = (this->distanceBetweenFrames * i) + 0.05;
			upper = (this->distanceBetweenFrames * (i + 1));

			for (double s = lower; s < upper; s = s + 0.05) {
				interpolationPoints.push_back(Point3f(sharedPoints[k].x, sharedPoints[k].y, s));
			}
			//-----------------------
			frameInterpolationPointsPoints.push_back(interpolationPoints);
			vector<Point3f>().swap(interpolationPoints);
		}

		if (sharedPoints.size() == 0) {
			frameSharedPoints.push_back(vector<Point3f>());
			frameInterpolationPointsPoints.push_back(vector<Point3f>());
		}

		this->notSharderPoints.push_back(Frame1Points3D);
		this->sharderPoints.push_back(frameSharedPoints);
		this->interpolatedPoints.push_back(frameInterpolationPointsPoints);

		imageFrameB.release();
		vector<Point3f>().swap(sharedPoints);
		vector<Point3f>().swap(Frame1Points3D);
		vector<Point3f>().swap(Frame2Points3D);
		vector<vector<Point3f>>().swap(frameSharedPoints);
		vector<vector<Point3f>>().swap(frameInterpolationPointsPoints);
		
	}
}

vector<Point3f> photoAcoustic::findPixelsBetweenThicknessAndBladder(Mat image, vector<Point2f> bladderContours, vector<Point2f> thicknessContours, int iter) {
	
	Mat black, result;
	threshold(image, image, 100, 255, THRESH_BINARY);
	black = Mat::zeros(image.size(), image.type());

	vector<Point>  rounded_thicknessContoursPoint;
	vector<Point>  rounded_bladderContoursPoint;
	for (int j = 0; j < bladderContours.size(); j++) {
		rounded_thicknessContoursPoint.push_back(Point(round(thicknessContours[j].x), round(thicknessContours[j].y)));
		rounded_bladderContoursPoint.push_back(Point(round(bladderContours[j].x), round(bladderContours[j].y)));
	}
	fillPoly(black, rounded_thicknessContoursPoint, Scalar(255, 255, 255));
	bitwise_and(black, image, result);
	fillPoly(result, rounded_bladderContoursPoint, Scalar(0, 0, 0));
	Mat Pixels;   // output, locations of non-zero pixels
	findNonZero(result, Pixels);
	vector<Point3f> Frame1Points3D;
	for (int k = 0; k < Pixels.total(); k++) {
		Frame1Points3D.push_back(Point3f((Pixels.at<Point>(k).x - this->imageCenter.x) * this->xspace, (Pixels.at<Point>(k).y - this->imageCenter.y) * this->yspace, this->distanceBetweenFrames * iter));
	}

	black.release();
	result.release();
	Pixels.release();
	vector<Point>().swap(rounded_thicknessContoursPoint);
	vector<Point>().swap(rounded_bladderContoursPoint);

	return Frame1Points3D;
}


void photoAcoustic::writeThicknessPoints() {
	Mat OXYImage;
	vector<Point2f> thicknessPoints;
	for (int k = this->initialFrame; k <= this->lastFrame; k++) {
		int frame = k - this->initialFrame;


		String Thicknesstxt = this->outputPointsDir + "/" + to_string(k) + ".txt";
		ofstream fileThickness;

		thicknessPoints = this->finalThicknessPoints[k - initialFrame];
		OXYImage = this->OXYimages[k].clone();

		for (int i = 0; i < thicknessPoints.size(); i++) {
			fileThickness << this->thicknessPoints[frame][i].x << " " << this->thicknessPoints[frame][i].y << endl;
			OXYImage.at<uchar>(thicknessPoints[i].y, thicknessPoints[i].x) = 255;
		}

		String pathbmp = this->outputSegmentedImagesDir + "/" + to_string(k) + ".bmp";
		imwrite(pathbmp, OXYImage);
		fileThickness.close();
	}

	OXYImage.release();
	vector<Point2f>().swap(thicknessPoints);
}


Point2f photoAcoustic::findCenterOfContour(vector<Point2f> contour) {
	Point2f center = accumulate(contour.begin(), contour.end(), Point2f(0.0, 0.0));
	center.x /= contour.size();
	center.y /= contour.size();

	return center;
}

void photoAcoustic::finalizeAllThicknessContours (vector<vector<Point2f>> thicknessContours) {
	vector<vector<Point2f>>().swap(finalThicknessPoints);
	for (int i = 0; i < thicknessContours.size(); i++) {
		//this->finalThicknessPoints.push_back(smoothContour(thicknessContours[i], 100));
		this->finalThicknessPoints.push_back(sortUsingPolarCoordinates(thicknessContours[i], 100));
	}	
}


vector<Point2f> photoAcoustic::smoothContour(vector<Point2f> contour, int num_spline) {

	vector<vector<vector<double>>> cl_3D;
	vector<vector<double>> temp_cl_3D;
	vector<double> xx;
	vector<double> yy;
	vector<double> zz;

	for (int i = 0; i < contour.size(); i++) {
		vector<double> temp_temp_cl_3D;
		temp_temp_cl_3D.push_back(contour[i].x);
		temp_temp_cl_3D.push_back(contour[i].y);
		temp_cl_3D.push_back(temp_temp_cl_3D);
	}
	cl_3D.push_back(temp_cl_3D);

	contour = smoothCurve(cl_3D, num_spline);

	//free some memory
	vector<vector<vector<double>>>().swap(cl_3D);
	vector<vector<double>>().swap(temp_cl_3D);
	vector<double>().swap(xx);
	vector<double>().swap(yy);
	vector<double>().swap(zz);

	return contour;
}

vector<Point2f> photoAcoustic::sortUsingPolarCoordinates(vector<Point2f> p, int num_spline) {

	Point2f center = accumulate(p.begin(), p.end(), Point2f(0.0, 0.0));
	center.x /= p.size();
	center.y /= p.size();

	transform(p.begin(), p.end(),
		p.begin(), std::bind2nd(std::minus<Point2f>(), center));

	Mat xpts(p.size(), 1, CV_32F, &p[0].x, 2 * sizeof(float));
	Mat ypts(p.size(), 1, CV_32F, &p[0].y, 2 * sizeof(float));

	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);

	vector<Point2f> polarXY;
	for (int i = 0; i < p.size(); i++) {
		polarXY.push_back(Point2f(magnitude.at<Float32>(i), angle.at<Float32>(i)));
	}

	sort(polarXY.begin(), polarXY.end(), sortPolar);

	Mat mag(polarXY.size(), 1, CV_32F, &polarXY[0].x, 2 * sizeof(float));
	Mat ang(polarXY.size(), 1, CV_32F, &polarXY[0].y, 2 * sizeof(float));

	Mat xnew, ynew, xnewSkin, ynewSkin;
	polarToCart(mag, ang, xnew, ynew);

	vector<Point2f> pp;

	for (int i = 0; i < p.size(); i++) {
		pp.push_back(Point2f(xnew.at<Float32>(i), ynew.at<Float32>(i)));
	}

	pp.push_back(pp[0]);

	vector<Point2f>().swap(p);
	p = smoothContour(pp, num_spline); //num_spline = 50;

	transform(p.begin(), p.end(), p.begin(), std::bind2nd(std::plus<Point2f>(), center));

	xpts.release();
	ypts.release();
	mag.release();
	ang.release();
	xnew.release();
	ynew.release();
	xnewSkin.release();
	ynewSkin.release();
	vector<Point2f>().swap(polarXY);
	vector<Point2f>().swap(pp);
	
	LoggerMessage("Sorting based polar coordinates was performed successfully!");

	return p;
}


// Cubic spline interpolation to smooth centerline
vector<Point2f> photoAcoustic::smoothCurve(vector<vector<vector<double>>> centerline, int num_spline) { //vector<vector<double>>
	vector<vector<vector<double>>> smoothCtr;
	for (int ii = 0; ii < centerline.size(); ++ii) {
		// Read centerline into vtk points
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		for (int i = 0; i < centerline[ii].size(); ++i) {
			points->InsertNextPoint(centerline[ii][i][0], centerline[ii][i][1], centerline[ii][i][1]); // ATTENTION TO X-Y ORDER
		}
		// Create polydata using the set of path points
		vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
		data->SetPoints(points);
		double resamplingParameter = 1; // 10 Sampling resolution
		double distance = 0; // distance from start to end of centerline
							 // Calculate distance
		for (int i = 0; i < data->GetNumberOfPoints() - 1; ++i) {
			double p1[3], p2[3];
			data->GetPoint(i, p1);
			data->GetPoint(i + 1, p2);
			double diff_x = pow(p2[0] - p1[0], 2.0);
			double diff_y = pow(p2[1] - p1[1], 2.0);
			double diff_z = pow(p2[2] - p1[2], 2.0);
			distance += sqrt(diff_x + diff_y + diff_z);
		}
		// Create arrays to store samples for each axis
		vtkSmartPointer<vtkDoubleArray> x, y, z;
		x = vtkSmartPointer<vtkDoubleArray>::New();
		x->SetNumberOfComponents(1);
		y = vtkSmartPointer<vtkDoubleArray>::New();
		y->SetNumberOfComponents(1);
		z = vtkSmartPointer<vtkDoubleArray>::New();
		z->SetNumberOfComponents(1);
		double p0[3];
		// Sample centerline and insert points to respectable arrays
		for (int i = 0; i < data->GetNumberOfPoints(); i = i + resamplingParameter) {
			data->GetPoint(i, p0);
			x->InsertNextValue(p0[0]);
			y->InsertNextValue(p0[1]);
			z->InsertNextValue(p0[2]);
		}
		// Create a spline object for each axis
		vtkSmartPointer<vtkCardinalSpline> xspline = vtkSmartPointer<vtkCardinalSpline>::New();
		xspline->ClosedOff(); // Curve is not closed
		vtkSmartPointer<vtkCardinalSpline> yspline = vtkSmartPointer<vtkCardinalSpline>::New();
		yspline->ClosedOff();
		vtkSmartPointer<vtkCardinalSpline> zspline = vtkSmartPointer<vtkCardinalSpline>::New();
		zspline->ClosedOff();
		double sampledPoints = (double)x->GetNumberOfTuples();
		// Spline interpolation
		// t represents the position of each point on the curve
		// e.g.: 0: start point, 1: end point, 0.5: mid point
		for (int i = 0; i < sampledPoints - 1; ++i) {
			double t = (double)i / sampledPoints;
			double xpoint = x->GetTuple1(i);
			xspline->AddPoint(t, xpoint); // x axis
			double ypoint = y->GetTuple1(i);
			yspline->AddPoint(t, ypoint); // y axis
			double zpoint = z->GetTuple1(i);
			zspline->AddPoint(t, zpoint); // z axis
		}
		// Set final point for interpolation
		double xpoint = x->GetTuple1(sampledPoints - 1);
		xspline->AddPoint(1, xpoint);
		double ypoint = y->GetTuple1(sampledPoints - 1);
		yspline->AddPoint(1, ypoint);
		double zpoint = z->GetTuple1(sampledPoints - 1);
		zspline->AddPoint(1, zpoint);
		double evaluationPoints = 10000, dt = 1 / evaluationPoints; // 10000
		// Create an array to store interpolated points
		double** p = new double* [evaluationPoints];
		for (int i = 0; i < evaluationPoints; ++i) {
			p[i] = new double[3];
		}
		int i = 0;
		double t = 0;
		// Evaluate splines at small intervals and store points
		while (i < evaluationPoints) {
			p[i][0] = xspline->Evaluate(t);
			p[i][1] = yspline->Evaluate(t);
			p[i][2] = zspline->Evaluate(t);
			i++; t = t + dt;
		}
		// Calculate distance between the spline points
		double smoothed_distance = 0;
		double* newdst = new double[evaluationPoints - 1];
		for (int i = 0; i < evaluationPoints - 1; ++i) {
			double diff_x = pow(p[i + 1][0] - p[i][0], 2.0);
			double diff_y = pow(p[i + 1][1] - p[i][1], 2.0);
			double diff_z = pow(p[i + 1][2] - p[i][2], 2.0);
			smoothed_distance += sqrt(diff_x + diff_y + diff_z);
			newdst[i] = sqrt(diff_x + diff_y + diff_z);
		}
		//cout << smoothed_distance << " ";
		double planeDistance;
		if (num_spline == 0) {
			planeDistance = 13; // sampling distance 54
		}
		else {
			planeDistance = ceil(smoothed_distance) / num_spline; //ceil(smoothed_distance) / num_spline
		}
		//double planeDistance = ceil(smoothed_distance) / num_spline;//0.5; // sampling distance 54
		//double planeDistance = ceil(dims[0] * 15 / 1024);

		double interval = 0;
		vector<int> ind; // vector to store indexes of equidistant points
		for (int i = 0; i < evaluationPoints - 1; ++i) {
			interval += newdst[i];
			if (interval >= planeDistance - 0.01) { // select points with equal distance on curve
				interval = 0;
				ind.push_back(i);
			}
		}
		// Store equidistant points
		vector<vector<double>> equidistantPoints;
		for (int i = 0; i < ind.size(); ++i) {
			vector<double> point;
			point.push_back(p[ind[i]][0]); point.push_back(p[ind[i]][1]); point.push_back(p[ind[i]][2]);
			equidistantPoints.push_back(point);
		}
		smoothCtr.push_back(equidistantPoints); // store interpolated centerline

		vector<vector<double>>().swap(equidistantPoints);
		vector<int>().swap(ind);
	}


	//convert point from double to Point (opencv)
	vector<Point2f> t;

	t.push_back(Point(centerline[0][0][0], centerline[0][0][1])); //push the first point into the vector

	for (int i = 0; i < smoothCtr[0].size(); i++) {
		t.push_back(Point2f(smoothCtr[0][i][0], smoothCtr[0][i][1]));
	}

	//free memory
	vector<vector<vector<double>>>().swap(smoothCtr);

	return t;// smoothCtr[0]; //return only the first centerline 
}