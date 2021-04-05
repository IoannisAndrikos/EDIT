#pragma once

#include "ultrasound.h"

#include "morphsnakes.h"

#include "messages.h"

struct sortBasedYCoordinate {
	bool operator() (cv::Point2f pt1, cv::Point2f pt2) { return (pt1.y  < pt2.y ); }
} sortBasedYCoordinate;


struct EuclideanDistance
{
	EuclideanDistance(const Point2f& _p) : p(_p) {}

	bool operator()(const Point2f& lhs, const Point2f& rhs) const
	{
		return sqrt(pow(lhs.x - p.x, 2) + pow(lhs.y - p.y, 2)) < sqrt(pow(rhs.x - p.x, 2) + pow(rhs.y - p.y, 2));
	}

private:
	Point2f p;
};

struct linedToCenter
{
	linedToCenter(const Point2f& _p) : p(_p) {}

	bool operator()(const Point2f& lhs) const
	{
		return lhs.x == p.x;
	}

private:
	Point2f p;
};


bool Longest(const std::vector<Point>& lhs, const std::vector<Point>& rhs)
{
	return lhs.size() < rhs.size();
}

ultrasound::ultrasound() {}
ultrasound::~ultrasound() {}

string ultrasound::exportImages(string dicomPath) {

	if (experimental::filesystem::path(dicomPath).filename().extension().generic_string()!= ".dcm") {
		return warningMessages::cannotReadTheDicomFile;
	}
	
	this->dicomPath = dicomPath;
	string nameAndExt = experimental::filesystem::path(dicomPath).filename().generic_string();

	size_t lastindex = nameAndExt.find_last_of(".");
	this->filename = nameAndExt.substr(0, lastindex);
	//create dir for all outputs

	wstring dcm_str(dicomPath.length(), L' ');
	copy(dicomPath.begin(), dicomPath.end(), dcm_str.begin());

	//clear folder of images
	for (const auto& entry : std::experimental::filesystem::directory_iterator(this->outputImagesDir))
		std::experimental::filesystem::remove_all(entry.path());
	
	//clear memory of images
	vector<Mat>().swap(images);

	CDicomReader* reader = new CDicomReader();
	
	this->tags = reader->GetDicomInfo(dcm_str.c_str()); //read dicom tags
	if (this->tags.size() < 6) {
		return warningMessages::cannotGetAllNecessaryDicomTags;
	}
		
	this->images = reader->dcmimage_Mat(dcm_str.c_str(), this->outputImagesDir, tags[2], tags[3], tags[4], tags[5]);
	//reader->dcmimage_bmp(dcm_str.c_str(), out_str.c_str()); // read the given dicom and save images into the given folder
	if (this->images.size() == 0) {
		return warningMessages::cannotReadTheDicomFile;
	}
	
	reader->~CDicomReader();

	return this->success;
}


void ultrasound::creatDirectories() {
	//Create all neccessary directories
	this->studyDir = this->mainOutputDirectory; //this->filename;
	_mkdir(studyDir.c_str());

	this->outputImagesDir = studyDir + separator() + "ultrasound_images";
	_mkdir(outputImagesDir.c_str());

	this->outputSegmentedImagesDir = studyDir + separator() +  "ultrasound_segmented_images";
	_mkdir(outputSegmentedImagesDir.c_str());

	this->outputPointsDir = studyDir + separator() + "ultrasound_points";
	_mkdir(outputPointsDir.c_str());

	/*this->loggertxt = this->studyDir + "/logger.txt";
	if(this->enableLogging) this->logFile.open(this->loggertxt);*/
}



template<class T> morphsnakes::NDImage<T, 2> cimg2ndimage(CImg<T>& img)
{
	morphsnakes::Shape<2> shape = { img.height(), img.width() };
	morphsnakes::Stride<2> stride = { img.width() * sizeof(T), sizeof(T) };

	return morphsnakes::NDImage<T, 2>(img.data(), shape, stride);
}

string ultrasound::processing(int initialFrame, int lastFrame, Point clickPoint) {

	vector<vector<Point2f>>().swap(this->lumenPoints);
	vector<vector<Point2f>>().swap(this->skinPoints);
	vector<double>().swap(this->lumenArea);

	this->initialFrame = initialFrame;
	this->lastFrame = lastFrame;

	for (int i = initialFrame; i <= lastFrame; i++) { //initialFrame - 1
		LoggerMessage(string("Working on Frame: " + to_string(i)));
		Mat lumen;
		if (applyEqualizeHist) {
			equalizeHist(this->images[i], lumen);
		}
		else {
			lumen = this->images[i].clone();
		}

		CImg<double>  imgC = *cvImgToCImg(lumen);
		
		//Morphological ACWE
		if (i == this->initialFrame || i >=this->lastFrame-5) { // maybe 5 has to be configurable from gui
			//in the first as well as in the final 5 frames we apply circle levelset image and 200 repeats
			CImg<unsigned char> embedding = circle_levelset(imgC.height(), imgC.width(), { clickPoint.y , clickPoint.x }, this->levelsetSize);
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
			for (int i = 0; i < 200; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((embedding * 255));
			this->levelsetImage = embedding;
		}
		else {
			//while in the rest images we use the segmentation of the previous image and 50 repeats
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(this->levelsetImage), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
			for (int i = 0; i < this->repeats; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((this->levelsetImage * 255));
		}

		//Canny(lumen, lumen, 0, 0, 3);
	   //threshold(lumen, lumen, 100, 255, THRESH_BINARY);

		//---------------------------------------------Keep the longest contour--------------------------------------
		threshold(lumen, lumen, 100, 255, THRESH_BINARY);
		vector<vector<Point>> contours;
		vector<Vec4i> hierarchy;
		Canny(lumen, lumen, 0, 0, 3);//with or without, explained later.
		findContours(lumen, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		Mat drawing = Mat::zeros(lumen.size(), lumen.type());

		Scalar color(255, 0, 0);
		
		//Make sure that the longest contour is the 0
		drawContours(drawing, contours, findLongestVector(contours), color);

		lumen = drawing;

		//find white pixels and calculte their center and store points into vector!
		vector<Point2f>* lumenPoints = new vector<Point2f>();
		Point2f* lumenCenter = new Point2f();
		Point* highestWhitePixel = new Point();

		try {
			ResultOfProcess poitsWereFound = centerAndPointsOfContour(lumen, lumenPoints, lumenCenter, highestWhitePixel);
			//update the click point! This point is used for the segmentation of the final (almost 5) frames
			clickPoint = { (int)lumenCenter->x, (int)lumenCenter->y };


			if (poitsWereFound == ResultOfProcess::SUCCESS) {
				ResultOfProcess poitsWereSorted = sortClockwise(lumenPoints, lumenCenter, i);
				if (poitsWereSorted == ResultOfProcess::FAILURE) {
					return warningMessages::firtFrameSegmentationFailed + to_string(i);
				}
			}
			else if (poitsWereFound == ResultOfProcess::FAILURE && i == this->initialFrame) {
				return warningMessages::firtFrameSegmentationFailed;
			}
		}
		catch (exception e) {
			return warningMessages::frameSegmentationFailed + to_string(i);
		}
		
		vector<Point2f>().swap(*lumenPoints);
		vector<vector<Point>>().swap(contours);
		vector<Vec4i>().swap(hierarchy);
		drawing.release();
		lumen.release();
	}
	return this->success;
}

string ultrasound::recalculate(int frame, vector<Point2f> points) {

	LoggerMessage(string("Working (Recalculate) on Frame: " + to_string(frame)));
	Mat lumen;
	if (applyEqualizeHist) {
		equalizeHist(this->images[frame], lumen);
	}
	else {
		lumen = this->images[frame].clone();
	}

	Point2f contourClickPoint = getCenterOfMass(points);

	CImg<double>  imgC = *cvImgToCImg(lumen);

	//in the first as well as in the final 5 frames we apply circle levelset image and 200 repeats
	CImg<unsigned char> embedding = circle_levelset(imgC.height(), imgC.width(), { (int)contourClickPoint.y , (int)contourClickPoint.x }, this->levelsetSize);
	morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
	for (int i = 0; i < 200; ++i) {
		macwe.step();
	}
	//Back to Mat
	lumen = CImgtoCvImg((embedding * 255));


	//Canny(lumen, lumen, 0, 0, 3);
   //threshold(lumen, lumen, 100, 255, THRESH_BINARY);

	//---------------------------------------------Keep the longest contour--------------------------------------
	threshold(lumen, lumen, 100, 255, THRESH_BINARY);
	vector<vector<Point>> contours;
	vector<Vec4i> hierarchy;
	Canny(lumen, lumen, 0, 0, 3);//with or without, explained later.
	findContours(lumen, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

	Mat drawing = Mat::zeros(lumen.size(), lumen.type());

	Scalar color(255, 0, 0);

	//Make sure that the longest contour is the 0
	drawContours(drawing, contours, findLongestVector(contours), color);

	lumen = drawing;

	//find white pixels and calculte their center and store points into vector!
	vector<Point2f>* lumenPoints = new vector<Point2f>();
	Point2f* lumenCenter = new Point2f();
	Point* highestWhitePixel = new Point();

	try {
		ResultOfProcess poitsWereFound = centerAndPointsOfContour(lumen, lumenPoints, lumenCenter, highestWhitePixel);

		if (poitsWereFound == ResultOfProcess::SUCCESS) {
			ResultOfProcess poitsWereSorted = sortClockwise(lumenPoints, lumenCenter, frame, true);
			if (poitsWereSorted == ResultOfProcess::FAILURE) {
				return warningMessages::firtFrameSegmentationFailed;
			}
		}
		else {
			return warningMessages::firtFrameSegmentationFailed;
		}
	}
	catch (exception e) {
		return warningMessages::frameSegmentationFailed + to_string(frame);
	}

	vector<Point2f>().swap(*lumenPoints);
	vector<vector<Point>>().swap(contours);
	vector<Vec4i>().swap(hierarchy);
	drawing.release();
	lumen.release();
	
	return this->success;
}



//-------------------------------------------------------------------------------------------------------


string ultrasound::extactTumor2D(Point clickPoint, vector<vector<Point2f>> lumen2DPoints, vector<vector<Point2f>> thickness2DPoints) {

	vector<Mat>().swap(this->tumorImages);

	CImg<unsigned char> embedding;

	for (int i = 0; i < thickness2DPoints.size(); i++) { //initialFrame - 1

		String imagePath = String("C:/Users/Legion Y540/desktop/remove_tumor/" + to_string(i + this->initialFrame)) + ".bmp";

		Mat lumen = imread(imagePath, 0);

		CImg<double>  imgC = *cvImgToCImg(lumen);

		//Morphological ACWE
		if (i == 0) { // maybe 5 has to be configurable from gui
			//in the first as well as in the final 5 frames we apply circle levelset image and 200 repeats
			CImg<unsigned char> embedding = circle_levelset(imgC.height(), imgC.width(), { clickPoint.y , clickPoint.x }, this->levelsetSize);
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), 5, 5, 5);
			for (int i = 0; i < 200; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((embedding * 255));
			this->levelsetImage = embedding;
		}
		else {
			//while in the rest images we use the segmentation of the previous image and 50 repeats
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(this->levelsetImage), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
			for (int i = 0; i < this->repeats; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((this->levelsetImage * 255));
		}

		vector<Point>  rounded_bladderContoursPoint, rounded_thicknessContoursPoint;
		for (int j = 0; j < lumen2DPoints[i].size(); j++) {
			rounded_bladderContoursPoint.push_back(Point(round(lumen2DPoints[i][j].x), round(lumen2DPoints[i][j].y)));
		}
		fillPoly(lumen, rounded_bladderContoursPoint, Scalar(0, 0, 0));
		Mat element = getStructuringElement(MORPH_RECT, Size(15, 15), Point(-1, -1)); // kernel performing drode 
		erode(lumen, lumen, element);

		for (int j = 0; j < thickness2DPoints[i].size(); j++) {
			rounded_thicknessContoursPoint.push_back(Point(round(thickness2DPoints[i][j].x), round(thickness2DPoints[i][j].y)));
		}

		Mat black = Mat::zeros(lumen.size(), lumen.type());
		fillPoly(black, rounded_thicknessContoursPoint, Scalar(255, 255, 255));
		bitwise_and(black, lumen, lumen);

		/*namedWindow("window", 1);
		imshow("window", lumen);
		waitKey(0);*/

		this->tumorImages.push_back(lumen);

		lumen.release();
	}
	return this->success;
}

vector<vector<Point2f>> ultrasound::getTumorBorders() {
	vector<vector<Point2f>> tumorBorders;

	vector<vector<Point>> contours;
	vector<Vec4i> hierarchy;
	int index, minDistIndex;
	vector<Point2f> tumorContour;
	double dist;
	Mat black, whitePixels;
	Scalar color(255, 0, 0);
	for (int i = 0; i < this->tumorImages.size(); i++) {
		vector<Point2f>().swap(tumorContour);

		black = Mat::zeros(this->tumorImages[0].size(), this->tumorImages[0].type());

		findContours(this->tumorImages[i], contours, hierarchy, 3, 2, Point(0, 0));
		if (!contours.empty()) {
			index = findLongestVector(contours);

			//Canny(black, black, 0, 0, 3);//with or without, explained later.
			drawContours(black, contours, index, color);
			/*namedWindow("window", 1);
			imshow("window", black);
			waitKey(0); */
		
			////TO DO
			//// output, locations of non-zero pixels

			findNonZero(black, whitePixels);

			for (int i = 0; i < whitePixels.total(); i++) {
				Point lumen_p = Point(whitePixels.at<Point>(i).x, whitePixels.at<Point>(i).y);
				tumorContour.push_back(Point2f(whitePixels.at<Point>(i).x, whitePixels.at<Point>(i).y));
			}

			//transform(contours[index].begin(), contours[index].end(), back_inserter(tumorContour), [](const Point& p) { return (Point2f)p; });
			tumorBorders.push_back(smoothContour(sortBasedEuclideanDistance(tumorContour),0,true));
		}
		else {
			tumorBorders.push_back(vector<Point2f>()); //push back empty
		}
	}
	return tumorBorders;
}

//for all sequence
string ultrasound::fixArtifact(Point clickPoint, vector<vector<Point2f>> points) {

	vector<vector<Point2f>>().swap(this->lumenPoints);
	vector<vector<Point2f>>().swap(this->skinPoints);
	vector<double>().swap(this->lumenArea);

	this->initialFrame = initialFrame;
	this->lastFrame = lastFrame;

	for (int i = initialFrame; i <= lastFrame; i++) { //initialFrame - 1
		LoggerMessage(string("Working on Frame: " + to_string(i)));
		Mat lumen;
		if (applyEqualizeHist) {
			equalizeHist(this->images[i], lumen);
		}
		else {
			lumen = this->images[i].clone();
		}

		//TO DO 
		//-------------------------------------------
		vector<Point2f> bladder, convexBladder;
		bladder = points[i - this->initialFrame];
		Point2f center = getCenterOfGravity(bladder);
		convexHull(bladder, convexBladder, true, false);
		//convexBladder = interpolateConvexPoints(convexBladder, interpolationMethod::AKIMA);
		convexBladder = smoothContour(convexBladder, 200, true);
		Mat black = Mat::zeros(lumen.size(), lumen.type());
		vector<Point>  rounded_bladder, rounded_convexBladder;
		for (int j = 0; j < bladder.size(); j++) {
			rounded_bladder.push_back(Point(round(bladder[j].x), round(bladder[j].y)));
			rounded_convexBladder.push_back(Point(round(convexBladder[j].x), round(convexBladder[j].y)));
		}
		fillPoly(black, rounded_convexBladder, Scalar(255, 0, 0));
		fillPoly(black, rounded_bladder, Scalar(0, 0, 0));


		vector<vector<Point>> artifact_Contours;
		vector<Vec4i> artifact_hierarchy;
		findContours(black, artifact_Contours, artifact_hierarchy, 3, 2, Point(0, 0));


		vector<Point2f> centers;

		for (int i = 0; i < artifact_Contours.size(); i++) {
			centers.push_back(getCenterOfGravity(artifact_Contours[i]));
		}

		transform(centers.begin(), centers.end(), centers.begin(), std::bind2nd(std::minus<Point2f>(), center));
		Mat magnitude, angle;
		Mat xpts(centers.size(), 1, CV_32F, &centers.at(0).x, 2 * sizeof(float));
		Mat ypts(centers.size(), 1, CV_32F, &centers.at(0).y, 2 * sizeof(float));

		double degrees;

		fillPoly(lumen, rounded_bladder, Scalar(0, 0, 0));

		cartToPolar(xpts, ypts, magnitude, angle);
		vector<Point2f> polarContourCenters;
		for (int i = 0; i < centers.size(); i++) {
			polarContourCenters.push_back(Point2f(magnitude.at<float>(i), angle.at<float>(i)));
			degrees = 180 * polarContourCenters[i].y / CV_PI;

			if (degrees > 220 && degrees < 320) {
				fillPoly(lumen, artifact_Contours[i], Scalar(0, 0, 0));
			}
		}

		vector<Point2f>().swap(bladder);
		vector<Point2f>().swap(convexBladder);
		vector<vector<Point>>().swap(artifact_Contours);
		vector<Vec4i>().swap(artifact_hierarchy);
		vector<Point>().swap(rounded_bladder);
		vector<Point>().swap(rounded_convexBladder);
		black.release();
		magnitude.release();
		angle.release();
		xpts.release();
		ypts.release();
		//-------------------------------------------
		

		CImg<double>  imgC = *cvImgToCImg(lumen);


		//Morphological ACWE
		if (i == this->initialFrame || i >= this->lastFrame - 5) { // maybe 5 has to be configurable from gui
			//in the first as well as in the final 5 frames we apply circle levelset image and 200 repeats
			CImg<unsigned char> embedding = circle_levelset(imgC.height(), imgC.width(), { clickPoint.y , clickPoint.x }, this->levelsetSize);
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
			for (int i = 0; i < 200; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((embedding * 255));
			this->levelsetImage = embedding;
		}
		else {
			//while in the rest images we use the segmentation of the previous image and 50 repeats
			morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(this->levelsetImage), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
			for (int i = 0; i < this->repeats; ++i) {
				macwe.step();
			}
			//Back to Mat
			lumen = CImgtoCvImg((this->levelsetImage * 255));
		}

		//Canny(lumen, lumen, 0, 0, 3);
	   //threshold(lumen, lumen, 100, 255, THRESH_BINARY);

		//---------------------------------------------Keep the longest contour--------------------------------------
		threshold(lumen, lumen, 100, 255, THRESH_BINARY);
		vector<vector<Point>> contours;
		vector<Vec4i> hierarchy;
		Canny(lumen, lumen, 0, 0, 3);//with or without, explained later.
		findContours(lumen, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		Mat drawing = Mat::zeros(lumen.size(), lumen.type());

		Scalar color(255, 0, 0);

		//Make sure that the longest contour is the 0
		drawContours(drawing, contours, findLongestVector(contours), color);

		lumen = drawing;

		//find white pixels and calculte their center and store points into vector!
		vector<Point2f>* lumenPoints = new vector<Point2f>();
		Point2f* lumenCenter = new Point2f();
		Point* highestWhitePixel = new Point();

		try {
			ResultOfProcess poitsWereFound = centerAndPointsOfContour(lumen, lumenPoints, lumenCenter, highestWhitePixel);
			//update the click point! This point is used for the segmentation of the final (almost 5) frames
			clickPoint = { (int)lumenCenter->x, (int)lumenCenter->y };


			if (poitsWereFound == ResultOfProcess::SUCCESS) {
				ResultOfProcess poitsWereSorted = sortClockwise(lumenPoints, lumenCenter, i);
				if (poitsWereSorted == ResultOfProcess::FAILURE) {
					return warningMessages::firtFrameSegmentationFailed + to_string(i);
				}
			}
			else if (poitsWereFound == ResultOfProcess::FAILURE && i == this->initialFrame) {
				return warningMessages::firtFrameSegmentationFailed;
			}
		}
		catch (exception e) {
			return warningMessages::frameSegmentationFailed + to_string(i);
		}

		vector<Point2f>().swap(*lumenPoints);
		vector<vector<Point>>().swap(contours);
		vector<Vec4i>().swap(hierarchy);
		drawing.release();
		lumen.release();
	}
	return this->success;
}

//for one frame
string ultrasound::fixArtifact(int frame, vector<Point2f> points) {

	LoggerMessage(string("Working (Recalculate) on Frame: " + to_string(frame)));
	Mat lumen;
	if (applyEqualizeHist) {
		equalizeHist(this->images[frame], lumen);
	}
	else {
		lumen = this->images[frame].clone();
	}

	//TO DO 
	//-------------------------------------------
	vector<Point2f> bladder, convexBladder;
	bladder = points;
	Point2f center = getCenterOfGravity(bladder);
	convexHull(bladder, convexBladder, true, false);
	//convexBladder = interpolateConvexPoints(convexBladder, interpolationMethod::AKIMA);
	convexBladder = smoothContour(convexBladder, 200, true);
	Mat black = Mat::zeros(lumen.size(), lumen.type());
	vector<Point>  rounded_bladder, rounded_convexBladder;
	for (int j = 0; j < bladder.size(); j++) {
		rounded_bladder.push_back(Point(round(bladder[j].x), round(bladder[j].y)));
		rounded_convexBladder.push_back(Point(round(convexBladder[j].x), round(convexBladder[j].y)));
	}
	fillPoly(black, rounded_convexBladder, Scalar(255, 0, 0));
	fillPoly(black, rounded_bladder, Scalar(0, 0, 0));


	vector<vector<Point>> artifact_Contours;
	vector<Vec4i> artifact_hierarchy;
	findContours(black, artifact_Contours, artifact_hierarchy, 3, 2, Point(0, 0));


	vector<Point2f> centers;

	for (int i = 0; i < artifact_Contours.size(); i++) {
		centers.push_back(getCenterOfGravity(artifact_Contours[i]));
	}

	transform(centers.begin(), centers.end(), centers.begin(), std::bind2nd(std::minus<Point2f>(), center));
	Mat magnitude, angle;
	Mat xpts(centers.size(), 1, CV_32F, &centers.at(0).x, 2 * sizeof(float));
	Mat ypts(centers.size(), 1, CV_32F, &centers.at(0).y, 2 * sizeof(float));

	double degrees;

	fillPoly(lumen, rounded_bladder, Scalar(0, 0, 0));

	cartToPolar(xpts, ypts, magnitude, angle);
	vector<Point2f> polarContourCenters;
	for (int i = 0; i < centers.size(); i++) {
		polarContourCenters.push_back(Point2f(magnitude.at<float>(i), angle.at<float>(i)));
		degrees = 180 * polarContourCenters[i].y / CV_PI;

		if (degrees > 220 && degrees < 320) {
			fillPoly(lumen, artifact_Contours[i], Scalar(0, 0, 0));
		}
	}

	vector<Point2f>().swap(bladder);
	vector<Point2f>().swap(convexBladder);
	vector<vector<Point>>().swap(artifact_Contours);
	vector<Vec4i>().swap(artifact_hierarchy);
	vector<Point>().swap(rounded_bladder);
	vector<Point>().swap(rounded_convexBladder);
	black.release();
	magnitude.release();
	angle.release();
	xpts.release();
	ypts.release();
	//-------------------------------------------


	CImg<double>  imgC = *cvImgToCImg(lumen);

	Point2f clickPointFrame = getCenterOfMass(points);

	//Morphological ACWE
	//in the first as well as in the final 5 frames we apply circle levelset image and 200 repeats
	CImg<unsigned char> embedding = circle_levelset(imgC.height(), imgC.width(), { (int)clickPointFrame.y , (int)clickPointFrame.x }, this->levelsetSize);
	morphsnakes::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), this->smoothing, this->lamda1, this->lamda2);
	for (int i = 0; i < 200; ++i) {
		macwe.step();
	}
	//Back to Mat
	lumen = CImgtoCvImg((embedding * 255));
			
	//Canny(lumen, lumen, 0, 0, 3);
    //threshold(lumen, lumen, 100, 255, THRESH_BINARY);

	//---------------------------------------------Keep the longest contour--------------------------------------
	threshold(lumen, lumen, 100, 255, THRESH_BINARY);
	vector<vector<Point>> contours;
	vector<Vec4i> hierarchy;
	Canny(lumen, lumen, 0, 0, 3);//with or without, explained later.
	findContours(lumen, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

	Mat drawing = Mat::zeros(lumen.size(), lumen.type());

	Scalar color(255, 0, 0);

	//Make sure that the longest contour is the 0
	drawContours(drawing, contours, findLongestVector(contours), color);

	lumen = drawing;

	//find white pixels and calculte their center and store points into vector!
	vector<Point2f>* lumenPoints = new vector<Point2f>();
	Point2f* lumenCenter = new Point2f();
	Point* highestWhitePixel = new Point();

	try {
		ResultOfProcess poitsWereFound = centerAndPointsOfContour(lumen, lumenPoints, lumenCenter, highestWhitePixel);

		if (poitsWereFound == ResultOfProcess::SUCCESS) {
			ResultOfProcess poitsWereSorted = sortClockwise(lumenPoints, lumenCenter, frame, true);
			if (poitsWereSorted == ResultOfProcess::FAILURE) {
				return warningMessages::frameSegmentationFailed;
			}
		}
	}
	catch (exception e) {
		return warningMessages::frameSegmentationFailed;
	}

	vector<Point2f>().swap(*lumenPoints);
	vector<vector<Point>>().swap(contours);
	vector<Vec4i>().swap(hierarchy);
	drawing.release();
	lumen.release();
	return this->success;
}




//---------------------------------------------------------------------------------------



void ultrasound::finalizeAllBladderContours(vector<vector<Point2f>> points) {
	vector<vector<Point2f>>().swap(lumenPoints);
	for (int i = 0; i < points.size(); i++) {
		if (!IsClockwise(points[i])) {
			reverse(points[i].begin(), points[i].end());
		}

		points[i] = smoothContour(points[i], 100, true);

		Point2f center = getCenterOfMass(points[i]);

		transform(points[i].begin(), points[i].end(), points[i].begin(), std::bind2nd(std::minus<Point2f>(), center));

		Mat xpts(points[i].size(), 1, CV_32F, &points[i][0].x, 2 * sizeof(float));
		Mat ypts(points[i].size(), 1, CV_32F, &points[i][0].y, 2 * sizeof(float));
		Mat magnitude, angle;
		cartToPolar(xpts, ypts, magnitude, angle);
		vector<Point2f> polarXY, sortedPolarXY;
		vector<double> polarAngles;
		for (int j = 0; j < points[i].size(); j++) {
			polarXY.push_back(Point2f(magnitude.at<Float32>(j), angle.at<Float32>(j)));
			polarAngles.push_back(angle.at<Float32>(j));
		}
		auto min = std::min_element(polarAngles.begin(), polarAngles.end());

		int minIndex = distance(polarAngles.begin(), min);

		if (minIndex != 0) {
			sortedPolarXY.insert(sortedPolarXY.end(), polarXY.begin() + minIndex, polarXY.end());
			polarXY.erase(polarXY.begin() + minIndex, polarXY.end());
		}
		sortedPolarXY.insert(sortedPolarXY.end(), polarXY.begin(), polarXY.end());

		Mat mag(sortedPolarXY.size(), 1, CV_32F, &sortedPolarXY[0].x, 2 * sizeof(float));
		Mat ang(sortedPolarXY.size(), 1, CV_32F, &sortedPolarXY[0].y, 2 * sizeof(float));
		Mat xnew, ynew;
		polarToCart(mag, ang, xnew, ynew);
		vector<Point2f> pp;

		for (int j = 0; j < points[i].size(); j++) {
			pp.push_back(Point2f(xnew.at<Float32>(j), ynew.at<Float32>(j)));
		}

		transform(pp.begin(), pp.end(), pp.begin(), std::bind2nd(std::plus<Point2f>(), center));
		this->lumenPoints.push_back(pp);

		xpts.release();
		ypts.release();
		mag.release();
		ang.release();
		xnew.release();
		ynew.release();
		vector<double>().swap(polarAngles);
		vector<Point2f>().swap(polarXY);
		vector<Point2f>().swap(sortedPolarXY);
		vector<Point2f>().swap(pp);
			
	}
	vector<vector<Point2f>>().swap(points);
}

void ultrasound::extractSkinPoints(vector<vector<Point2f>> bladderPoints) {

	vector<vector<Point2f>>().swap(skinPoints);
	finalizeAllBladderContours(bladderPoints);

	this->skinPoints = this->getlumenPoints();
	int imageCount = this->initialFrame;

	for (int i = 0; i< skinPoints.size(); i++) {
		 
		Point2f center = getCenterOfMass(skinPoints[i]);

		transform(skinPoints[i].begin(), skinPoints[i].end(), skinPoints[i].begin(), std::bind2nd(std::minus<Point2f>(), center));

		Mat xpts(skinPoints[i].size(), 1, CV_32F, &skinPoints[i][0].x, 2 * sizeof(float));
		Mat ypts(skinPoints[i].size(), 1, CV_32F, &skinPoints[i][0].y, 2 * sizeof(float));
		Mat magnitude, angle;
		cartToPolar(xpts, ypts, magnitude, angle);

		vector<Point2f> SkinPolarXY;
		int index = 0;
		for (int j = 0; j < skinPoints[i].size(); j++) {
			if (abs((angle.at<Float32>(j) * 180 / CV_PI) - 270) < 5) {
				break;
			}
			index++;
		}

		for (int j = 0; j < skinPoints[i].size(); j++) {
			SkinPolarXY.push_back(Point2f(magnitude.at<Float32>(j), angle.at<Float32>(j)));
		}

	Mat image = this->images[imageCount].clone();
	medianBlur(image, image, 9);
	GaussianBlur(image, image, Size(9, 9), 0);
	threshold(image, image, 100, 255, THRESH_BINARY);
	
	int distBladderSkin = 0;

	for (int j = 0; j < skinPoints[i][index].y + center.y; j++) {
		if ((int)image.at<uchar>(j, skinPoints[i][index].x + center.x) == 0) {
			distBladderSkin++;
		}
		else {
			break; //find the first white (255) pixel
		}
	}

	int dist = skinPoints[i][index].y + center.y - distBladderSkin;


	for (int j = 0; j < SkinPolarXY.size(); j++) {
		SkinPolarXY[j].x += dist;
	}

	Mat mag(SkinPolarXY.size(), 1, CV_32F, &SkinPolarXY[0].x, 2 * sizeof(float));
	Mat ang(SkinPolarXY.size(), 1, CV_32F, &SkinPolarXY[0].y, 2 * sizeof(float));
	Mat xnew, ynew;
	polarToCart(mag, ang, xnew, ynew);
	vector<Point2f> pp;


	for (int j = 0; j < SkinPolarXY.size(); j++) {
		pp.push_back(Point2f(xnew.at<Float32>(j) + center.x, ynew.at<Float32>(j) + center.y));
	}

	vector<Point2f>().swap(this->skinPoints[i]);

	this->skinPoints[i] = pp;


	imageCount++;

	//free memory
	image.release();
	xpts.release();
	ypts.release();
	mag.release();
	ang.release();
	xnew.release();
	ynew.release();
	vector<Point2f>().swap(SkinPolarXY);
	vector<Point2f>().swap(pp);
	}
}


void ultrasound::writePointsAndImages() {
	Mat img;
	for (int i = 0; i <= (this->lastFrame - this->initialFrame); i++) {
		String lumenhtxt = this->outputPointsDir + "/" + to_string(i + this->initialFrame) + ".txt";
		ofstream filelumen;
		filelumen.open(lumenhtxt);
		img = this->images[i + this->initialFrame];
		for (int j = 0; j < this->lumenPoints[i].size(); j++) { //p.size()
				filelumen << this->lumenPoints[i][j].x << " " << this->lumenPoints[i][j].y << endl;
				img.at<uchar>(round(this->lumenPoints[i][j].y), round(this->lumenPoints[i][j].x)) = 255;
		}

		for (int j = 0; j < this->skinPoints[i].size(); j++) { //p.size()
			img.at<uchar>(round(this->skinPoints[i][j].y), round(this->skinPoints[i][j].x)) = 255;
		}

		String pathbmp = this->outputSegmentedImagesDir + "/" + to_string(i+ this->initialFrame) + ".bmp";
		imwrite(pathbmp, img);
		filelumen.close();
	}
	img.release();
}

//--------------------------------------------------C O R E - F U N C T I O N S--------------------------------------------------------------


ultrasound::ResultOfProcess ultrasound::sortUsingPolarCoordinates(vector<Point2f>* p, int iter, Point2f* center, Mat image, int skinDistance) {

	transform(p->begin(), p->end(), p->begin(), std::bind2nd(std::minus<Point2f>(), *center));

	Mat xpts(p->size(), 1, CV_32F, &p->at(0).x, 2 * sizeof(float));
	Mat ypts(p->size(), 1, CV_32F, &p->at(0).y, 2 * sizeof(float));

	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);

	vector<Point2f> polarXY;
	for (int i = 0; i < p->size(); i++) {
		polarXY.push_back(Point2f(magnitude.at<Float32>(i), angle.at<Float32>(i)));
	}

	sort(polarXY.begin(), polarXY.end(), sortBasedYCoordinate);

	/*for (int i = 0; i < polarXY.size(); i++) {
		cout << "Angle =  " << polarXY[i].y * 180 / CV_PI << " R = " << polarXY[i].x << endl;
	}*/


	Mat mag(polarXY.size(), 1, CV_32F, &polarXY[0].x, 2 * sizeof(float));
	Mat ang(polarXY.size(), 1, CV_32F, &polarXY[0].y, 2 * sizeof(float));

	Mat xnew, ynew;
	polarToCart(mag, ang, xnew, ynew);

	vector<Point2f> pp;

	for (int i = 0; i < p->size(); i++) {
		pp.push_back(Point2f(xnew.at<Float32>(i), ynew.at<Float32>(i)));
	}

	vector<Point2f> smoothed = smoothContour(sortBasedEuclideanDistance(pp), 0); //num_spline = 50;

	if (!smoothed.empty()) {
		transform(smoothed.begin(), smoothed.end(), smoothed.begin(), std::bind2nd(std::plus<Point2f>(), *center));
		this->lumenPoints.push_back(smoothed);
		this->levelsetSize = round(smoothed.size() / 3);//-----> re-estimate levelsetSize
	}
	else if (iter > this->initialFrame) {
		this->lumenPoints.push_back(this->lumenPoints.back());
	}
	else if (iter == this->initialFrame) {
		return ResultOfProcess::FAILURE;
	}

	LoggerMessage("Bladder borders were detected successfully!");

	xpts.release();
	ypts.release();
	mag.release();
	ang.release();
	xnew.release();
	ynew.release();

	vector<Point2f>().swap(polarXY);
	vector<Point2f>().swap(pp);

	return ResultOfProcess::SUCCESS;
}


ultrasound::ResultOfProcess ultrasound::sortClockwise(vector<Point2f>* p, Point2f* center, int iter, bool recalculate) {

	transform(p->begin(), p->end(), p->begin(), std::bind2nd(std::minus<Point2f>(), *center));

	Mat xpts(p->size(), 1, CV_32F, &p->at(0).x, 2 * sizeof(float));
	Mat ypts(p->size(), 1, CV_32F, &p->at(0).y, 2 * sizeof(float));

	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);

	vector<Point2f> polarXY, sortedPolarXY;
	for (int i = 0; i < p->size(); i++) {
		polarXY.push_back(Point2f(magnitude.at<Float32>(i), angle.at<Float32>(i)));
	}
	sort(polarXY.begin(), polarXY.end(), sortBasedYCoordinate);

	////--------------store the first point--------------
	//investigate if there are points with greater radius than polarXY[0] close to the polarXY[0]
	int count = 0;
	for (int i = 1; i < 5; i++) {
		if (polarXY[i].x - polarXY[count].x > 4) {
			count = i;
		}
	}

	sortedPolarXY.push_back(polarXY[count]);
	if (count == 0) {
		
		polarXY.erase(polarXY.begin());
	}
	else {
		rotate(polarXY.begin(), polarXY.begin()+count-1, polarXY.end());
	}

	////--------------store the next 3 points--------------
	vector<double> nextSetPoints;
	int minIndex;

	//--------------store the rest points--------------
	while (!polarXY.empty()) {
		for (int i = 0; i < polarXY.size(); i++) {
			//calculate euclidean distance in polar ccordinates
			nextSetPoints.push_back(sqrt(pow(polarXY[i].x, 2) + pow(sortedPolarXY.back().x, 2)  - 2*polarXY[i].x*sortedPolarXY.back().x * cos(polarXY[i].y - sortedPolarXY.back().y)));
		}
		auto min = std::min_element(nextSetPoints.begin(), nextSetPoints.end());
		minIndex = distance(nextSetPoints.begin(), min);
		sortedPolarXY.push_back(polarXY[minIndex]);
		polarXY.erase(polarXY.begin() + minIndex);
		vector<double>().swap(nextSetPoints);
	}

	Mat mag(sortedPolarXY.size(), 1, CV_32F, &sortedPolarXY[0].x, 2 * sizeof(float));
	Mat ang(sortedPolarXY.size(), 1, CV_32F, &sortedPolarXY[0].y, 2 * sizeof(float));

	Mat xnew, ynew;
	polarToCart(mag, ang, xnew, ynew);

	vector<Point2f> pp;

	for (int i = 0; i < p->size(); i++) {
		pp.push_back(Point2f(xnew.at<Float32>(i), ynew.at<Float32>(i)));
	}

	vector<Point2f> smoothed = smoothContour(pp, 0, true); //num_spline = 50;
	transform(smoothed.begin(), smoothed.end(), smoothed.begin(), std::bind2nd(std::plus<Point2f>(), *center));

	if (!recalculate) {
		if (!smoothed.empty()) {
			this->lumenPoints.push_back(smoothed);
			this->levelsetSize = round(smoothed.size() / 3);//-----> re-estimate levelsetSize
		}
		else if (iter > this->initialFrame) {
			this->lumenPoints.push_back(this->lumenPoints.back());
		}
		else if (iter == this->initialFrame) {
			return ResultOfProcess::FAILURE;
		}
	}
	else {
		vector<Point2f>().swap(this->contourForFix);
		smoothed.swap(contourForFix);
	}
	

	LoggerMessage("Bladder borders were detected successfully!");

	xpts.release();
	ypts.release();
	mag.release();
	ang.release();
	xnew.release();
	ynew.release();
	vector<Point2f>().swap(polarXY);
	vector<Point2f>().swap(sortedPolarXY);
	vector<Point2f>().swap(pp);
	vector<Point2f>().swap(smoothed);

	return ResultOfProcess::SUCCESS;
}


Point2f ultrasound::getCenterOfGravity(vector<Point2f> points) {

	Mat matContour = Mat(points);
	vector<Point2f> pp;
	convexHull(points, pp, false);

	cv::Point2f Coord;
	cv::Moments mm = cv::moments(pp, false);
	double moment10 = mm.m10;
	double moment01 = mm.m01;
	double moment00 = mm.m00;
	Coord.x = moment10 / moment00;
	Coord.y = moment01 / moment00;

	vector<Point2f>().swap(pp);
	vector<Point2f>().swap(points);
	matContour.release();
	return Coord;
}


Point2f ultrasound::getCenterOfGravity(vector<Point> points) {

	Mat matContour = Mat(points);
	vector<Point> pp;
	convexHull(points, pp, false);
	cv::Moments mm = cv::moments(pp, false);

	vector<Point>().swap(pp);
	vector<Point>().swap(points);
	matContour.release();
	return Point2f(mm.m10 / mm.m00, mm.m01 / mm.m00);
}


Point2f ultrasound::getCenterOfMass(vector<Point2f> points) {

	Point2f center = accumulate(points.begin(), points.end(), Point2f(0.0, 0.0));
	center.x /= points.size();
	center.y /= points.size();

	vector<Point2f>().swap(points);
	return center;
}


ultrasound::ResultOfProcess ultrasound::centerAndPointsOfContour(Mat processed, vector<Point2f> *points, Point2f *center, Point *highestWhitePixel)
{
	Mat whitePixels;   // output, locations of non-zero pixels
	
	findNonZero(processed, whitePixels);

	if (whitePixels.empty()) {
		LoggerMessage("No Bladder points were detected");
		return ResultOfProcess::FAILURE;
	} 

	highestWhitePixel->x = whitePixels.at<Point>(0).x;
	highestWhitePixel->y = whitePixels.at<Point>(0).y;

	center->x = 0.0;
	center->y = 0.0;

	for (int i = 0; i < whitePixels.total(); i++) {
		Point lumen_p = Point(whitePixels.at<Point>(i).x, whitePixels.at<Point>(i).y);
			points->push_back(Point2f(whitePixels.at<Point>(i).x, whitePixels.at<Point>(i).y));
			center->x += whitePixels.at<Point>(i).x;
			center->y += whitePixels.at<Point>(i).y;
	}
	if (whitePixels.total() != 0) {
		*center = getCenterOfMass(*points);
	}
	else {
		return ResultOfProcess::FAILURE;
	}


	whitePixels.release();

	return ResultOfProcess::SUCCESS;
}



CImg<unsigned char>* ultrasound::cvImgToCImg(Mat &cvImg)
{
	CImg<unsigned char> * result = new CImg<unsigned char>(cvImg.cols, cvImg.rows);

	for (int x = 0; x < cvImg.cols; ++x)
		for (int y = 0; y < cvImg.rows; ++y)
			(*result)(x, y) = cvImg.at<uchar>(y, x);

	return result;
}

Mat ultrasound::CImgtoCvImg(CImg<unsigned char> img)
{
	Mat result = Mat(img.height(), img.width(), CV_8U);

	for (int x = 0; x < result.cols; ++x)
		for (int y = 0; y < result.rows; ++y)
			result.at<uchar>(y, x) = (img)(x, y);

	return result;
}


CImg<unsigned char> ultrasound::circle_levelset(int height, int width, const array<int, 2>& center, double radius, double scalerow)
{
	CImg<unsigned char> res(width, height);
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int diffy = (i - center[0]);
			int diffx = (j - center[1]);
			res(j, i) = (radius*radius - (diffx*diffx + diffy * diffy)) > 0;
		}
	}
	return res;
}

//sorted.push_back(sorted[0]);

vector<Point2f> ultrasound::sortBasedEuclideanDistance(vector<Point2f> points) {
	vector<Point2f> sorted;

	//keep the first point
	sorted.push_back(points[0]);
	points.erase(points.begin());
	
	Point2f po;
	while (!points.empty()) {
		po = sorted.back();
		sort(points.begin(), points.end(), EuclideanDistance(po));
		if (sqrt(pow(points[0].x - sorted.back().x, 2) + pow(points[0].y - sorted.back().y, 2)) < 20) {
			sorted.push_back(points[0]);
		}
		points.erase(points.begin());
	}
	
	return sorted;
}

int ultrasound::findLongestVector(vector<vector<Point>> vec) {
	int max = 0;
	for (int i = 1; i < vec.size(); i++) {
		if (vec[i].size() > vec[max].size()) {
			max = i;
		}
	}
	return max;
}

bool ultrasound::IsClockwise(vector<Point2f> points)
{/*
	double sum = 0.0;
	for (int i = 0; i < points.size()-1; i++) {
		sum += (points[i+1].x - points[i].x) * (points[i + 1].y + points[i].y);
	}
	return sum > 0.0;*/


	double area = 0;
	int j;
	for (int i = 0; i < points.size(); i++) {
		j = (i + 1) % points.size();
		area += points[i].x * points[j].y;
		area -= points[j].x * points[i].y;
	}
	vector<Point2f>().swap(points);
	return (area/2)>0;
}


vector<Point2f> ultrasound::smoothContour(vector<Point2f> contour, int num_spline, bool closedContour) {

	if (closedContour) contour.push_back(contour[0]);

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

	vector<Point2f>().swap(contour);
	contour = smoothCurve(cl_3D, num_spline);

	//free some memory
	vector<vector<vector<double>>>().swap(cl_3D);
	vector<vector<double>>().swap(temp_cl_3D);
	vector<double>().swap(xx);
	vector<double>().swap(yy);
	vector<double>().swap(zz);

	return contour;

}


vector<Point2f> ultrasound::interpolateConvexPoints(vector<Point2f> p, interpolationMethod method) {

	vector<double> diff;
	for (int i = 0; i < p.size() - 2; i++) {
		diff.push_back(abs((p[i].y - p[i + 1].y) / (p[i].x - p[i + 1].x) - (p[i + 1].y - p[i + 2].y) / (p[i + 1].x - p[i + 2].x)));
	}
	auto min = min_element(diff.begin(), diff.end());
	int minIndex = distance(diff.begin(), min);
	if (minIndex > 0) {
		rotate(p.begin(), p.begin() + minIndex + 1, p.end());
	}

	vector<double>().swap(diff);


	Point2f center = getCenterOfGravity(p);

	transform(p.begin(), p.end(), p.begin(), std::bind2nd(std::minus<Point2f>(), center));

	Mat xpts(p.size(), 1, CV_32F, &p.at(0).x, 2 * sizeof(float));
	Mat ypts(p.size(), 1, CV_32F, &p.at(0).y, 2 * sizeof(float));
	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);
	vector<Point2f> polarXY, interpolated;
	for (int i = 0; i < p.size(); i++) {
		polarXY.push_back(Point2f(magnitude.at<float>(i), angle.at<float>(i)));
	}

	sort(polarXY.begin(), polarXY.end(), sortBasedYCoordinate);

	polarXY.push_back(Point2f(polarXY[0].x, polarXY[0].y + 2 * CV_PI));
	double minAngle = polarXY[0].y;
	double maxAngle = polarXY.back().y;


	//------------------------LIBRARY 1--------------------
	//Curve* curve = new BSpline();
	////Curve* curve = new Bezier();
	////Curve* curve = new CatmullRom();
	//curve->set_steps(100); // generate 100 interpolate points between the last 4 way points

	//for (int i = 0; i < polarXY.size(); i++) {
	//	curve->add_way_point(Vector(polarXY[i].x, polarXY[i].y, 0));
	//}
	////curve->add_way_point(Vector(polarXY[0].x, polarXY[0].y, 0));

	//for (int i = 0; i < curve->node_count(); ++i) {
	//	interpolated.push_back(Point2f(curve->node(i).x, curve->node(i).y));
	//}
	//delete curve;


	//------------------------LIBRARY 2--------------------
	/*vector<double> X, Y;
	for (int i = 0; i < polarXY.size(); i++) {
		X.push_back(polarXY[i].y);
		Y.push_back(polarXY[i].x);
	}

	tk::spline s;
	s.set_points(X, Y);

	for (double i = 0; i <= 2 * CV_PI; i = i + 0.02) {
		interpolated.push_back(Point2f(s(i), i));
	}*/

	//------------------------LIBRARY 3--------------------
	string xx = "[" + to_string(polarXY[0].x);
	string yy = "[" + to_string(polarXY[0].y);

	for (int i = 1; i < polarXY.size(); i++) {
		xx += "," + to_string(polarXY[i].x);
		yy += "," + to_string(polarXY[i].y);
	}
	xx += "]";
	yy += "]";

	alglib::real_1d_array X = xx.c_str();
	alglib::real_1d_array Y = yy.c_str();
	alglib::ae_int_t n = 100;
	alglib::spline1dinterpolant s;

	switch (method) {
	case interpolationMethod::LINEAR:
		alglib::spline1dbuildlinear(Y, X, s);
		break;
	case interpolationMethod::AKIMA:
		alglib::spline1dbuildakima(Y, X, s);
		break;
	case interpolationMethod::CATMULLROM:
		alglib::spline1dbuildcatmullrom(Y, X, s);
		break;
	case interpolationMethod::CUBIC:
		alglib::spline1dbuildcubic(Y, X, s);
		break;
	case interpolationMethod::MONOTONE:
		alglib::spline1dbuildmonotone(Y, X, s);
		break;
	}

	for (double i = minAngle; i <= maxAngle; i = i + 0.08) {
		interpolated.push_back(Point2f(spline1dcalc(s, i), i));
	}
	//-----------------------------------------------------


	Mat mag(interpolated.size(), 1, CV_32F, &interpolated[0].x, 2 * sizeof(float));
	Mat ang(interpolated.size(), 1, CV_32F, &interpolated[0].y, 2 * sizeof(float));

	Mat xnew, ynew;
	polarToCart(mag, ang, xnew, ynew);

	vector<Point2f> pp;

	for (int i = 0; i < interpolated.size(); i++) {
		pp.push_back(Point2f(xnew.at<float>(i), ynew.at<float>(i)));
	}

	transform(pp.begin(), pp.end(), pp.begin(), std::bind2nd(std::plus<Point2f>(), center));


	//vector<double>().swap(X); //if use libray 2
	//vector<double>().swap(Y);
	vector<Point2f>().swap(polarXY);
	vector<Point2f>().swap(interpolated);
	xpts.release();
	ypts.release();
	magnitude.release();
	angle.release();
	mag.release();
	ang.release();
	return pp;
}


// Cubic spline interpolation to smooth centerline
vector<Point2f> ultrasound::smoothCurve(vector<vector<vector<double>>> centerline, int num_spline) { //vector<vector<double>>
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
			planeDistance = 10; // sampling distance 54
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

	//t.push_back(Point(centerline[0][0][0], centerline[0][0][1])); //push the first point into the vector

	for (int i = 0; i < smoothCtr[0].size(); i++) {
		t.push_back(Point2f(smoothCtr[0][i][0], smoothCtr[0][i][1]));
	}

	//vector<Point2f>().swap(t);

	//free memory
	vector<vector<vector<double>>>().swap(smoothCtr);

	return t;// smoothCtr[0]; //return only the first centerline 
}