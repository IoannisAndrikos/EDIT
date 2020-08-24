//// openCVInstallation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
//
//#include <stdlib.h> 
//#include <opencv2/core/core.hpp>
//#include "opencv2/imgproc/imgproc.hpp"
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/opencv.hpp>
//#include <iostream>
//#include <string>
//#include <opencv2/opencv.hpp>
//#include <iostream>
//#include <numeric>
//#include <fstream>
//#include <iostream>
//#include <math.h>
//#include <experimental/filesystem>
//
//
//// ----------------------------VTK Triangulation----------------------
//#include "vtk-9.0/vtkSmartPointer.h"
//#include "vtk-9.0/vtkCardinalSpline.h"
//#include "vtk-9.0/vtkPoints.h"
//#include "vtk-9.0/vtkPolyData.h"
//#include "vtk-9.0/vtkDoubleArray.h"
////---------------------------------------------------------------------
//
//
//#include "morphsnakes.h"
//#include "CImg.h"
//#include "DicomReader.h"
//#include "spline.h"
//
//
//
//using namespace std;
//using namespace cv;
//namespace ms = morphsnakes;
//using namespace cimg_library;
//using namespace tk;
//
//
//Point userPoint = { 0,0 };
//
//
////---------------------------------------------------------
//// Cubic spline interpolation to smooth centerline
//vector<Point> smoothCenterline(vector<vector<vector<double>>> centerline) { //vector<vector<double>>
//	vector<vector<vector<double>>> smoothCtr;
//	for (int ii = 0; ii < centerline.size(); ++ii) {
//		// Read centerline into vtk points
//		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//		for (int i = 0; i < centerline[ii].size(); ++i) {
//			points->InsertNextPoint(centerline[ii][i][0], centerline[ii][i][1], centerline[ii][i][1]); // ATTENTION TO X-Y ORDER
//		}
//		// Create polydata using the set of path points
//		vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
//		data->SetPoints(points);
//		double resamplingParameter = 1; // 10 Sampling resolution
//		double distance = 0; // distance from start to end of centerline
//							 // Calculate distance
//		for (int i = 0; i < data->GetNumberOfPoints() - 1; ++i) {
//			double p1[3], p2[3];
//			data->GetPoint(i, p1);
//			data->GetPoint(i + 1, p2);
//			double diff_x = pow(p2[0] - p1[0], 2.0);
//			double diff_y = pow(p2[1] - p1[1], 2.0);
//			double diff_z = pow(p2[2] - p1[2], 2.0);
//			distance += sqrt(diff_x + diff_y + diff_z);
//		}
//		// Create arrays to store samples for each axis
//		vtkSmartPointer<vtkDoubleArray> x, y, z;
//		x = vtkSmartPointer<vtkDoubleArray>::New();
//		x->SetNumberOfComponents(1);
//		y = vtkSmartPointer<vtkDoubleArray>::New();
//		y->SetNumberOfComponents(1);
//		z = vtkSmartPointer<vtkDoubleArray>::New();
//		z->SetNumberOfComponents(1);
//		double p0[3];
//		// Sample centerline and insert points to respectable arrays
//		for (int i = 0; i < data->GetNumberOfPoints(); i = i + resamplingParameter) {
//			data->GetPoint(i, p0);
//			x->InsertNextValue(p0[0]);
//			y->InsertNextValue(p0[1]);
//			z->InsertNextValue(p0[2]);
//		}
//		// Create a spline object for each axis
//		vtkSmartPointer<vtkCardinalSpline> xspline = vtkSmartPointer<vtkCardinalSpline>::New();
//		xspline->ClosedOff(); // Curve is not closed
//		vtkSmartPointer<vtkCardinalSpline> yspline = vtkSmartPointer<vtkCardinalSpline>::New();
//		yspline->ClosedOff();
//		vtkSmartPointer<vtkCardinalSpline> zspline = vtkSmartPointer<vtkCardinalSpline>::New();
//		zspline->ClosedOff();
//		double sampledPoints = (double)x->GetNumberOfTuples();
//		// Spline interpolation
//		// t represents the position of each point on the curve
//		// e.g.: 0: start point, 1: end point, 0.5: mid point
//		for (int i = 0; i < sampledPoints - 1; ++i) {
//			double t = (double)i / sampledPoints;
//			double xpoint = x->GetTuple1(i);
//			xspline->AddPoint(t, xpoint); // x axis
//			double ypoint = y->GetTuple1(i);
//			yspline->AddPoint(t, ypoint); // y axis
//			double zpoint = z->GetTuple1(i);
//			zspline->AddPoint(t, zpoint); // z axis
//		}
//		// Set final point for interpolation
//		double xpoint = x->GetTuple1(sampledPoints - 1);
//		xspline->AddPoint(1, xpoint);
//		double ypoint = y->GetTuple1(sampledPoints - 1);
//		yspline->AddPoint(1, ypoint);
//		double zpoint = z->GetTuple1(sampledPoints - 1);
//		zspline->AddPoint(1, zpoint);
//		double evaluationPoints = 10000, dt = 1 / evaluationPoints; // 10000
//		// Create an array to store interpolated points
//		double **p = new double *[evaluationPoints];
//		for (int i = 0; i < evaluationPoints; ++i) {
//			p[i] = new double[3];
//		}
//		int i = 0;
//		double t = 0;
//		// Evaluate splines at small intervals and store points
//		while (i < evaluationPoints) {
//			p[i][0] = xspline->Evaluate(t);
//			p[i][1] = yspline->Evaluate(t);
//			p[i][2] = zspline->Evaluate(t);
//			i++; t = t + dt;
//		}
//		// Calculate distance between the spline points
//		double smoothed_distance = 0;
//		double * newdst = new double[evaluationPoints - 1];
//		for (int i = 0; i < evaluationPoints - 1; ++i) {
//			double diff_x = pow(p[i + 1][0] - p[i][0], 2.0);
//			double diff_y = pow(p[i + 1][1] - p[i][1], 2.0);
//			double diff_z = pow(p[i + 1][2] - p[i][2], 2.0);
//			smoothed_distance += sqrt(diff_x + diff_y + diff_z);
//			newdst[i] = sqrt(diff_x + diff_y + diff_z);
//		}
//		//cout << smoothed_distance << " ";
//
//		double planeDistance = ceil(smoothed_distance) / 100;//0.5; // sampling distance 54
//		//double planeDistance = ceil(dims[0] * 15 / 1024);
//
//		double interval = 0;
//		vector<int> ind; // vector to store indexes of equidistant points
//		for (int i = 0; i < evaluationPoints - 1; ++i) {
//			interval += newdst[i];
//			if (interval >= planeDistance - 0.01) { // select points with equal distance on curve
//				interval = 0;
//				ind.push_back(i);
//			}
//		}
//		// Store equidistant points
//		vector<vector<double>> equidistantPoints;
//		for (int i = 0; i < ind.size(); ++i) {
//			vector<double> point;
//			point.push_back(p[ind[i]][0]); point.push_back(p[ind[i]][1]); point.push_back(p[ind[i]][2]);
//			equidistantPoints.push_back(point);
//		}
//		smoothCtr.push_back(equidistantPoints); // store interpolated centerline
//	}
//
//
//	//convert point from double to Point (opencv)
//	vector<Point> t;
//
//	t.push_back(Point(centerline[0][0][0], centerline[0][0][1])); //push the first point into the vector
//
//	for (int i = 0; i < smoothCtr[0].size(); i++) {
//		t.push_back(Point(smoothCtr[0][i][0], smoothCtr[0][i][1]));
//	}
//
//	return t;// smoothCtr[0]; //return only the first centerline 
//}
//
//
//
//struct sortclass {
//	bool operator() (cv::Point2f pt1, cv::Point2f pt2) { return (pt1.y < pt2.y); }
//} sortObject;
//
//
//
//cimg_library::CImg<unsigned char> *cvImgToCImg(cv::Mat &cvImg)
//{
//	cimg_library::CImg<unsigned char> * result = new cimg_library::CImg<unsigned char>(cvImg.cols, cvImg.rows);
//
//	for (int x = 0; x < cvImg.cols; ++x)
//		for (int y = 0; y < cvImg.rows; ++y)
//			(*result)(x, y) = cvImg.at<uchar>(y, x);
//
//	return result;
//}
//
//cv::Mat CImgtoCvImg(cimg_library::CImg<unsigned char> img)
//{
//	Mat result = cv::Mat(img.height(), img.width(), CV_8U);
//
//	for (int x = 0; x < result.cols; ++x)
//		for (int y = 0; y < result.rows; ++y)
//			result.at<uchar>(y, x) = (img)(x, y);
//
//	return result;
//}
//
//
//
//
//CImg<unsigned char> circle_levelset(int height, int width,
//	const std::array<int, 2>& center,
//	double radius,
//	double scalerow = 1.0)
//{
//	CImg<unsigned char> res(width, height);
//	for (int i = 0; i < height; ++i)
//	{
//		for (int j = 0; j < width; ++j)
//		{
//			int diffy = (i - center[0]);
//			int diffx = (j - center[1]);
//			res(j, i) = (radius*radius - (diffx*diffx + diffy * diffy)) > 0;
//		}
//	}
//
//	return res;
//}
//
//vector<vector<vector<double>>> sortBasedEuclideanDistance(vector<Point2f> points) {
//
//	vector<Point2f> sorted;
//	for (int i = 0; i < 3; i++) {
//		sorted.push_back(points[i]);
//		points[i] = { 256, 256 };
//
//	}
//
//
//	while (sorted.size() < points.size()) {
//		int index;
//		vector<double> d;
//		vector<double>::iterator iter;
//		Point2f po = sorted.back();
//		for (int i = 0; i < points.size(); i++) {
//			 double dist =  sqrt(pow(po.x - points[i].x, 2) + pow(po.y - points[i].y, 2));
//			 d.push_back(dist);
//		}
//		double min_el = *min_element(d.begin(), d.end());
//		iter = find(d.begin(), d.end(), min_el);
//		index = std::distance(d.begin(), iter);
//		sorted.push_back(points[index]);
//		points[index] = { 256, 256 };
//	}
//
//	sorted.push_back(sorted[0]);
//
//
//	vector<vector<vector<double>>> cl_3D;
//	vector<vector<double>> temp_cl_3D;
//	vector<double> xx;
//	vector<double> yy;
//	vector<double> zz;
//
//	for (int i = 0; i < sorted.size(); i++) {
//		vector<double> temp_temp_cl_3D;
//		temp_temp_cl_3D.push_back(sorted[i].x);
//		temp_temp_cl_3D.push_back(sorted[i].y);
//		temp_cl_3D.push_back(temp_temp_cl_3D);
//	}
//
//	cl_3D.push_back(temp_cl_3D);
//
//
//	return cl_3D;
//}
//
//void sortUsingPolarCoordinates(vector<Point> p1, int iter, Point center) {
//
//	vector<Point2f> p;
//	for (int i = 0; i < p1.size(); i++) {
//		//move points to (0,0) -> click point on image
//		p.push_back(Point2f(p1[i].x-center.x, p1[i].y - center.y));
//	}
//
//	Mat xpts(p.size(), 1, CV_32F, &p[0].x, 2 * sizeof(float));
//	Mat ypts(p.size(), 1, CV_32F, &p[0].y, 2 * sizeof(float));
//
//	Mat magnitude, angle;
//	cartToPolar(xpts, ypts, magnitude, angle);
//
//	vector<Point2f> polarXY;
//	for (int i = 0; i < p.size(); i++) {
//		polarXY.push_back(Point2f(magnitude.at<Float32>(i),angle.at<Float32>(i)));
//	}
//
//	sort(polarXY.begin(), polarXY.end(), sortObject);
//
//	Mat mag(polarXY.size(), 1, CV_32F, &polarXY[0].x, 2 * sizeof(float));
//	Mat ang(polarXY.size(), 1, CV_32F, &polarXY[0].y, 2 * sizeof(float));
//
//	Mat xnew, ynew;
//	polarToCart(mag, ang, xnew, ynew);
//
//	ofstream myfile2;
//	myfile2.open(String("C:/Users/Legion Y540/Desktop/j2/" + to_string(iter)) + ".txt");
//
//	vector<Point2f> pp;
//
//	for (int i = 0; i < p.size(); i++) {
//		pp.push_back(Point2f(xnew.at<Float32>(i), ynew.at<Float32>(i)));
//	}
//
//
//	/*vector<vector<vector<double>>> cl_3D;
//	vector<vector<double>> temp_cl_3D;
//	vector<double> xx;
//	vector<double> yy;
//	vector<double> zz;
//
//	for (int i = 0; i < p.size(); i++) {
//		vector<double> temp_temp_cl_3D;
//		temp_temp_cl_3D.push_back(xnew.at<Float32>(i));
//		temp_temp_cl_3D.push_back(ynew.at<Float32>(i));
//		temp_cl_3D.push_back(temp_temp_cl_3D);
//	}
//
//	cl_3D.push_back(temp_cl_3D);*/
//
//	vector<Point> smothed = smoothCenterline(sortBasedEuclideanDistance(pp));
//	//vector<Point> smothed = smoothCenterline(cl_3D);
//	
//	for (int i = 0; i < smothed.size(); i++) { //p.size()
//		//myfile2 << round(xnew.at<Float32>(i)) + center.x << " " << round(ynew.at<Float32>(i)) + center.y << endl;
//		myfile2 << smothed[i].x + center.x << " " << smothed[i].y + center.y << endl;
//	}
//
//	myfile2.close();
//}
//
//
//template<class T> ms::NDImage<T, 2> cimg2ndimage(CImg<T>& img)
//{
//	ms::Shape<2> shape = { img.height(), img.width() };
//	ms::Stride<2> stride = { img.width() * sizeof(T), sizeof(T) };
//
//	return ms::NDImage<T, 2>(img.data(), shape, stride);
//}
//
//
//
//vector<Point> lakes(Mat img, Point point, int iter)
//{
//	equalizeHist(img, img);
//	CImg<double>  imgC = *cvImgToCImg(img);
//
//
//	//Initialize embedding function
//	auto embedding = circle_levelset(imgC.height(), imgC.width(), { userPoint.y , userPoint.x }, 40); //{ 256 , 256 }
//	String ss = String("C:/Users/Legion Y540/Desktop/click_points/" + to_string(iter)) + ".png";
//	(embedding * 255).save_png("click_points");
//
//	Mat click_point = CImgtoCvImg((embedding * 255));
//	imwrite(ss, click_point);
//
//	//Morphological ACWE
//	ms::MorphACWE<double, 2> macwe(cimg2ndimage(embedding), cimg2ndimage(imgC), 3);
//	for (int i = 0; i < 200; ++i)
//		macwe.step();
//
//	// Save results
//	//(embedding * 255).save_png("filtered_image.png");
//
//	Mat mlk = CImgtoCvImg((embedding * 255));
//
//	Canny(mlk, mlk, 0, 0, 3);
//	threshold(mlk, mlk, 100, 255, THRESH_BINARY);
//
//	/*Mat a = toPolar(mlk, point);
//	String contourspath = String("C:/Users/Legion Y540/Desktop/contours/" + to_string(iter)) + ".png";
//	imwrite(contourspath, a);*/
//
//	//--------------------------------------------POLAR for skin---------------------------------------------------------
//	Mat contour;
//	Point2f center(point.x, point.y);
//	int flags = INTER_LINEAR + WARP_FILL_OUTLIERS;
//	linearPolar(mlk, contour, center, 256, flags );
//	int pix = 180;
//	Mat out = Mat::zeros(contour.size(), contour.type());
//	contour(Rect(0, 0, contour.cols -pix, contour.rows)).copyTo(out(Rect(pix, 0, contour.cols -pix, contour.rows)));
//	Mat outCartesian;
//	linearPolar(out, outCartesian, center, 256, flags + WARP_INVERSE_MAP);
//	threshold(out, out, 100, 255, THRESH_BINARY);
//
//	String contourspath = String("C:/Users/Legion Y540/Desktop/contours/" + to_string(iter)) + ".png";
//	imwrite(contourspath, outCartesian); //outCartesian
//	//-----------------------------------------------------------------------------------------------------------
//
//	Mat locations;   // output, locations of non-zero pixels
//	findNonZero(mlk, locations);
//
//
//	vector<Point> lumen_xy;
//
//	for (int i = 0; i < locations.total(); i++) {
//		Point lumen_p = Point(locations.at<Point>(i).x, locations.at<Point>(i).y);
//		lumen_xy.push_back(lumen_p);
//	}
//
//	sortUsingPolarCoordinates(lumen_xy, iter, point);
//
//	return lumen_xy;
//}
//
//void initiatePoint(int event, int x, int y, int, void* imgptr) {
//
//	Mat & img = (*(Mat*)imgptr); // first cast, then deref
//	Point pt1 = Point(x, y);
//	if (event == EVENT_LBUTTONDOWN) {
//		userPoint = pt1;
//		circle(img, pt1, 2, Scalar(0, 0, 200), 3, 8, 0);
//	}
//	imshow("window", img);
//	waitKey(1);
//}
//
//
//Point getMeanPoint(vector<Point> points, Mat img, int iter, Point userPoint) {
//
//	int imCenter = 0; // floor(img.size().height / 2); // very important --> shift center of image to (0,0)
//
//	//sortPoints(points, iter);
//
//	ofstream myfile;
//	myfile.open(String("C:/Users/Legion Y540/Desktop/output_points/" + to_string(iter)) + ".txt");
//	
//
//	int sum_x =0, sum_y=0;
//	int s = points.size();
//
//	Mat black = Mat::zeros(img.size(), img.type());
//
//	for (int i = 0; i < s; i++) {
//		int x = points[i].x - imCenter;
//		int y = points[i].y - imCenter;
//
//		myfile << x  << " " <<y  << "\n"; // write points to .txt
//
//		//------------------------only for view-----------------------
//		sum_x += points[i].x; //do not remove imCenter
//		sum_y += points[i].y; //do not remove imCenter
//		img.at<uchar>(points[i].y, points[i].x) = 255;
//		//------------------------------------------------------------
//	}
//	myfile.close();
//
//	imwrite(String("C:/Users/Legion Y540/Desktop/segmented_images/" + to_string(iter)) + ".bmp", img) ;
//
//
//	if (s != 0) {
//		Point p = { sum_x / s, sum_y / s };
//
//		return p;
//	}
//	else {
//		return userPoint;
//	}
//	
//}
//
//
//
//
//
//int main(int argv, char* argc)
//{
//
//
//	//---------------------------------------------READ DICOM-----------------------------------
//
//	string name = experimental::filesystem::path("C:/Users/Legion Y540/Desktop/EDIT Data/US_R73_2020-06-23-17-07-10.dcm").filename().generic_string();
//	size_t lastindex = name.find_last_of(".");
//	string rawname = name.substr(0, lastindex);
//	cout << rawname << endl;
//
//
//	//C:/Users/Legion Y540/Desktop/EDIT Data/US_R74_2020-06-23-13-54-24.dcm
//	wstring dcmpath = (L"C:/Users/Legion Y540/Desktop/EDIT Data/US_R73_2020-06-23-17-07-10.dcm"); //insert DICOM file
//	wstring bmppath = (L"C:/Users/Legion Y540/Desktop/output_images/");
//
//	CDicomReader *reader = new CDicomReader();
//	reader->dcmimage_bmp(dcmpath.c_str(), bmppath.c_str()); // read the given dicom and save images into the given folder
//
//	vector<double> Tags;
//	Tags = reader->GetDicomInfo(dcmpath.c_str()); //save dicom tags
//
//	vector<Mat> dcmImages = reader->dcmimage_Mat(dcmpath.c_str(), Tags[2], Tags[3], Tags[4], Tags[5]);
//
//
//
//	reader->~CDicomReader();
//	//-----------------------------------------------------------------------------------------
//
//	cout << "PhysicalDetaX = " << to_string(Tags[0]) << endl;
//	cout << "PhysicalDetaY = " << to_string(Tags[1]) << endl;
//	cout << "PhysicalDetaY = " << to_string(Tags[2]) << endl;
//	cout << "PhysicalDetaY = " << to_string(Tags[3]) << endl;
//	cout << "PhysicalDetaY = " << to_string(Tags[4]) << endl;
//	cout << "PhysicalDetaY = " << to_string(Tags[5]) << endl;
//
//
//	int firstImage = 6;
//	int lastImage = 76;
//
//	vector<Point> segm;
//
//	for (int i = firstImage; i <= lastImage; i++) {
//
//		cout << "I am starting " << i << endl;
//		String path = "C:/Users/Legion Y540/Desktop/output_images/" + to_string(i);
//		path += ".bmp";
//		cout << path << endl;
//		Mat img = imread(path, IMREAD_GRAYSCALE); // IMREAD_GRAYSCALE); //
//
//		//Mat img = dcmImages[i]; // IMREAD_GRAYSCALE); //
//
//		//----------------------------Crop--------------------
//		int x = round(Tags[3]) - round(Tags[2]);
//		int y = round(Tags[5]) - round(Tags[4]);
//		img = img(Rect(round(Tags[2]), round(Tags[4]), x, y));
//		imwrite(path, img);
//		//----------------------------------------------------
//
//		Mat inputforPoint = img;
//
//		if (i == firstImage) {
//			namedWindow("window", 1);
//			setMouseCallback("window", initiatePoint, &inputforPoint);
//			imshow("window", inputforPoint);
//			waitKey(0);
//			//swap(userPoint.x, userPoint.y);
//		}
//		segm = lakes(img, userPoint, i);
//		userPoint = getMeanPoint(segm, img, i, userPoint);
//
//		cout << "I found borders of " << i << endl;
//	}
//}