#pragma once

#include "tumor.h"
#include "messages.h"

tumor::tumor() {}
tumor::~tumor() {}

void tumor::sortClockwise(vector<vector<vector<Point2f>>> points) {
	vector<vector<Point2f>>().swap(tumorPoints);

	points[0][0] = smoothContour(points[0][0], 100, true);
	Point2f center = getCenterOfMass(points[0][0]);

	transform(points[0][0].begin(), points[0][0].end(), points[0][0].begin(), std::bind2nd(std::minus<Point2f>(), center));

	Mat xpts(points[0][0].size(), 1, CV_32F, &points[0][0][0].x, 2 * sizeof(float));
	Mat ypts(points[0][0].size(), 1, CV_32F, &points[0][0][0].y, 2 * sizeof(float));
	Mat magnitude, angle;
	cartToPolar(xpts, ypts, magnitude, angle);
	vector<Point2f> polarXY, sortedPolarXY;
	vector<double> polarAngles;
	for (int j = 0; j < points[0][0].size(); j++) {
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

	for (int j = 0; j < points[0][0].size(); j++) {
		pp.push_back(Point2f(xnew.at<Float32>(j), ynew.at<Float32>(j)));
	}
	transform(pp.begin(), pp.end(), pp.begin(), std::bind2nd(std::plus<Point2f>(), center));

	this->tumorPoints.push_back(pp);
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


	vector<Point2f> temp;
	for (int i = 1; i < points.size(); i++) {

		vector<double> dist;
		temp = smoothContour(points[i][0], 100, true);

		if (!IsClockwise(temp)) {
			reverse(temp.begin(), temp.end());
		}


		Point2f referencePoint = this->tumorPoints.back().back();
		for (int j = 0; j < temp.size(); j++) {
			dist.push_back(sqrt(pow(referencePoint.x - temp[j].x, 2) + pow(referencePoint.y - temp[j].y, 2)));
		}
		min = std::min_element(dist.begin(), dist.end());
		minIndex = distance(dist.begin(), min);
		std::rotate(temp.begin(), temp.begin() + minIndex, temp.end());

		this->tumorPoints.push_back(temp);
		
	}
	vector<Point2f>().swap(temp);

	vector<vector<vector<Point2f>>>().swap(points);
}

Point2f tumor::getCenterOfMass(vector<Point2f> points) {

	Point2f center = accumulate(points.begin(), points.end(), Point2f(0.0, 0.0));
	center.x /= points.size();
	center.y /= points.size();

	vector<Point2f>().swap(points);
	return center;
}

Point2f tumor::getCenterOfGravity(vector<Point2f> points) {

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

bool tumor::IsClockwise(vector<Point2f> points)
{
	double area = 0;
	int j;
	for (int i = 0; i < points.size(); i++) {
		j = (i + 1) % points.size();
		area += points[i].x * points[j].y;
		area -= points[j].x * points[i].y;
	}
	vector<Point2f>().swap(points);
	return (area / 2) > 0;
}

vector<Point2f> tumor::smoothContour(vector<Point2f> contour, int num_spline, bool closedContour) {

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

// Cubic spline interpolation to smooth centerline
vector<Point2f> tumor::smoothCurve(vector<vector<vector<double>>> centerline, int num_spline) { //vector<vector<double>>
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