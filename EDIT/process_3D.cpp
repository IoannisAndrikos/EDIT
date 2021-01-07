#include "process_3D.h"

#include "messages.h"

process_3D::process_3D() {
}


process_3D::~process_3D() {
}

vector<vector<Point3f>> process_3D::fix3D(vector<vector<Point2f>> point_cloud) {

	double thres = 1.7 * this->distanceBetweenFrames;
	vector<vector<Point3f>> interpolatedPointCloud, finalPointCloud;
	vector<Point3f> frame3D;
	vector<double> difference;
	int interpolationSize;
	double dist, x, y, z;

	String lumenhtxt = "C:/Users/Legion Y540/Desktop/mlk/points.txt";
	ofstream filelumen;
	filelumen.open(lumenhtxt);

	for (int i = 0; i < point_cloud.size(); i++) {
		for (int j = 0; j < point_cloud[i].size(); j++) {
			frame3D.push_back(Point3f((point_cloud[i][j].x - this->imageCenter.x) * this->xspace,
				(point_cloud[i][j].y - this->imageCenter.y) * this->yspace,
				i * this->distanceBetweenFrames));

			filelumen << (point_cloud[i][j].x - this->imageCenter.x) * this->xspace << " " << (point_cloud[i][j].y - this->imageCenter.y) * this->yspace << " " << i * this->distanceBetweenFrames << endl;
		}
		interpolatedPointCloud.push_back(frame3D);
		vector<Point3f>().swap(frame3D);
	}

	for (int i = 0; i < interpolatedPointCloud.size() - 1; i++) {
		finalPointCloud.push_back(interpolatedPointCloud[i]);
		for (int j = 0; j < interpolatedPointCloud[i].size(); j++) {
			difference.push_back(sqrt(pow(interpolatedPointCloud[i][j].x - interpolatedPointCloud[i + 1][j].x, 2) +
				pow(interpolatedPointCloud[i][j].y - interpolatedPointCloud[i + 1][j].y, 2) +
				pow(interpolatedPointCloud[i][j].z - interpolatedPointCloud[i + 1][j].z, 2)));
		}

		//store data those we have interpolated
		auto max = max_element(difference.begin(), difference.end());
		int maxIndex = distance(difference.begin(), max);
		cout << difference[maxIndex] << endl;

		if (difference[maxIndex] > thres) {
			interpolationSize = (int)(difference[maxIndex] / thres);
			for (int k = 1; k <= interpolationSize; k++) {
				for (int j = 0; j < interpolatedPointCloud[i].size(); j++) {
					dist = difference[j] / (interpolationSize + 1);
					cout << "Mphke gia i: " << i << endl;
					x = interpolatedPointCloud[i][j].x + ((k * dist) / difference[j]) * (interpolatedPointCloud[i + 1][j].x - interpolatedPointCloud[i][j].x);
					y = interpolatedPointCloud[i][j].y + ((k * dist) / difference[j]) * (interpolatedPointCloud[i + 1][j].y - interpolatedPointCloud[i][j].y);
					z = interpolatedPointCloud[i][j].z + ((k * dist) / difference[j]) * (interpolatedPointCloud[i + 1][j].z - interpolatedPointCloud[i][j].z);
					frame3D.push_back(Point3f(x, y, z));
					filelumen << x << " " << y << " " << z << endl;

				}
				finalPointCloud.push_back(frame3D);
				vector<Point3f>().swap(frame3D);
			}
		}

		vector<double>().swap(difference);
	}

	filelumen.close();

	finalPointCloud.push_back(interpolatedPointCloud.back());

	return finalPointCloud;
}



string process_3D::triangulation(vector<vector<Point2f>> point_cloud, STLType type) {
	if (point_cloud.size() < 2) {
		return warningMessages::cannotProduce3D;
	}


	vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < point_cloud.size(); i++) {
		for (int j = 0; j < point_cloud[i].size(); j++) {
			points->InsertNextPoint((point_cloud[i][j].x - this->imageCenter.x) * this->xspace,
				(point_cloud[i][j].y - this->imageCenter.y) * this->yspace,
				i * this->distanceBetweenFrames);
		}
	}





	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	

	int ss = point_cloud[0].size();
	int k = ss;
	for (int i = 0; i < point_cloud.size() - 1; i++) {
		vtkSmartPointer<vtkTriangle> triangle_4 = vtkSmartPointer<vtkTriangle>::New();
		triangle_4->GetPointIds()->SetId(0, i * k);
		triangle_4->GetPointIds()->SetId(1, (i + 1) * k - 1);
		triangle_4->GetPointIds()->SetId(2, (i + 1) * k);
		triangles->InsertNextCell(triangle_4);

		for (int j = 0; j < point_cloud[i].size() - 1; j++) {
			vtkSmartPointer<vtkTriangle> triangle_1 = vtkSmartPointer<vtkTriangle>::New();
			triangle_1->GetPointIds()->SetId(2, i * k + j);
			triangle_1->GetPointIds()->SetId(1, i * k + j + 1);
			triangle_1->GetPointIds()->SetId(0, i * k + j + ss);

			vtkSmartPointer<vtkTriangle> triangle_2 = vtkSmartPointer<vtkTriangle>::New();
			triangle_2->GetPointIds()->SetId(0, i * k + j + ss);
			triangle_2->GetPointIds()->SetId(1, i * k + j + ss + 1);
			triangle_2->GetPointIds()->SetId(2, i * k + j + 1);

			triangles->InsertNextCell(triangle_1);
			triangles->InsertNextCell(triangle_2);
		}
		vtkSmartPointer<vtkTriangle> triangle_3 = vtkSmartPointer<vtkTriangle>::New();
		triangle_3->GetPointIds()->SetId(2, i * k + ss - 1);
		triangle_3->GetPointIds()->SetId(1, i * k + ss);
		triangle_3->GetPointIds()->SetId(0, i * k + 2 * ss - 1);
		triangles->InsertNextCell(triangle_3);
	}

	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	poly->SetPolys(triangles);
	poly->SetPoints(points);

	// write the polydata to a file
	//string filename_stl = this->outputObjectsDir + "/no_smoothed_" + getSTLName(type) + ".stl";
	//vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
	//writer->SetFileName(filename_stl.c_str());
	//writer->SetInputData(poly);
	////writer->SetInputData(poly);
	//writer->Write();
	//if (type == STLType::BLADDER) {
	//	this->bladderGeometry = filename_stl;
	//}
	//else if (type == STLType::THICKNESS) {
	//	this->thicknessGeometry = filename_stl;
	//}
	//saveGeometryPath(filename_stl, type);
	//return this->success;

	return surface_smoothing(poly, type);

	process_3D::LoggerMessage("Initial triangularion was produced successfully");

}


string process_3D::surface_smoothing(vtkSmartPointer<vtkPolyData> surface, STLType type){
	
	try {
		
		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputData(surface);
		cleanPolyData->Update();


		//smooth surface 
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInputData(cleanPolyData->GetOutput());
		smoothFilter->SetNumberOfIterations(10); //10
		smoothFilter->SetRelaxationFactor(0.3); //0.3
		smoothFilter->FeatureEdgeSmoothingOff();
		smoothFilter->BoundarySmoothingOn();
		smoothFilter->Update();

		vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter;
		if (this->fillHoles) {
			//close holes of surface
			fillHolesFilter = vtkSmartPointer<vtkFillHolesFilter>::New();
			fillHolesFilter->CAN_PRODUCE_SUB_EXTENT;
			fillHolesFilter->SetHoleSize(10000); //1E6
			fillHolesFilter->SetInputConnection(smoothFilter->GetOutputPort());
			fillHolesFilter->Update();
		}
		
		/*vtkSmartPointer<vtkPolyDataMapper> filledMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		filledMapper->SetInputData(fillHolesFilter->GetOutput());
		filledMapper->Update();*/

		// Update normals on newly smoothed polydata
		vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
		if (this->fillHoles) {
			normalGenerator->SetInputConnection(fillHolesFilter->GetOutputPort());
		}
		else {
			normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
		}
		normalGenerator->FlipNormalsOff();
		normalGenerator->ComputePointNormalsOn();
		normalGenerator->ComputeCellNormalsOn();

		//normalGenerator->get

		normalGenerator->Update();

		

		//remesh
		//vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
		////dynamic_cast<vtkLinearSubdivisionFilter*> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(2);
		//subdivisionFilter->SetInputData(normalGenerator->GetOutput());
		//subdivisionFilter->Update();



		process_3D::LoggerMessage("Smoothed model of the surface was produced successfully");

		string filename_stl = this->outputObjectsDir + separator() + "smoothed_" + getSTLName(type) + ".stl";
			
		//write polydata to .stl file
		vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
		writer->SetFileName(filename_stl.c_str());
		writer->SetInputConnection(normalGenerator->GetOutputPort());
		writer->Write();

		saveGeometryPath(filename_stl, type);

		/*return filename_stl;*/
	}
	catch (Exception e) {
		LoggerMessage("A problem was occured during the stl smoothing process");
		return warningMessages::problemDuringTriangulation;
	}	

	return this->success;
}


string process_3D::findPixelsArePlacedIntoGeometries(vector<vector<vector<Point3f>>> sharderPixels, vector<vector<vector<Point3f>>> interpolatedPixels, STLType type) {

	vector<Point3f> final_Points;

	// bladder
	vtkSmartPointer <vtkSTLReader> readerThickness = vtkSmartPointer<vtkSTLReader>::New();
	readerThickness->SetFileName((this->thicknessGeometry).c_str());
	readerThickness->Update();
	vtkSmartPointer<vtkFillHolesFilter> fillHolesFilterThickness;
	fillHolesFilterThickness = vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHolesFilterThickness->SetHoleSize(10000); //1E6
	fillHolesFilterThickness->SetInputConnection(readerThickness->GetOutputPort());
	fillHolesFilterThickness->Update();
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints_thickness = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selectEnclosedPoints_thickness->Initialize(fillHolesFilterThickness->GetOutput());
	//selectEnclosedPoints_thickness->Update();


	vtkSmartPointer <vtkSTLReader> readerBladder = vtkSmartPointer<vtkSTLReader>::New();
	readerBladder->SetFileName((this->bladderGeometry).c_str());
	readerBladder->Update();
	vtkSmartPointer<vtkFillHolesFilter> fillHolesFilterBladder;
	fillHolesFilterBladder = vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHolesFilterBladder->SetHoleSize(10000); //1E6
	fillHolesFilterBladder->SetInputConnection(readerBladder->GetOutputPort());
	fillHolesFilterBladder->Update();
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints_bladder = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selectEnclosedPoints_bladder->Initialize(fillHolesFilterBladder->GetOutput());


	Point3f point1, point2;
	int con1, con2;

	for (int i = 0; i < sharderPixels.size(); i++) {
		for (int j = 0; j < sharderPixels[i].size(); j++) {
			/*if (!selectEnclosedPoints_bladder->IsInsideSurface(notSharderPixels[i][j].x, notSharderPixels[i][j].y, notSharderPixels[i][j].z) && selectEnclosedPoints_thickness->IsInsideSurface(notSharderPixels[i][j].x, notSharderPixels[i][j].y, notSharderPixels[i][j].z)) {
				final_Points.push_back(notSharderPixels[i][j]);
			}*/
			if (!sharderPixels[i][j].empty()) {
				point1 = sharderPixels[i][j][0];
				point2 = sharderPixels[i][j][1];
				con1 = (!selectEnclosedPoints_bladder->IsInsideSurface(point1.x, point1.y, point1.z) && selectEnclosedPoints_thickness->IsInsideSurface(point1.x, point1.y, point1.z));
				con2 = (!selectEnclosedPoints_bladder->IsInsideSurface(point2.x, point2.y, point2.z) && selectEnclosedPoints_thickness->IsInsideSurface(point2.x, point2.y, point2.z));
				if (con1 && con2) {
					final_Points.insert(final_Points.end() ,interpolatedPixels[i][j].begin(), interpolatedPixels[i][j].end());
					//final_Points.push_back(point1);
					//final_Points.push_back(point2);
				}
			}
		}
	}

	//some free memory
	vector<vector<vector<Point3f>>>().swap(sharderPixels);
	vector<vector<vector<Point3f>>>().swap(interpolatedPixels);


	String txtFilename = this->outputObjectsDir + separator() + getSTLName(type) + ".txt";
	ofstream txtfile;
	txtfile.open(txtFilename);
	for (int i = 0; i < final_Points.size(); i++) {
		txtfile << final_Points[i].x << " " << final_Points[i].y << " " << final_Points[i].z << endl;
	}

	txtfile.close();

	vector<Point3f>().swap(final_Points);
	//vector<Point3f>().swap(skinGeometry);

	saveGeometryPath(txtFilename, type);


	return this->success;
}




string process_3D::writeTumor(vector<vector<vector<Point3f>>> sharderPixels, vector<vector<vector<Point3f>>> interpolatedPixels, STLType type) {

	vector<Point3f> final_Points;

	Point3f point1, point2;
	int con1, con2;

	for (int i = 0; i < sharderPixels.size(); i++) {
		for (int j = 0; j < sharderPixels[i].size(); j++) {

			final_Points.insert(final_Points.end(), interpolatedPixels[i][j].begin(), interpolatedPixels[i][j].end());
			
		}
	}

	cout << "ok 3" << endl;

	//some free memory
	vector<vector<vector<Point3f>>>().swap(sharderPixels);
	vector<vector<vector<Point3f>>>().swap(interpolatedPixels);


	String txtFilename = this->outputObjectsDir + separator() + getSTLName(type) + ".txt";
	ofstream txtfile;
	txtfile.open(txtFilename);
	for (int i = 0; i < final_Points.size(); i++) {
		txtfile << final_Points[i].x << " " << final_Points[i].y << " " << final_Points[i].z << endl;
	}

	txtfile.close();

	vector<Point3f>().swap(final_Points);
	//vector<Point3f>().swap(skinGeometry);

	saveGeometryPath(txtFilename, type);


	return this->success;
}



void process_3D::saveGeometryPath(string  filename, STLType type) {

	switch (type)
	{
	case process_3D::BLADDER:
		this->bladderGeometry = filename;
		break;
	case process_3D::SKIN:
		this->skinGeometry = filename;
		break;
	case process_3D::THICKNESS:
		this->thicknessGeometry = filename;
		break;
	case process_3D::OXY:
		this->OXYGeometry = filename;
		break;
	case process_3D::DeOXY:
		this->DeOXYGeometry = filename;
		break;
	case process_3D::Tumor:
		this->TumorGeometry = filename;
		break;
	default:
		break;
	}
}


string process_3D::getSTLName(STLType type) {

	switch (type)
	{
	case process_3D::BLADDER:
		return "Bladder";
		break;
	case process_3D::SKIN:
		return "Skin";
		break;
	case process_3D::THICKNESS:
		return "Thickness";
		break;
	case process_3D::OXY:
		return "OXY";
		break;
	case process_3D::DeOXY:
		return "DeOXY";
		break;
	case process_3D::Tumor:
		return "Tumor";
		break;

	}
}