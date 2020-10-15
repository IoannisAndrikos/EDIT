#include "process_3D.h"



process_3D::process_3D() {
}


process_3D::~process_3D() {
}

string process_3D::triangulation(vector<vector<Point2f>> point_cloud, STLType type) {

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < point_cloud.size(); i++) {
		for (int j = 0; j < point_cloud[i].size(); j++) {
			points->InsertNextPoint((point_cloud[i][j].x-this->imageCenter.x)*this->xspace,
									(point_cloud[i][j].y-this->imageCenter.y)*this->yspace, 
									i*this->distanceBetweenFrames);
		}
	}

	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	int ss = point_cloud[0].size();
	int k = ss;
	for (int i = 0; i < point_cloud.size() - 1; i++) {

		vtkSmartPointer<vtkTriangle> triangle_4 = vtkSmartPointer<vtkTriangle>::New();
		triangle_4->GetPointIds()->SetId(0, i*k);
		triangle_4->GetPointIds()->SetId(1, (i + 1)*k - 1);
		triangle_4->GetPointIds()->SetId(2, (i + 1)*k);

		triangles->InsertNextCell(triangle_4);

		for (int j = 0; j < point_cloud[i].size() - 1; j++) {
			vtkSmartPointer<vtkTriangle> triangle_1 = vtkSmartPointer<vtkTriangle>::New();
			triangle_1->GetPointIds()->SetId(2, i*k + j);
			triangle_1->GetPointIds()->SetId(1, i*k + j + 1);
			triangle_1->GetPointIds()->SetId(0, i*k + j + ss);

			triangles->InsertNextCell(triangle_1);

			vtkSmartPointer<vtkTriangle> triangle_2 = vtkSmartPointer<vtkTriangle>::New();
			triangle_2->GetPointIds()->SetId(0, i*k + j + ss);
			triangle_2->GetPointIds()->SetId(1, i*k + j + ss + 1);
			triangle_2->GetPointIds()->SetId(2, i*k + j + 1);
			triangles->InsertNextCell(triangle_2);
		}

		vtkSmartPointer<vtkTriangle> triangle_3 = vtkSmartPointer<vtkTriangle>::New();
		triangle_3->GetPointIds()->SetId(2, i*k + ss - 1);
		triangle_3->GetPointIds()->SetId(1, i*k + ss);
		triangle_3->GetPointIds()->SetId(0, i*k + 2 * ss - 1);

		triangles->InsertNextCell(triangle_3);

	}

	// Create a polydata object
	vtkSmartPointer<vtkPolyData> trianglePolyData = vtkSmartPointer<vtkPolyData>::New();

	// Add the geometry and topology to the polydata
	trianglePolyData->SetPoints(points);
	trianglePolyData->SetPolys(triangles);


	// write the polydata to a file
	/*string filename_stl = this->outputObjectsDir + "/no_smoothed_" + getSTLName(type) + ".stl";
	vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
	writer->SetFileName(filename_stl.c_str());
	writer->SetInputData(trianglePolyData);
	writer->Write();*/

	process_3D::LoggerMessage("Initial triangularion was produced successfully");

	
	string filename_stl = surface_smoothing(trianglePolyData, type);

	return filename_stl;
}


string process_3D::surface_smoothing(vtkSmartPointer<vtkPolyData> surface, STLType type){
	
	try {
		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputData(surface);
		cleanPolyData->Update();


		//smooth surface 
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInputData(cleanPolyData->GetOutput());
		smoothFilter->SetNumberOfIterations(10); //15
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

		if (type == STLType::BLADDER) {
			this->bladderGeometry = filename_stl;
		}
		else if (type == STLType::THICKNESS) {
			this->thicknessGeometry = filename_stl;
		}
			

		//write polydata to .stl file
		vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
		writer->SetFileName(filename_stl.c_str());
		writer->SetInputConnection(normalGenerator->GetOutputPort());
		writer->Write();

		vtkSmartPointer<vtkPolyData> smoothed_branch;

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
		vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();


		smoothed_branch = normalGenerator->GetOutput();


		return filename_stl;
	}
	catch (Exception e) {
		LoggerMessage("A problem was occured during the stl smoothing process");
	}	
}


string process_3D::findPixelsArePlacedIntoGeometries(vector<vector<Point3f>> pixelsFrames, STLType type) {
	
	vector<vector<Point3f>> final_Points;
	vector<int> pointSizePerFrame;

	vector<Point3f> pixels3D;
	for (int i = 0; i < pixelsFrames.size(); i++) {
		pixels3D.insert(pixels3D.end(), pixelsFrames[i].begin(), pixelsFrames[i].end());
		pointSizePerFrame.push_back(pixelsFrames[i].size());
	}

	// bladder
	vtkSmartPointer <vtkSTLReader> readerThickness = vtkSmartPointer<vtkSTLReader>::New();
	readerThickness->SetFileName((this->thicknessGeometry).c_str());
	readerThickness->Update();
	vtkSmartPointer<vtkFillHolesFilter> fillHolesFilterThickness;
	fillHolesFilterThickness = vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHolesFilterThickness->SetHoleSize(10000); //1E6
	fillHolesFilterThickness->SetInputConnection(readerThickness->GetOutputPort());
	fillHolesFilterThickness->Update();
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints_thickness  = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
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
	//selectEnclosedPoints_bladder->Update();

	//cout << "pointSizePerFrame Size: " << pointSizePerFrame.size() << endl;
	//
	//vector<vector<Point3f>> allframes;
	//vector<Point3f> frame;

	//int count = 0;
	//int k = 0;
	//for (int i = 0; i < pixels3D.size(); i++) {

	////	cout << "k value: " << k << endl;

	//	if(!selectEnclosedPoints_bladder->IsInsideSurface(pixels3D[i].x, pixels3D[i].y, pixels3D[i].z) && selectEnclosedPoints_thickness->IsInsideSurface(pixels3D[i].x, pixels3D[i].y, pixels3D[i].z)){
	//		final_Points.push_back(pixels3D[i]);
	//		frame.push_back(pixels3D[i]);
	//	}
	//	if (count == pointSizePerFrame[k]) {
	//		allframes.push_back(frame);
	//		vector<Point3f>().swap(frame);
	//		count = 0;
	//		k++;
	//	}
	//	count++;
	//}


	vector<vector<Point3f>> final_FramePoints;
	vector<Point3f> FramePoints;
	for (int i = 0; i < pixelsFrames.size(); i++) {
		for (int j = 0; j < pixelsFrames[i].size(); j++) {
			if (!selectEnclosedPoints_bladder->IsInsideSurface(pixelsFrames[i][j].x, pixelsFrames[i][j].y, pixelsFrames[i][j].z) && selectEnclosedPoints_thickness->IsInsideSurface(pixelsFrames[i][j].x, pixelsFrames[i][j].y, pixelsFrames[i][j].z)) {
				FramePoints.push_back(pixelsFrames[i][j]);
			}
		}
		final_Points.push_back(FramePoints);
		final_FramePoints.push_back(FramePoints);
		vector<Point3f>().swap(FramePoints);
	}


	cout << "final_FramePoints: " << final_FramePoints.size() << endl;



	double precision = 0.0000001;

	double threshold = 0.0500000001111111111111;
	//vector<Point3f> interpolated_Points(final_Points.begin(), final_Points.end());
	
	//vector<Point3f> extrapoints;
	//double x1, y1, z1;
	//double x2, y2, z2;
	//for (int i = 0; i < final_FramePoints.size()-1; i++) {
	//	for (int j = 0; j < final_FramePoints[i].size(); j++) {
	//		x1 = final_FramePoints[i][j].x;
	//		y1 = final_FramePoints[i][j].y;
	//		z1 = final_FramePoints[i][j].z;

	//		x2 = final_FramePoints[i + 1][j].x;
	//		y2 = final_FramePoints[i + 1][j].y;
	//		z2 = final_FramePoints[i + 1][j].z;

	//		/*double dist = sqrt(pow(x2 - x1, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));*/
	//		//if ((dist+ precision )< threshold) {
	//			extrapoints.push_back(Point3f(x1, y1, (z1 + z2) / 2));
	//		//}
	//		

	//		/*cout <<"ith " << final_FramePoints[i][j].x << " " << final_FramePoints[i][j].y << " " << final_FramePoints[i][j].z << endl;
	//		cout << "ith+1 "<< final_FramePoints[i+1][j].x << " " << final_FramePoints[i+1][j].y << " " << final_FramePoints[i+1][j].z << endl;*/

	//		//cout << x << " " << y << " " << z << endl;

	//		//if (x == 0 && y == 0) {
	//		//	extrapoints.push_back(Point3f(final_FramePoints[i][j].x, final_FramePoints[i][j].y, (final_FramePoints[i][j].z + final_FramePoints[i+1][j].z) / 2));
	//		//	//cout << "MPHKE" << endl;
	//		//}
	//		////cout << "EDW" << i << j << endl;

	//	}
	//}
	//final_Points.push_back(extrapoints);
	
	String txtFilename = this->outputObjectsDir + separator() + getSTLName(type) + ".txt";
	ofstream txtfile;
	txtfile.open(txtFilename);
	for (int i = 0; i < final_Points.size(); i++) {
		for (int j = 0; j < final_Points[i].size(); j++) {
			txtfile << final_Points[i][j].x << " " << final_Points[i][j].y << " " << final_Points[i][j].z << endl;
		}
	}

	/*for (int k = 0; k < interpolated_Points.size(); k++) {
		txtfile << interpolated_Points[k].x << " " << interpolated_Points[k].y << " " << interpolated_Points[k].z << endl;
	}*/

	txtfile.close();

	vector<vector<Point3f>>().swap(final_Points);
	//vector<Point3f>().swap(interpolated_Points);

	return txtFilename;
}


string process_3D::findPixelsArePlacedIntoGeometries2(vector<vector<vector<Point3f>>> sharderPixels, vector<vector<vector<Point3f>>> interpolatedPixels, STLType type) {


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
	//vector<Point3f>().swap(interpolated_Points);

	return txtFilename;
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

	}
}