#include "process_3D.h"

#include "messages.h"

process_3D::process_3D() {
}


process_3D::~process_3D() {
}

string process_3D::triangulation(vector<vector<Point2f>> point_cloud, STLType type) {
	if (point_cloud.size() < 2) {
		return warningMessages::cannotProduce3D;
	}

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
	writer->Write();
	if (type == STLType::BLADDER) {
			this->bladderGeometry = filename_stl;
		}
		else if (type == STLType::THICKNESS) {
			this->thicknessGeometry = filename_stl;
		}
	return warningMessages::success;
	*/

	process_3D::LoggerMessage("Initial triangularion was produced successfully");


	return surface_smoothing(trianglePolyData, type);
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