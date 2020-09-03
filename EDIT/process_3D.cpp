#include "process_3D.h"


process_3D::process_3D() {
}


process_3D::~process_3D() {
}

vtkSmartPointer<vtkPolyData> process_3D::triangulation(vector<vector<Point2f>> point_cloud, STLType type) {

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
	for (int i = 0; i < point_cloud.size(); i++) {
		for (int j = 0; j < point_cloud[i].size(); j++) {
			points->InsertNextPoint((point_cloud[i][j].x-this->imageCenter.x)*this->xspace,
									(point_cloud[i][j].y-this->imageCenter.y)*this->yspace, 
									i*this->distanceBetweenFrames);
		}
	}

	vector<Point2f> firstContour = point_cloud[0];
	vector<Point2f> lastContour = point_cloud[point_cloud.size()-1];

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

	
	surface_smoothing(trianglePolyData, type, firstContour, lastContour);

	return trianglePolyData;
}


vtkSmartPointer<vtkPolyData> process_3D::surface_smoothing(vtkSmartPointer<vtkPolyData> surface, STLType type, vector<Point2f> firstContour, vector<Point2f> lastContour){
	
	try {
		//smooth surface 
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInputData(surface);
		smoothFilter->SetNumberOfIterations(10); //15
		smoothFilter->SetRelaxationFactor(0.3); //0.3
		smoothFilter->FeatureEdgeSmoothingOff();
		smoothFilter->BoundarySmoothingOn();
		smoothFilter->Update();

		//close holes of surface
		vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter = vtkSmartPointer<vtkFillHolesFilter>::New();
		fillHolesFilter->CAN_PRODUCE_SUB_EXTENT;
		fillHolesFilter->SetHoleSize(10000); //1E6
		fillHolesFilter->SetInputConnection(smoothFilter->GetOutputPort());
		fillHolesFilter->Update();
		/*vtkSmartPointer<vtkPolyDataMapper> filledMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		filledMapper->SetInputData(fillHolesFilter->GetOutput());
		filledMapper->Update();*/

		// Update normals on newly smoothed polydata
		vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalGenerator->SetInputConnection(fillHolesFilter->GetOutputPort());
		normalGenerator->ComputePointNormalsOn();
		normalGenerator->ComputeCellNormalsOn();
		normalGenerator->Update();

		//remesh
		//vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
		////dynamic_cast<vtkLinearSubdivisionFilter*> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(2);
		//subdivisionFilter->SetInputData(normalGenerator->GetOutput());
		//subdivisionFilter->Update();



		process_3D::LoggerMessage("Smoothed model of the surface was produced successfully");

		string filename_stl = this->outputObjectsDir + "/smoothed_" + getSTLName(type) + ".stl";


		//write polydata to .stl file
		vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
		writer->SetFileName(filename_stl.c_str());
		writer->SetInputConnection(normalGenerator->GetOutputPort());
		writer->Write();

		vtkSmartPointer<vtkPolyData> smoothed_branch;

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer< vtkPoints >::New();
		vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();


		smoothed_branch = normalGenerator->GetOutput();


		return smoothed_branch;
	}
	catch (Exception e) {
		LoggerMessage("A problem was occured during the stl smoothing process");
	}

	
}

string process_3D::getSTLName(STLType type) {

	switch (type)
	{
	case process_3D::BLADDER:
		return "Bladder";
		
	case process_3D::SKIN:
		return "Skin";

	case process_3D::THICKNESS:
		return "Thickness";
	}
}