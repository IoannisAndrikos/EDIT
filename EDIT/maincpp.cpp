#include "ultrasound.h"
#include "process_3D.h"
#include "photoAcoustic.h"

#include <direct.h>

Point userPoint = { 0,0 };

void initiatePoint(int event, int x, int y, int, void* imgptr) {

	Mat & img = (*(Mat*)imgptr); // first cast, then deref
	Point pt1 = Point(x, y);
	if (event == EVENT_LBUTTONDOWN) {
		userPoint = pt1;
		circle(img, pt1, 2, Scalar(0, 0, 200), 3, 8, 0);
	}
	imshow("window", img);
	waitKey(1);
}

int main(int argv, char* argc) {
	string mainPath = "C:/Users/Legion Y540/Desktop/EDIT_STUDIES";

	string dcmPath = "C:/Users/Legion Y540/Desktop/EDIT Data/month 0/R73/R73_US.dcm";
	//string outputPath = "C:/Users/Legion Y540/Desktop/trial";

	//ultrasound *ultr = new ultrasound(dcmPath, outputPath);
	ultrasound *ultr = new ultrasound();
	ultr->setMainOutputDirectory(mainPath);
	ultr->exportImages(dcmPath);
		

	vector<double> Tags = ultr->getTags();
	vector<Mat> images = ultr->getImages();

	string filename = ultr->getFilename();

	string studyPath = mainPath + "/" + filename;
	_mkdir(studyPath.c_str());



	int iniFrame = 53;
	int lastFrame = 70;


	ultr->repeats = 20;
	ultr->smoothing = 3;
	ultr->lamda1 = 1;
	ultr->lamda2 = 1;
	ultr->levelsetSize = 40;
	ultr->applyEqualizeHist = false;

	Mat inputforPoint = images[iniFrame];
	namedWindow("window", 1);
	setMouseCallback("window", initiatePoint, &inputforPoint);
	imshow("window", inputforPoint);
	waitKey(0);

	ultr->processing(iniFrame, lastFrame, userPoint);


	//vector<vector<Point2f>> lumenPoints = ultr->getlumenPoints();
	//vector<vector<Point2f>> skinPoints = ultr->getSkinPoints();
	ultr->finalizeAllBladderContours(ultr->getlumenPoints());
	vector<vector<Point2f>> lumenPoints = ultr->getlumenPoints();

	//directory for images
	string segmented_images = studyPath + "/segmented_images";
	_mkdir(segmented_images.c_str());
	//directory for points to txt
	string points_txt = studyPath + "/points_txt";
	_mkdir(points_txt.c_str());


	//for (int j = 0; j < lumenPoints.size(); j++) {
	//	Mat img = images[j+iniFrame-1];

	//	
	//	//imwrite(pathbmp, img);


	//	String pathtxt = points_txt + "/" + to_string(j + iniFrame) + ".txt";
	//	ofstream myfilelumen;
	//	myfilelumen.open(pathtxt);



	//	for (int i = 0; i < lumenPoints[j].size(); i++) {
	//		img.at<uchar>(lumenPoints[j][i].y, lumenPoints[j][i].x) = 255;
	//		myfilelumen << lumenPoints[j][i].x << " " << lumenPoints[j][i].y << "\n"; // write points to .txt
	//		//------------------------------------------------------------
	//	}
	//	myfilelumen.close();

	//	String pathbmp = segmented_images + "/" + to_string(j + iniFrame) + ".bmp";
	//	imwrite(pathbmp, img);

	//}
	//0.203

	/*process_3D *proc = new process_3D();


	
	proc->xspace = Tags[0] * 10;
	proc->yspace = Tags[1] * 10;
	proc->distanceBetweenFrames = 0.203;
	proc->setStudyDir(ultr->getStudyDir());

	string stlpath = studyPath + "/stl_objects";
	_mkdir(stlpath.c_str());
	proc->triangulation(lumenPoints, process_3D::STLType::BLADDER);

	

	ultr->extractSkinPoints();

	vector<vector<Point2f>> skinPoints = ultr->getSkinPoints();
	cout << skinPoints.size() << endl;
	cout << skinPoints[0].size() << endl;*/


	photoAcoustic* photo = new photoAcoustic();
	photo->setMainOutputDirectory(mainPath);
	photo->exportOXYImages("C:/Users/Legion Y540/Desktop/EDIT Data/month 0/R73/R73_OXY.dcm");



	photo->setInitialFrame(ultr->getInitialFrame());
	photo->setLastFrame(ultr->getLastFrame());
	photo->setlumenPoints(ultr->getlumenPoints());
	photo->thicknessExtraction();


	//proc->triangulation(skinPoints, process_3D::STLType::SKIN);

	//ultr->writePointsAndImages();
	
}