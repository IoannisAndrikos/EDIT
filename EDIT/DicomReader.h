#ifndef DICOMREADER_H
#define DICOMREADER_H
#pragma once
#include"dctk.h"

#include <opencv2/opencv.hpp>
#include <opencv2/cv.h>
#include <opencv2/core.hpp>
#include "opencv2/core/types_c.h"
#include "opencv2/core/utility.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include "dcmtk/dcmimage/diregist.h"

using namespace cv;
using namespace std;


typedef struct _cdicomInfo
{

	double x, y, z;
	double pz;
	double impx, impy, impz;
	char patientName[5000];
	char studyDate[5000];
}CDicomInfo;


/* output file types */
enum E_FileType
{
	EFT_RawPNM,
	EFT_8bitPNM,
	EFT_16bitPNM,
	EFT_NbitPNM,
	EFT_BMP,
	EFT_8bitBMP,
	EFT_24bitBMP,
	EFT_JPEG,
	EFT_TIFF,
	EFT_PNG
};
class CDicomReader
{
public:
	CDicomReader(void);
	~CDicomReader(void);
	CDicomInfo GetDicomInfoD(const wchar_t * inor1, const wchar_t * inor2);
	vector<double> GetDicomInfo(const wchar_t * inor);
	void dcmimage_bmp(const wchar_t * inor, const wchar_t * outor);
	vector<Mat> dcmimage_Mat(const wchar_t * inor, double xmin, double xmax, double ymin, double ymax);
	vector<Mat> dcmimage_Mat(const wchar_t * inor, String outputDir, double xmin, double xmax, double ymin, double ymax);
	int GetInstanceNumber(const wchar_t * inor);
	void hu_conversion(const wchar_t * inor, const wchar_t * outor);
	void DecompressImage(const wchar_t * inor);
	int GetTransferSyntax(const wchar_t * inor);

private:


	int getNumber(char *s1, const char* patientsName, int off)
	{

		int i = off;
		int j = 0;
		do {

			if (patientsName[i] != ' ')
			{
				s1[j] = patientsName[i];
				j++;
			}
			i++;

		} while (i < 100 && patientsName[i] != '\\'&& patientsName[i] != '\0');


		return i;
	}


};
#endif
