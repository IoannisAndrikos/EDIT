#include "StdAfx.h"
#include "DicomReader.h"

CDicomReader::CDicomReader(void)
{
}

CDicomReader::~CDicomReader(void)
{
}



void CDicomReader::dcmimage_bmp(const wchar_t * inor, const wchar_t * outor)
{

	E_FileReadMode      opt_readMode = ERM_autoDetect;    /* default: fileformat or dataset */
	E_TransferSyntax    opt_transferSyntax = EXS_Unknown; /* default: xfer syntax recognition */
	unsigned long       opt_compatibilityMode = CIF_MayDetachPixelData | CIF_TakeOverExternalDataset;
	/* default: pixel data may detached if no longer needed */
	OFCmdUnsignedInt    opt_frame = 1;                    /* default: first frame */
	OFCmdUnsignedInt    opt_frameCount = 0;               /* default: one frame */
	int                 opt_multiFrame = 1;               /* default: no multiframes */
	int                 opt_convertToGrayscale = 0;       /* default: color or grayscale */

	OFCmdUnsignedInt    opt_quality = 90;                 /* default: 90% JPEG quality */
	E_SubSampling       opt_sampling = ESS_422;           /* default: 4:2:2 sub-sampling */
	E_DecompressionColorSpaceConversion opt_decompCSconversion = EDC_photometricInterpretation;

	int                 opt_verboseMode = 1;              /* default: be more or less quiet */
	int                 opt_debugMode = 0;                /* default: no debug */
	int                 opt_suppressOutput = 0;           /* default: create output */
	E_FileType          opt_fileType = EFT_RawPNM;        /* default: 8-bit PGM/PPM */
														  /* (binary for file output and ASCII for stdout) */
	OFCmdUnsignedInt    opt_fileBits = 0;                 /* default: 0 */
	const char *        opt_ifname = NULL;
	const char *        opt_ofname = NULL;


	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char in[newsize];
	char out[newsize];
	wcstombs_s(&convertedChars, in, origsize, inor, _TRUNCATE);

	origsize = wcslen(outor) + 1;

	convertedChars = 0;

	wcstombs_s(&convertedChars, out, origsize, outor, _TRUNCATE);
	opt_ifname = in;
	opt_ofname = out;
	opt_readMode = ERM_autoDetect;
	opt_fileType = EFT_BMP;

	DcmRLEDecoderRegistration::registerCodecs(OFFalse /*pCreateSOPInstanceUID*/, (opt_debugMode ? OFTrue : OFFalse));
	DJDecoderRegistration::registerCodecs(opt_decompCSconversion, EUC_default, EPC_default, (opt_debugMode ? OFTrue : OFFalse));

	DcmFileFormat *dfile = new DcmFileFormat();
	OFCondition cond = dfile->loadFile(opt_ifname, opt_transferSyntax, EGL_withoutGL, DCM_MaxReadLength, opt_readMode);

	E_TransferSyntax xfer = dfile->getDataset()->getOriginalXfer();

	DicomImage *di = new DicomImage(dfile, xfer, opt_compatibilityMode, opt_frame - 1, opt_frameCount);

	if (!opt_suppressOutput)
	{
		if ((opt_convertToGrayscale) && (!di->isMonochrome()))
		{
			//    if (opt_verboseMode > 1)
			//      CERR << "converting image to grayscale." << endl;

			DicomImage *newimage = di->createMonochromeImage();
			if (newimage == NULL)
			{

			}
			else if (newimage->getStatus() != EIS_Normal)
			{

			}
			else
			{
				delete di;
				di = newimage;
			}
		}

		int result = 0;
		FILE *ofile = NULL;
		char ofname[255];
		//unsigned int fcount = OFstatic_cast(unsigned int, ((opt_frameCount > 0) && (opt_frameCount <= di->getFrameCount())) ? opt_frameCount : di->getFrameCount());
		unsigned int fcount = di->getFrameCount();
		const char *ofext = NULL;

		switch (opt_fileType)
		{
		case EFT_BMP:
		case EFT_8bitBMP:
		case EFT_24bitBMP:
			ofext = "bmp";
			break;
		case EFT_JPEG:
			ofext = "jpg";
			break;
		case EFT_TIFF:
			ofext = "tif";
			break;
		case EFT_PNG:
			ofext = "png";
			break;
		default:
			if (di->isMonochrome()) ofext = "pgm"; else ofext = "ppm";
			break;
		}

		for (unsigned int frame = 0; frame < fcount; frame++)
		{
			if (opt_ofname)
			{
				if (opt_multiFrame)
					sprintf(ofname, "%s%d.%s", opt_ofname, frame, ofext);
				else
					strcpy(ofname, opt_ofname);
				//        if (opt_verboseMode > 1)
				//                     CERR << "writing frame " << (opt_frame + frame) << " to " << ofname << endl;
				ofile = fopen(ofname, "wb");
				if (ofile == NULL)
				{
					OFOStringStream oss;
					oss << "cannot create file " << ofname << OFStringStream_ends;
					OFSTRINGSTREAM_GETSTR(oss, tmpString)
						OFSTRINGSTREAM_FREESTR(tmpString)
				}
			}
			else {
				ofile = stdout;
				//      if (opt_verboseMode > 1)
				//           CERR << "writing frame " << (opt_frame + frame) << " to stdout" << endl;
			}

			switch (opt_fileType)
			{
			case EFT_RawPNM:
				result = di->writeRawPPM(ofile, 8, frame);
				break;
			case EFT_8bitPNM:
				result = di->writePPM(ofile, 8, frame);
				break;
			case EFT_16bitPNM:
				result = di->writePPM(ofile, 16, frame);
				break;
			case EFT_NbitPNM:
				result = di->writePPM(ofile, OFstatic_cast(int, opt_fileBits), frame);
				break;
			case EFT_BMP:
				di->setMinMaxWindow(0);
				result = di->writeBMP(ofile, 0, frame);
				break;
			case EFT_8bitBMP:
				result = di->writeBMP(ofile, 8, frame);
				break;
			case EFT_24bitBMP:
				result = di->writeBMP(ofile, 24, frame);
				break;

			case EFT_JPEG:
			{
				DiJPEGPlugin plugin;
				plugin.setQuality(OFstatic_cast(unsigned int, opt_quality));
				plugin.setSampling(opt_sampling);
				result = di->writePluginFormat(&plugin, ofile, frame);
			}
			break;

			default:
				if (opt_ofname)
					result = di->writeRawPPM(ofile, 8, frame);
				else
					result = di->writePPM(ofile, 8, frame);
				break;
			}

			if (opt_ofname)
				fclose(ofile);

			if (!result)
			{

			}
		}
	}
	//
	//if (opt_verboseMode > 1)
	//   CERR << "cleaning up memory." << endl;
	delete di;

	DcmRLEDecoderRegistration::cleanup();
	DJDecoderRegistration::cleanup();

}



vector<Mat> CDicomReader::dcmimage_Mat(const wchar_t * inor, double xmin, double xmax, double ymin, double ymax)
{

	E_FileReadMode      opt_readMode = ERM_autoDetect;    /* default: fileformat or dataset */
	E_TransferSyntax    opt_transferSyntax = EXS_Unknown; /* default: xfer syntax recognition */
	unsigned long       opt_compatibilityMode = CIF_MayDetachPixelData | CIF_TakeOverExternalDataset;
	/* default: pixel data may detached if no longer needed */
	OFCmdUnsignedInt    opt_frame = 1;                    /* default: first frame */
	OFCmdUnsignedInt    opt_frameCount = 0;               /* default: one frame */
	int                 opt_multiFrame = 1;               /* default: no multiframes */
	int                 opt_convertToGrayscale = 0;       /* default: color or grayscale */

	OFCmdUnsignedInt    opt_quality = 90;                 /* default: 90% JPEG quality */
	E_SubSampling       opt_sampling = ESS_422;           /* default: 4:2:2 sub-sampling */
	E_DecompressionColorSpaceConversion opt_decompCSconversion = EDC_photometricInterpretation;

	int                 opt_verboseMode = 1;              /* default: be more or less quiet */
	int                 opt_debugMode = 0;                /* default: no debug */
	int                 opt_suppressOutput = 0;           /* default: create output */
	E_FileType          opt_fileType = EFT_RawPNM;        /* default: 8-bit PGM/PPM */
														  /* (binary for file output and ASCII for stdout) */
	OFCmdUnsignedInt    opt_fileBits = 0;                 /* default: 0 */
	const char *        opt_ifname = NULL;
	const char *        opt_ofname = NULL;


	vector<Mat> images;

	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char in[newsize];
	char out[newsize];
	wcstombs_s(&convertedChars, in, origsize, inor, _TRUNCATE);


	convertedChars = 0;

	opt_ifname = in;

	DcmRLEDecoderRegistration::registerCodecs(OFFalse /*pCreateSOPInstanceUID*/, (opt_debugMode ? OFTrue : OFFalse));
	DJDecoderRegistration::registerCodecs(opt_decompCSconversion, EUC_default, EPC_default, (opt_debugMode ? OFTrue : OFFalse));

	DcmFileFormat *dfile = new DcmFileFormat();
	OFCondition cond = dfile->loadFile(opt_ifname, opt_transferSyntax, EGL_withoutGL, DCM_MaxReadLength, opt_readMode);

	E_TransferSyntax xfer = dfile->getDataset()->getOriginalXfer();

	DicomImage *di = new DicomImage(dfile, xfer, opt_compatibilityMode, opt_frame - 1, opt_frameCount);

	if (!opt_suppressOutput)
	{
		if ((opt_convertToGrayscale) && (!di->isMonochrome()))
		{
			//    if (opt_verboseMode > 1)
			//      CERR << "converting image to grayscale." << endl;

			DicomImage *newimage = di->createMonochromeImage();
			if (newimage == NULL)
			{

			}
			else if (newimage->getStatus() != EIS_Normal)
			{

			}
			else
			{
				delete di;
				di = newimage;
			}
		}


		for (unsigned int frame = 0; frame < di->getFrameCount(); frame++)
		{
			di->setMinMaxWindow(0);
			cv::Mat cvimage(int(di->getHeight()), int(di->getWidth()), CV_MAKETYPE(di->getDepth(), 3), (Uint8*)di->getOutputData(di->getDepth(), frame, 0));
			cvtColor(cvimage, cvimage, COLOR_BGR2GRAY); //COLOR_BGR2RGB COLOR_BGR2GRAY

			//some cropping

			int x = round(xmax) - round(xmin);
			int y = round(ymax) - round(ymin);
			cvimage = cvimage(Rect(round(xmin), round(ymin), x, y));

			//imwrite

			images.push_back(cvimage);
		}

	}
	delete di;

	DcmRLEDecoderRegistration::cleanup();
	DJDecoderRegistration::cleanup();

	return images;

}

//third one 
vector<Mat> CDicomReader::dcmimage_Mat(const wchar_t * inor, String outputDir, double xmin, double xmax, double ymin, double ymax, ImageChannel channel)
{

	E_FileReadMode      opt_readMode = ERM_autoDetect;    /* default: fileformat or dataset */
	E_TransferSyntax    opt_transferSyntax = EXS_Unknown; /* default: xfer syntax recognition */
	unsigned long       opt_compatibilityMode = CIF_MayDetachPixelData | CIF_TakeOverExternalDataset;
	/* default: pixel data may detached if no longer needed */
	OFCmdUnsignedInt    opt_frame = 1;                    /* default: first frame */
	OFCmdUnsignedInt    opt_frameCount = 0;               /* default: one frame */
	int                 opt_multiFrame = 1;               /* default: no multiframes */
	int                 opt_convertToGrayscale = 0;       /* default: color or grayscale */

	OFCmdUnsignedInt    opt_quality = 90;                 /* default: 90% JPEG quality */
	E_SubSampling       opt_sampling = ESS_422;           /* default: 4:2:2 sub-sampling */
	E_DecompressionColorSpaceConversion opt_decompCSconversion = EDC_photometricInterpretation;

	int                 opt_verboseMode = 1;              /* default: be more or less quiet */
	int                 opt_debugMode = 0;                /* default: no debug */
	int                 opt_suppressOutput = 0;           /* default: create output */
	E_FileType          opt_fileType = EFT_RawPNM;        /* default: 8-bit PGM/PPM */
														  /* (binary for file output and ASCII for stdout) */
	OFCmdUnsignedInt    opt_fileBits = 0;                 /* default: 0 */
	const char *        opt_ifname = NULL;
	const char *        opt_ofname = NULL;


	vector<Mat> images;

	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char in[newsize];
	char out[newsize];
	wcstombs_s(&convertedChars, in, origsize, inor, _TRUNCATE);


	convertedChars = 0;

	opt_ifname = in;

	DcmRLEDecoderRegistration::registerCodecs(OFFalse /*pCreateSOPInstanceUID*/, (opt_debugMode ? OFTrue : OFFalse));
	DJDecoderRegistration::registerCodecs(opt_decompCSconversion, EUC_default, EPC_default, (opt_debugMode ? OFTrue : OFFalse));

	DcmFileFormat *dfile = new DcmFileFormat();
	OFCondition cond = dfile->loadFile(opt_ifname, opt_transferSyntax, EGL_withoutGL, DCM_MaxReadLength, opt_readMode);

	E_TransferSyntax xfer = dfile->getDataset()->getOriginalXfer();

	DicomImage *di = new DicomImage(dfile, xfer, opt_compatibilityMode, opt_frame - 1, opt_frameCount);

	if (!opt_suppressOutput)
	{
		if ((opt_convertToGrayscale) && (!di->isMonochrome()))
		{
			//    if (opt_verboseMode > 1)
			//      CERR << "converting image to grayscale." << endl;

			DicomImage *newimage = di->createMonochromeImage();
			if (newimage == NULL)
			{

			}
			else if (newimage->getStatus() != EIS_Normal)
			{

			}
			else
			{
				delete di;
				di = newimage;
			}
		}


		for (unsigned int frame = 0; frame < di->getFrameCount(); frame++)
		{
			di->setMinMaxWindow(0);
			cv::Mat cvimage(int(di->getHeight()), int(di->getWidth()), CV_MAKETYPE(di->getDepth(), 3), (Uint8*)di->getOutputData(di->getDepth(), frame, 0));
			cvtColor(cvimage, cvimage, COLOR_BGR2RGB); //COLOR_BGR2RGB COLOR_BGR2GRAY

			//some cropping

			int x = round(xmax) - round(xmin);
			int y = round(ymax) - round(ymin);
			cvimage = cvimage(Rect(round(xmin), round(ymin), x, y));
			String pathbmp = outputDir + "/" + to_string(frame) + ".bmp";
			imwrite(pathbmp, cvimage);
			
			switch (channel)
			{
			case CDicomReader::GRAYSCALE:
				cvtColor(cvimage, cvimage, cv::COLOR_RGB2GRAY);
				images.push_back(cvimage);
				break;
			case CDicomReader::RED:
				Mat planes[3];
				split(cvimage, planes);
				images.push_back(planes[2]);
				planes->release();
				break;
			/*case CDicomReader::GREEN:
				break;
			case CDicomReader::BLUE:
				break;
			default:
				break;*/
			}
			
			cvimage.release();
		}

	}
	delete di;

	DcmRLEDecoderRegistration::cleanup();
	DJDecoderRegistration::cleanup();

	return images;

}






vector<double> CDicomReader::GetDicomInfo(const wchar_t * inor)
{

	vector<double> Tags;


	//int instanceNumber=-1;
	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char dicomname[newsize];

	wcstombs_s(&convertedChars, dicomname, origsize, inor, _TRUNCATE);

	DcmFileFormat *fileformat = new DcmFileFormat();
	OFString patientsName;

	DcmItem *ditem = NULL;


	//OFCondition status = fileformat->loadFile("dicom\\IM-0001-0001.dcm");  findAndGetUint8Array
	OFCondition status = fileformat->loadFile(dicomname);

	if (status.good())
	{
		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if(ditem->findAndGetOFString(DCM_PhysicalDeltaX, patientsName).good())
			Tags.push_back(atof(patientsName.data()));
		}


		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if (ditem->findAndGetOFString(DCM_PhysicalDeltaY, patientsName).good())
				Tags.push_back(atof(patientsName.data()));
		}

		//starting pixel X (for crop)
		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if (ditem->findAndGetOFString(DCM_RegionLocationMinX0, patientsName).good())
				Tags.push_back(atof(patientsName.data()));
		}

		//ending pixel X (for crop) 
		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if (ditem->findAndGetOFString(DCM_RegionLocationMaxX1, patientsName).good())
				Tags.push_back(atof(patientsName.data()));
		}

		//starting pixel Y (for crop) 
		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if (ditem->findAndGetOFString(DCM_RegionLocationMinY0, patientsName).good())
				Tags.push_back(atof(patientsName.data()));
		}

		//ending pixel Y (for crop) 
		if (fileformat->getDataset()->findAndGetSequenceItem(DCM_SequenceOfUltrasoundRegions, ditem).good()) {
			if (ditem->findAndGetOFString(DCM_RegionLocationMaxY1, patientsName).good())
				Tags.push_back(atof(patientsName.data()));
		}
	}

	return Tags;
}


CDicomInfo CDicomReader::GetDicomInfoD(const wchar_t * inor1, const wchar_t * inor2)
{

	double z1;
	double z2;

	CDicomInfo inf;
	//int instanceNumber=-1;
	size_t origsize1 = wcslen(inor1) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char dicomname1[newsize];
	char dicomname2[newsize];

	wcstombs_s(&convertedChars, dicomname1, origsize1, inor1, _TRUNCATE);


	size_t origsize2 = wcslen(inor2) + 1;
	const size_t newsize2 = 500;
	size_t convertedChars2 = 0;
	wcstombs_s(&convertedChars2, dicomname2, origsize2, inor2, _TRUNCATE);


	DcmFileFormat *fileformat1 = new DcmFileFormat();
	DcmFileFormat *fileformat2 = new DcmFileFormat();
	OFString patientsName;


	//OFCondition status = fileformat->loadFile("dicom\\IM-0001-0001.dcm");
	OFCondition status1 = fileformat1->loadFile(dicomname1);
	OFCondition status2 = fileformat2->loadFile(dicomname2);

	if (status1.good())
	{



		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0028, 0x0030), patientsName).good())
		{


			char s1[100];
			char s2[100];
			const char *pixels = patientsName.data();


			inf.x = atof(pixels);
			inf.y = atof(pixels);






		}
		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0010, 0x0010), patientsName).good())
		{


			char s1[100];
			char s2[100];
			const char *pixels = patientsName.data();


			strcpy(inf.patientName, pixels);






		}

		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0008, 0x0030), patientsName).good())
		{


			char s1[100];
			char s2[100];
			const char *pixels = patientsName.data();


			strcpy(inf.studyDate, pixels);






		}

		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x1041), patientsName).good())
		{

			inf.pz = atof(patientsName.data());


		}




		//(0020, 1041)	SliceLocation - 1395.9
		//DcmFileFormat *fileformat1 = new DcmFileFormat();
		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x1041), patientsName).good())
		{


			z1 = atof(patientsName.data());





		}


		//(0020, 0032)	ImagePositionPatient - 68.83203125\ - 263.83203125\1395.875
		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x0032), patientsName, 0).good())
		{

			char s1[100];
			int i = 0;
			inf.impx = atof(patientsName.data());
		}
		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x0032), patientsName, 1).good())
		{

			char s1[100];
			int i = 0;
			inf.impy = atof(patientsName.data());
		}

		if (fileformat1->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x0032), patientsName, 2).good())
		{

			char s1[100];
			int i = 0;
			inf.impz = atof(patientsName.data());
		}
	}


	if (status2.good())
	{
		//(0020, 1041)	SliceLocation - 1395.9
		//DcmFileFormat *fileformat1 = new DcmFileFormat();
		if (fileformat2->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x1041), patientsName).good())
		{


			z2 = atof(patientsName.data());

			inf.pz = z2 - z1;



		}


	}


	return inf;

}







int CDicomReader::GetInstanceNumber(const wchar_t * inor)
{


	int instanceNumber = -1;
	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char dicomname[newsize];

	wcstombs_s(&convertedChars, dicomname, origsize, inor, _TRUNCATE);

	DcmFileFormat *fileformat = new DcmFileFormat();
	OFString patientsName;


	//OFCondition status = fileformat->loadFile("dicom\\IM-0001-0001.dcm");
	OFCondition status = fileformat->loadFile(dicomname);

	if (status.good())
	{
		if (fileformat->getDataset()->findAndGetOFString(DcmTagKey(0x0020, 0x0013), patientsName).good())
		{

			instanceNumber = atoi(patientsName.data());
		}


	}


	return instanceNumber;

}





void CDicomReader::hu_conversion(const wchar_t * inor, const wchar_t * outor)
{



	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;


	char in[newsize];
	char out[newsize];
	wcstombs_s(&convertedChars, in, origsize, inor, _TRUNCATE);

	origsize = wcslen(outor) + 1;

	convertedChars = 0;

	wcstombs_s(&convertedChars, out, origsize, outor, _TRUNCATE);






	E_FileReadMode      opt_readMode = ERM_autoDetect;    /* default: fileformat or dataset */
	E_TransferSyntax    opt_transferSyntax = EXS_Unknown; /* default: xfer syntax recognition */
														  //unsigned long       opt_compatibilityMode = CIF_MayDetachPixelData | CIF_TakeOverExternalDataset;
														  //	unsigned long       opt_compatibilityMode = CIF_MayDetachPixelData | CIF_TakeOverExternalDataset | CIF_WrongPaletteAttributeTags;

	unsigned long       opt_compatibilityMode = CIF_AcrNemaCompatibility;
	/* default: pixel data may detached if no longer needed */
	OFCmdUnsignedInt    opt_frame = 1;                    /* default: first frame */
	OFCmdUnsignedInt    opt_frameCount = 1;               /* default: one frame */
														  //int                 opt_multiFrame = 0;               /* default: no multiframes */
														  //int                 opt_convertToGrayscale = 0;       /* default: color or grayscale */

														  //OFCmdUnsignedInt    opt_quality = 90;                 /* default: 90% JPEG quality */
														  //E_SubSampling       opt_sampling = ESS_422;           /* default: 4:2:2 sub-sampling */
	E_DecompressionColorSpaceConversion opt_decompCSconversion = EDC_photometricInterpretation;

	//int                 opt_verboseMode = 1;              /* default: be more or less quiet */
	int                 opt_debugMode = 0;                /* default: no debug */
														  //int                 opt_suppressOutput = 0;           /* default: create output */
	E_FileType          opt_fileType = EFT_RawPNM;        /* default: 8-bit PGM/PPM */
														  //													  /* (binary for file output and ASCII for stdout) */
														  //OFCmdUnsignedInt    opt_fileBits = 0;                 /* default: 0 */
														  //const char *        opt_ifname = NULL;
														  //const char *        opt_ofname = NULL;

														  //opt_ifname = "C:\\ARTreatData\\CT_Carotid\\BA1\\dicom\\1.dcm";
														  //opt_ofname = out;
	opt_readMode = ERM_autoDetect;
	opt_fileType = EFT_BMP;

	DcmRLEDecoderRegistration::registerCodecs(OFFalse /*pCreateSOPInstanceUID*/, (opt_debugMode ? OFTrue : OFFalse));
	DJDecoderRegistration::registerCodecs(opt_decompCSconversion, EUC_default, EPC_default, (opt_debugMode ? OFTrue : OFFalse));

	DcmFileFormat *dfile = new DcmFileFormat();
	OFCondition cond = dfile->loadFile(in, opt_transferSyntax, EGL_withoutGL, DCM_MaxReadLength, opt_readMode);

	E_TransferSyntax xfer = dfile->getDataset()->getOriginalXfer();

	DicomImage *di = new DicomImage(dfile, xfer, opt_compatibilityMode, opt_frame - 1, opt_frameCount);

	di->setMinMaxWindow(0);

	int iPos = 0;
	OFString str;

	DiPixel* dmp = const_cast<DiPixel*>(di->getInterData());
	void *pixelData = NULL;
	pixelData = (void *)dmp->getData();
	Sint16* p = (Sint16*)pixelData;

	//   dcm.getDataset()->findAndGetUint16Array(DCM_PixelData, p); 

	Sint16 tmp[512 * 512];

	//FILE *fp=fopen("hu.bin","wb");


	// ofstream myf;
	//myf.open("hu.bin",ios::out | ios::binary);

	int count = 0;
	for (int i = 0; i < 512; i++)
	{

		for (int j = 0; j < 512; j++)
		{
			Sint16 v = *p;
			if (v>60000)
				v = -3024;
			else
			{


			}
			// if ((i > 55) && (i < 456) && (j > 55) && (j < 456))
			{
				tmp[count++] = v;
			}
			p++;
		}

	}
	FILE *fp = fopen(out, "wb");
	//myf.write((char *) tmp, sizeof(Uint16)*512*512);
	fwrite(tmp, sizeof(Sint16), 512 * 512, fp);

	fclose(fp);

}

void CDicomReader::DecompressImage(const wchar_t * inor)
{

	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;

	char dicomname[newsize];
	wcstombs_s(&convertedChars, dicomname, origsize, inor, _TRUNCATE);

	DcmFileFormat *fileformat = new DcmFileFormat();
	OFString patientsName;

	const int ROWS = 512; const int COLUMNS = 512;
	const int SAMPLES_PER_PIXEL = 1;
	int frameSize = ROWS * COLUMNS * SAMPLES_PER_PIXEL;

	E_FileReadMode opt_readMode = ERM_autoDetect;    /* default: fileformat or dataset */
	E_TransferSyntax opt_transferSyntax = EXS_Unknown;

	unsigned long opt_compatibilityMode = CIF_AcrNemaCompatibility;
	/* default: pixel data may detached if no longer needed */
	OFCmdUnsignedInt opt_frame = 1;                    /* default: first frame */
	OFCmdUnsignedInt opt_frameCount = 1;
	E_DecompressionColorSpaceConversion opt_decompCSconversion = EDC_photometricInterpretation;

	int opt_debugMode = 0;
	E_FileType opt_fileType = EFT_BMP;

	DcmRLEDecoderRegistration::registerCodecs(OFFalse /*pCreateSOPInstanceUID*/, (opt_debugMode ? OFTrue : OFFalse));
	DJDecoderRegistration::registerCodecs(opt_decompCSconversion, EUC_default, EPC_default, (opt_debugMode ? OFTrue : OFFalse));

	DcmFileFormat *dfile = new DcmFileFormat();
	OFCondition cond = dfile->loadFile(dicomname, opt_transferSyntax, EGL_withoutGL, DCM_MaxReadLength, opt_readMode);

	E_TransferSyntax xfer = dfile->getDataset()->getOriginalXfer();

	DicomImage *di = new DicomImage(dfile, xfer, opt_compatibilityMode, opt_frame - 1, opt_frameCount);

	di->setMinMaxWindow(0);

	DiPixel* dmp = const_cast<DiPixel*>(di->getInterData());
	void *pixelData = NULL;
	pixelData = (void *)dmp->getData();
	Sint16* p = (Sint16*)pixelData;

	Sint16 tmp[ROWS * COLUMNS];
	int count = 0;
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLUMNS; j++)
		{
			Sint16 v = *p;
			tmp[count] = v;
			count++;
			p++;
		}
	}


	if (cond.good())
	{
		char uid[100];
		DcmFileFormat * fileformat = new DcmFileFormat();
		DcmDataset * copy = fileformat->getDataset();
		// Use DcmDataset's assignment operator to copy the complete dataset
		*copy = *dfile->getDataset();
		copy->putAndInsertUint16Array(DCM_PixelData, (const Uint16*)tmp, frameSize);

		OFCondition status = fileformat->saveFile(dicomname, EXS_LittleEndianImplicit);
		if (status.bad())
			std::cout << "Error: cannot write DICOM file (" << status.text() << ")" << std::endl;
	}
}

int CDicomReader::GetTransferSyntax(const wchar_t * inor)
{
	int transfersyntax = -1;

	size_t origsize = wcslen(inor) + 1;
	const size_t newsize = 500;
	size_t convertedChars = 0;

	char dicomname[newsize];

	wcstombs_s(&convertedChars, dicomname, origsize, inor, _TRUNCATE);

	DcmFileFormat *fileformat = new DcmFileFormat();
	OFCondition status = fileformat->loadFile(dicomname);

	if (status.good())
	{
		E_TransferSyntax xfer = fileformat->getDataset()->getOriginalXfer();
		transfersyntax = (int)xfer;
	}

	return transfersyntax;
}