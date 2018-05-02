#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <list>

#include <climits>
#include <ctime>
#include <cstdio>
#include <ctime>
#include <cfloat>

#include "stableStaffLineFinder.h"


using namespace std;
void staffLineRemoval1 (std::string img_path, std::vector <std::vector<Point2D> > &validStaves)
{
	IplImage* imgRem = cvLoadImage (img_path.c_str(), 0);// output image is always gray
	//cvThreshold (imgRem, imgRem, 128, 255, CV_THRESH_BINARY); // make it binary
	//cvThreshold (imgRem, imgRem, 0, 255, CV_THRESH_BINARY|CV_THRESH_OTSU); // make it binary

	int staffDistance, staffHeight;
	myIplImage *aux = new myIplImage;
	aux->height = imgRem->height;
	aux->imageData = (unsigned char*) imgRem->imageData;
	aux->imageDataOrigin = (unsigned char*) imgRem->imageDataOrigin;
	aux->width = imgRem->width;
	aux->widthStep = imgRem->widthStep;


	{
		stableStaffLineFinder slf (img_path.c_str()); 
		slf.newFindStaffHeightAndDistance (aux, validStaves, staffHeight, staffDistance);
	}

	delete aux;

	int threshold = 2*staffHeight;
	int tolerance = 1+ceil(staffHeight/3.0);
	printf ("tolerance %d\n", tolerance);
	for (int nvalid = 0; nvalid < validStaves.size(); nvalid++)
	{
		vector <Point2D> staff = validStaves[nvalid];
		for (int i = 0; i<staff.size(); i++)
		{
			int col = staff[i].x;
			int refRow = staff[i].y;
			int row = refRow;

			int run = 0;
			unsigned char pel = imgRem->imageData[row*imgRem->widthStep + col];
			if (pel == 255)
			{
				int dist1 = 0;
				while (pel == 255)
				{
					--row;
					++dist1;
					if (row<0)
						break;
					pel = imgRem->imageData[row*imgRem->widthStep + col];	
				}
				if (row<0)
					dist1 = imgRem->height-1;
				int dist2 = 0;
				row = refRow;
				pel = imgRem->imageData[row*imgRem->widthStep + col];
				while (pel == 255)
				{
					++row;
					++dist2;
					if (row>imgRem->height-1)
						break;
					pel = imgRem->imageData[row*imgRem->widthStep + col];	
				}
				if (row>imgRem->height-1)
					dist2 = imgRem->height-1;
				if (dist1 <= max(1, min (dist2, tolerance)))
				{
					refRow-=dist1;
				}
				else if (dist2 <= max(1, min (dist1, tolerance)))
				{
					refRow+=dist2;
				}
				else 
					continue;
			}
			row = refRow;
			pel = imgRem->imageData[row*imgRem->widthStep + col];	
			while (pel == 0)
			{
				++run;
				--row;
				if (row<0)
					break;
				pel = imgRem->imageData[row*imgRem->widthStep + col];	
			}
			row = refRow;
			pel = imgRem->imageData[row*imgRem->widthStep + col];	
			while (pel == 0)
			{
				++row;
				if (row>imgRem->height-1)
					break;
				pel = imgRem->imageData[row*imgRem->widthStep + col];	
				if(pel == 0)
					++run;
			}	
			if (run >= threshold)
				continue;
			// REMOVE
			row = refRow;
			pel = imgRem->imageData[row*imgRem->widthStep + col];	
			while (pel == 0)
			{
				imgRem->imageData[row*imgRem->widthStep + col] = 255;
				--row;
				if (row<0)
					break;
				pel = imgRem->imageData[row*imgRem->widthStep + col];	
			}
			row = refRow + 1;
			if (row > imgRem->height-1)
				continue;
			pel = imgRem->imageData[row*imgRem->widthStep + col];	
			while (pel == 0)
			{
				imgRem->imageData[row*imgRem->widthStep + col] = 255;
				++row;
				if (row> imgRem->height-1)
					break;
				pel = imgRem->imageData[row*imgRem->widthStep + col];	
			}
		}
	}

	string out_rem = img_path;
	out_rem = out_rem.substr(0, out_rem.rfind ('.'));
	out_rem += "-stable_removal1.bmp";
	cout << out_rem.rfind ('.') << endl;
	cvSaveImage(out_rem.c_str(), imgRem);
	cvReleaseImage(&imgRem);
}

using namespace std;
void saveDebugInfo(std::string img_path, std::vector <std::vector<Point2D> > &validStaves)
{
	IplImage* imgDebug = cvLoadImage (img_path.c_str(), 1);// output image is always color
	string out_textpath = img_path;
	out_textpath = out_textpath.substr(0, out_textpath.rfind ('.'));
	out_textpath += "-staffLines_vs01.txt";

	ofstream textfile (out_textpath.c_str());
	for (int nvalid = 0; nvalid < validStaves.size(); nvalid++)
	{
		vector <Point2D> staff = validStaves[nvalid];
		std::ostringstream outss;
		unsigned char color1 = 127+rand()%128;
		unsigned char color2 = 127+rand()%128;
		unsigned char color3 = 127+rand()%128;
		for (int i = 0; i<staff.size(); i++)
		{
			int col = staff[i].x;
			int row = staff[i].y;
			outss << col+1 << "," << row +1<< ",";
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 2] = color1;
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 1] = color2;
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 0] = color3;
		}
		textfile << outss.str() << endl;
	}
	string out_path = img_path;
	out_path = out_path.substr(0, out_path.rfind ('.'));
	out_path += "_staffLines_vs01.bmp";

	cvSaveImage(out_path.c_str(), imgDebug);

	printf("Output written to %s\n", out_path.c_str());
	cvReleaseImage(&imgDebug);
}

void LoadFileNames (string base, vector <string> &collection, const char * const filter)
{
	if (base.end () [-1] != '/' && base.end () [-1] != '\\')
		base += '/';

	WIN32_FIND_DATA fd;
	HANDLE h = FindFirstFile ((base + filter).c_str (), &fd);
	if (h != INVALID_HANDLE_VALUE) 
	{
		bool done = false;
		while (!done) 
		{
			if ((fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0)
				collection.push_back (base + fd.cFileName);
			done = !FindNextFile (h, &fd);
		}
		FindClose (h);
	}

	return ; // NOT IN SUBFOLDERS - uncomment to reverse 
	h = FindFirstFile ((base + "*").c_str (), &fd);

	if (h != INVALID_HANDLE_VALUE) 
	{
		bool done = false;
		while (!done) 
		{
			const string name = fd.cFileName;
			if (name != "." && name != "..") 
			{
				if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
					LoadFileNames (base + fd.cFileName, collection, filter);
			}
			done = !FindNextFile (h, &fd);
		}
		FindClose (h);
	}
}


int main(int argc, char *argv[])
{
#define ndeform 6
	vector<string> file_collection [ndeform];
	vector<double> time_vec;
	vector<string>::iterator theIterator;
	string curr_file;
	string basedir;

	basedir = argv[1];

	LoadFileNames (basedir.c_str (), file_collection[0], "*.png"); 


	for (int ndef = 0;ndef < ndeform; ndef++)
		std::sort (file_collection[ndef].begin(), file_collection[ndef].end());


	for (int ndef = 0;ndef < ndeform; ndef++)
	{
		time_t start = time (0);
		for (size_t i = 0; i < file_collection[ndef].size(); i++)
		{
			curr_file = file_collection[ndef][i];
			printf("\n\tProcessing FILE %s\n\n", curr_file.c_str());
			float threshold = -1;
			if (curr_file.find("bach")!=string::npos)
				threshold = 3;
			if (curr_file.find("baroque")!=string::npos)
				threshold = 4;
			if (curr_file.find("bellinzani")!=string::npos)
				threshold = 3;
			if (curr_file.find("brahms02")!=string::npos)
				threshold = 1.1;
			if (curr_file.find("bruckner01")!=string::npos)
				threshold = 1.1;
			if (curr_file.find("buxtehude")!=string::npos)
				threshold = 3;
			if (curr_file.find("carcassi")!=string::npos)
				threshold = 3;
			if (curr_file.find("dacrema")!=string::npos)
				threshold = 3;
			if (curr_file.find("dalitz03")!=string::npos)
				threshold = 3;
			if (curr_file.find("demoy")!=string::npos)
				threshold = 3;
			if (curr_file.find("derore01")!=string::npos)
				threshold = 3;
			if (curr_file.find("diabelli")!=string::npos)
				threshold = 3;
			if (curr_file.find("dowland")!=string::npos)
				threshold = 2;
			if (curr_file.find("gaultier01")!=string::npos)
				threshold = 4;
			if (curr_file.find("goess")!=string::npos)
				threshold = 2;
			if (curr_file.find("gregorian")!=string::npos)
				threshold = 2;
			if (curr_file.find("hotman")!=string::npos)
				threshold = 2;		
			if (curr_file.find("hurel")!=string::npos)
				threshold = 4;		
			if (curr_file.find("josquin04")!=string::npos)
				threshold = 3;		
			if (curr_file.find("mahler")!=string::npos)
				threshold = 2;
			if (curr_file.find("ockeghem")!=string::npos)
				threshold = 2;
			if (curr_file.find("pmw01")!=string::npos)
				threshold = 1.1;
			if (curr_file.find("pmw03")!=string::npos)
				threshold = 2;
			if (curr_file.find("pmw04")!=string::npos)
				threshold = 2;
			if (curr_file.find("rameau")!=string::npos)
				threshold = 3;
			if (curr_file.find("rossi")!=string::npos)
				threshold = 3;
			if (curr_file.find("schumann")!=string::npos)
				threshold = 3;
			if (curr_file.find("tye")!=string::npos)
				threshold = 3;
			if (curr_file.find("victoria09")!=string::npos)
				threshold = 3;
			if (curr_file.find("wagner")!=string::npos)
				threshold = 2;
			if (curr_file.find("williams")!=string::npos)
				threshold = 3;

			if (threshold < 0)
			{
				threshold = 3;
				printf ("UNKNOWN SCORE\n");
			}		
			

			vector <vector<Point2D> > validStaves;
			{
				stableStaffLineFinder slf1 (curr_file.c_str());
				slf1.stableStaffDetection(validStaves,curr_file.c_str());
				saveDebugInfo(curr_file, validStaves);
			}



			staffLineRemoval1 (curr_file, validStaves);
		}
		time_vec.push_back((time (0)-start));
		printf ("TOTAL TIME Deformation %d TIME %d\n", ndef, (time (0)-start));
	}
	for (int ndef = 0; ndef < time_vec.size(); ndef++)
	{
		printf ("\t Deformation %d TIME %f Number scores %d\n", ndef, time_vec[ndef], file_collection[ndef].size());
	}
	return 0;
}
