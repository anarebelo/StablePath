#include <iostream>
#include <string>
#include <vector> 
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <ctime>
#include <cmath>
#include <cfloat>

#include "external/opencv1.0/include/cv.h"		  // include core library interface
#include "external/opencv1.0/include/highgui.h"  // include GUI library interface
#include "FastAutoAdjust.h"
#include "assignmentoptimal.h"

//#define TWOTHREAD
//================================================***********************************
myIplImage * myCreateImage (int width, int height, unsigned char initValue)
{
	myIplImage * img = new myIplImage;
	img->height = height;
	img->width = width;
	img->widthStep = width;
	img->imageData = new unsigned char[img->widthStep*height];	
	img->imageDataOrigin = img->imageData;
	memset (img->imageData, initValue, img->widthStep*height);

	return img;
}
void myReleaseImage(myIplImage ** img)
{
	delete (*img)->imageDataOrigin;
	delete (*img);
	*img = 0;
}

void mySaveImage(const char * filename, myIplImage * img)
{
	IplImage * cvimg = cvCreateImage (cvSize (img->width, img->height), 8, 1);
	for(int row = 0; row < img->height; row++)
	{
        for (int col = 0; col < img->width; col++)
		{		
			unsigned char pel = img->imageData[row*img->widthStep + col];
			cvimg->imageData[row*cvimg->widthStep + col] = pel;
		}
	}
	cvSaveImage (filename, cvimg);
	cvReleaseImage (&cvimg);
}

void myVerticalErodeImage(myIplImage * img)
{
	for(int c = 0; c < img->width; c++)
	{
		unsigned char pel_prev = img->imageData[0*img->widthStep + c];
		unsigned char pel_curr = img->imageData[0*img->widthStep + c];
        for (int r = 0; r < img->height-1; r++)
		{		
			unsigned char pel_next = img->imageData[(r+1)*img->widthStep + c];
			unsigned char pel = get_min<unsigned char>(pel_prev, pel_curr, pel_next);
			img->imageData[r*img->widthStep + c] = pel;
			pel_prev = pel_curr;
			pel_curr = pel_next;
		}
		img->imageData[(img->height-1)*img->widthStep + c] = min(pel_prev, pel_curr);//last row
	}
}

myIplImage * myCloneImage (myIplImage * imgOri)
{
	myIplImage * img = new myIplImage;
	img->height = imgOri->height;
	img->width = imgOri->width;
	img->widthStep = imgOri->widthStep;
	img->imageData = new unsigned char[img->widthStep*img->height];
	memcpy (img->imageData, imgOri->imageData, img->widthStep*img->height);
	img->imageDataOrigin = img->imageData;
	return img;
}

//beginOf OPENCV dependence
myIplImage * myCloneImage (IplImage * imgOri)
{
	myIplImage * img = new myIplImage;
	img->height = imgOri->height;
	img->width = imgOri->width;
	img->widthStep = imgOri->widthStep;
	img->imageData = new unsigned char[img->widthStep*img->height];
	memcpy (img->imageData, imgOri->imageData, img->widthStep*img->height);
	img->imageDataOrigin = img->imageData;
	return img;
}
//endOf OPENCV dependence*/

//================================================***********************************

using namespace std;

#define BW_THRESHOLD (128+32)

const double staffLineFinder::MIN_BLACK_PER = 0.25;
const staffLineFinder::weight_t staffLineFinder::TOP_VALUE = static_cast<staffLineFinder::weight_t>(INT_MAX/2); 


int sumOfValuesInVector (vector<Point2D>& vec, myIplImage * img)
{	// WORKS FOR GRAY_SCALE IMAGES
	size_t len = vec.size();
	int sumOfValues = 0;
	int startCol = 0;
	int endCol = img->width-1; 
	for (int i = startCol; i <= endCol; i++)
	{
		int col = vec[i].x;
		int row = vec[i].y;
		unsigned char pel = img->imageData[(row)*img->widthStep + col];
		sumOfValues += pel;		
	}
	return sumOfValues;
}
//=============================================================================
template <class T> //HELPER FUNCTION
T findMostRepresentedValueOnSortedVector(vector<T>& vec)
{
	T run = vec[0];
	int freq = 0;
	int maxFreq = 0;
	T maxRun = run;

	for (unsigned i = 0; i< vec.size(); i++)
	{
		if (vec[i] == run)
			freq ++;
		if (freq > maxFreq)
		{
			maxFreq = freq;
			maxRun = run;
		}
		if (vec[i] != run)
		{
			freq = 0;
			run = vec[i];
		}
	}
	return maxRun;
}
//=============================================================================
double staffLineFinder::staffDissimilarity(vector<Point2D>& staff1, vector<Point2D>& staff2)
{
	double currAvg1 = 0;
	double currAvg2 = 0;
	for(size_t i = 0; i < staff1.size(); i++)
	{
		currAvg1 += staff1[i].y;
		currAvg2 += staff2[i].y;
	}
	currAvg1/=staff1.size();
	currAvg2/=staff2.size();
	double avgDiff = 0;
	for(size_t i = 0; i < staff1.size(); i++)
	{
		double curr_dif = abs(staff1[i].y-staff2[i].y - currAvg1+currAvg2);
		avgDiff += (curr_dif*curr_dif);
	}
	avgDiff /=staff1.size();
	avgDiff = sqrt (avgDiff);
	return avgDiff;
}

void staffLineFinder::computeMedianStaff (vector< vector <Point2D> > &staves, vector <Point2D> & medianStaff)
{	
	medianStaff.clear();
	for (int c = staves[0][0].x; c <= staves[0][staves[0].size()-1].x; c++)
	{
		vector <int> delta;
		for (size_t nvalid = 0; nvalid < staves.size(); nvalid++)
		{
			if (c != staves[0][0].x)
				delta.push_back(staves[nvalid][c].y-staves[nvalid][0].y);
			else
				delta.push_back(staves[nvalid][0].y);
		}
		sort(delta.begin(), delta.end());
		int y;
		if (c != staves[0][0].x) 
			y = medianStaff[0].y + delta[staves.size()/2];
		else
			y = delta[staves.size()/2];
		medianStaff.push_back(Point2D (c, y));
	}
}


void staffLineFinder::newFindStaffHeightAndDistance (myIplImage * img, vector<vector<Point2D> >&stablePaths, int& staffHeight, int& staffDistance)
{
	unsigned char WHITE = 255;
	vector<int> runs[2];

	for (int i=0; i <stablePaths.size(); i++)
	{
		vector<Point2D> &staff = stablePaths[i];
		for (int j=0; j <staff.size(); j++)
		{
			Point2D curr = staff[j];
			int col = curr.x;
			int row = curr.y;
			unsigned char pel = img->imageData[row*img->widthStep + col];
			if (pel == WHITE)
				continue;
			int runBlack = -1;
			while (pel == 0)
			{
				--row;
				++runBlack;
				if (row < 0)
					break;
				pel = img->imageData[row*img->widthStep + col];	
			}
			int runWHITE = 0;
			while (row >= 0 && (pel = img->imageData[row*img->widthStep + col]) == WHITE)
			{
				++runWHITE;
				--row;	
			}
			runs[1].push_back(runWHITE);
			row = curr.y;
			pel = img->imageData[row*img->widthStep + col];
			while (pel == 0)
			{
				++row;
				++runBlack;
				if (row>img->height-1)
					break;
				pel = img->imageData[row*img->widthStep + col];	
			}
			runWHITE = 0;
			while (row < img->height && (pel = img->imageData[row*img->widthStep + col]) == WHITE)
			{
				++runWHITE;
				++row;	
			}
			runs[1].push_back(runWHITE);

			runs[0].push_back(runBlack);

		}
	}
	// find most repeated
	sort(runs[0].begin(), runs[0].end());
	sort(runs[1].begin(), runs[1].end());

	staffDistance = -1;
	staffHeight = findMostRepresentedValueOnSortedVector<int>(runs[0]);
	staffDistance = findMostRepresentedValueOnSortedVector<int>(runs[1]);

	cout << "\n\t NEW staffHeight " << staffHeight << "  staff Distance " << staffDistance <<endl<< endl; 
	return;
}

//====================================================
//====================================================
inline staffLineFinder::weight_t staffLineFinder::weightFunction (int row1, int col1, int row2, int col2, staffLineFinder::NEIGHBOUR neigh)
{
	unsigned char pel1 = img->imageData[row1*img->widthStep + col1];
	unsigned char pel2 = img->imageData[row2*img->widthStep + col2];

	int d1 = verDistance[row1*img->width + col1];
	int d2 = verDistance[row2*img->width + col2];
	int vRun1 = verRun[row1*img->width + col1];
	int vRun2 = verRun[row2*img->width + col2];	

	int pel = min (pel1, pel2);
	int y0 = 4;
	int y255 = 8;
	if (neigh == staffLineFinder::NEIGHBOUR8) 
	{
		y0 = 6;
		y255 = 12;
	}
	int y = (pel == 0 ? y0:y255);
	if ( (pel == 0) && (min (vRun1, vRun2) <= staffLineHeight)) 
		--y;
	if (max(d1,d2) > 2*staffLineHeight+staffSpaceDistance) 
		++y; 
	return y;
}

bool staffLineFinder::tooMuchWhite (vector<Point2D>& vec, myIplImage * img, double minBlackPerc)
{	// WORKS FOR GRAY_SCALE IMAGES
	int sumOfValues = sumOfValuesInVector (vec, img);
	size_t len = vec.size();
	int startCol = 0;//0.1*len;
	int endCol = img->width-1; //0.9*len;
	int usedSize = endCol-startCol+1;
	if (sumOfValues > 255*(1-minBlackPerc)*(usedSize) )
	{
		//printf ("%f\n", 1- sumOfValues/(255*usedSize) );
		return true;
	}
	return false;
}

//=============================================================================
void trimPath (vector<unsigned char>& vec, int window, int&startX, int&endX)
{
	startX = 0;
	while(vec[startX] == 255 && startX < vec.size()/2)
	{
		startX++;
	}
	int i;	
	int sum = 0;
	for(i = startX; i< startX + window; i++)
	{
		sum += vec[i];
	}
	for( ; i< vec.size()/2; i++)
	{
		sum += vec[i];
		sum -= vec[i-window];
		if (sum > window*255*0.9) 
		{
			startX = i+1;
		}
	}

	endX = vec.size()-1;
	while(vec[endX] == 255 && endX > vec.size()/2)
	{
		endX--;
	}
	sum = 0;
	for(i = endX; i> endX - window; i--)
	{
		sum += vec[i];
	}
	for( ; i> vec.size()/2; i--)
	{
		sum += vec[i];
		sum -= vec[i+window];
		if (sum > window*255*0.9) 
		{
			endX = i-1;
		}
	}
}

//=============================================================================
void staffLineFinder::findStaffHeightAndDistance(myIplImage* img, int& staffHeight, int& staffDistance)
{
	unsigned char WHITE = 255;
	vector<int> runs[2];
	for (int c = 0; c < img->width; c++)
	{
		int run = 0;
		unsigned char val = WHITE;
		for(int r = 0; r < img->height; r++)
		{

			unsigned char pel = img->imageData[r*img->widthStep + c];
			//beginOf implicit horizontal dilate + erode to remove black noise
			unsigned char pel_left2 = WHITE;
			unsigned char pel_left1 = WHITE;
			unsigned char pel_right1 = WHITE;
			unsigned char pel_right2 = WHITE;
			if (c>0) pel_left1 = img->imageData[r*img->widthStep + c-1];
			if (c>1) pel_left2 = img->imageData[r*img->widthStep + c-2];
			if (c<img->width-1) pel_right1 = img->imageData[r*img->widthStep + c+1];
			if (c<img->width-2) pel_right2 = img->imageData[r*img->widthStep + c+2];
			//pel = get_min<unsigned char>(get_max<unsigned char> (pel_left2, pel_left1, pel),
			//	get_max<unsigned char> (pel_left1, pel, pel_right1),
			//	get_max<unsigned char> (pel, pel_right1, pel_right2) );
			pel = min(max(pel_left2, pel_left1), max (pel_left1, pel));
			//endOf implicit horizontal dilate + erode to remove black noise

			if (pel == val)
				run++;
			else
			{
				runs[val/WHITE].push_back(run);
				val = WHITE - val;
				run = 1;
			}
		}
		runs[val/WHITE].push_back(run);
	}
	// find most repeated
	sort(runs[0].begin(), runs[0].end());
	sort(runs[1].begin(), runs[1].end());

	staffHeight = findMostRepresentedValueOnSortedVector<int>(runs[0]);
	staffDistance = findMostRepresentedValueOnSortedVector<int>(runs[1]);

	cout << "staffHeight " << staffHeight << endl <<"staffDistance " << staffDistance << endl; 
	return;
}

staffLineFinder::staffLineFinder (const char * imgpath) : img_path (imgpath)
{
	globalStart = time (0);

	//beginOf OPENCV dependence
	IplImage* imgOri  = cvLoadImage (imgpath, 0);	// output image is always gray
	cvThreshold (imgOri, imgOri, BW_THRESHOLD, 255, CV_THRESH_BINARY); // make it binary
	img = myCloneImage(imgOri);
	cvReleaseImage(&imgOri);
	//endOf OPENCV dependence

	int img_width =img->width;
	int img_height=img->height;

	staffLineHeight = 0;
	staffSpaceDistance = img_height;

	graphPath = new NODE[img_width * img_height];
#ifdef TWOTHREAD
	t_graphPath = new NODE[img_width * img_height];
#else
	t_graphPath = 0;
#endif
	graphWeight = new NODEGRAPH[img_width * img_height];

	verRun = new int[img_height*img->width];
	verDistance = new int[img_height*img->width];
	memset (verDistance, 0, sizeof(int)*img_height*img_width);
}

staffLineFinder::~staffLineFinder () 
{
	myReleaseImage(&img);
	delete graphPath;
	delete t_graphPath;
	delete graphWeight;
	delete img;
	delete verRun;
	delete verDistance;
	printf ("\tGLOBAL TIME %d\n\n",	time (0)-globalStart);
}
//==================================================
int staffLineFinder::findShortestPath(myIplImage* img, int startCol, int endCol, std::vector<Point2D> &shortestPath)
{
	shortestPath.clear();
	int width = img->width;

	for (int row = 0; row < img->height; row++)
	{			
		graphPath[row*width + startCol].weight = static_cast<weight_t>(0.0);
		graphPath[row*width + startCol].start = Point2D(startCol, row);
	}
	for (int c = startCol+1; c <= endCol; c++)
	{
		int col = c;
		for (int row = 0; row < img->height; row ++)
		{	
			weight_t weight10, weight20, weight30;
			weight_t value1, value2, value3;
			weight10 = weight20 = weight30 = TOP_VALUE;
			value1 = value2 = value3  = TOP_VALUE; 

			if(row > 0)
			{
				weight10 = graphWeight[(row-1)*width + col-1].weight_down; 
				value1 = weight10 + graphPath[(row-1)*width + (col-1)].weight;
			}
			weight20 = graphWeight[row*width + col-1].weight_hor; 
			value2 = weight20 + graphPath[(row+0)*width + (col-1)].weight;
			if (row < img->height-1)
			{
				weight30 = graphWeight[(row+1)*width + col-1].weight_up; 
				value3 = weight30 + graphPath[(row+1)*width + (col-1)].weight;
			}
			//
			if ((value3)<= (value2) && (value3)<= (value1))
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row+1);
				graphPath[(row)*width + (col)].weight = value3;
				graphPath[(row)*width + (col)].start = graphPath[(row+1)*width + (col-1)].start;
			}
			else if ((value2)<= (value1) && (value2)<= (value3))
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row);
				graphPath[(row)*width + (col)].weight = value2;
				graphPath[(row)*width + (col)].start = graphPath[(row+0)*width + (col-1)].start;
			}
			else 
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row-1);
				graphPath[(row)*width + (col)].weight = value1;
				graphPath[(row)*width + (col)].start = graphPath[(row-1)*width + (col-1)].start;
			}
		}
	}
	int minIdx= 0;
	weight_t minCost = graphPath[0*width + endCol].weight;

	for (int j = 0; j < img->height; j++)
	{
		if (graphPath[j*width + endCol].weight  < minCost)
		{
			minIdx = j;
			minCost = graphPath[j*width + endCol].weight;
		}
	}

	Point2D p (endCol, minIdx);
	shortestPath.resize(endCol-startCol + 1);
	int pos = endCol-startCol;
	shortestPath[pos] = p;
	pos--;
	while(p.x != startCol)
	{
		p = graphPath[p.y*width + p.x].previous;
		shortestPath[pos] = p;
		pos--;
	}
	return minCost;
}

//==================================================

myIplImage* t_img;
int t_startCol;
int t_endCol;
vector<int> *pstartRow_i;
DWORD threadFunc (LPVOID lpParam) 
{
	staffLineFinder * obj = (staffLineFinder*)lpParam;
	obj->svc(t_img, t_startCol, t_endCol, *pstartRow_i);
	return 0;
}

int staffLineFinder::svc(myIplImage* img, int startCol, int endCol, vector<int>& startRow_i)
{
	int width = img->width;
	int endCol_i = img->width-1-startCol;
	int startCol_i = img->width-1-endCol;
	for (int row = 0; row < img->height; row++)
	{			
		t_graphPath[row*width + startCol_i].weight = static_cast<weight_t>(0);
		t_graphPath[row*width + startCol_i].start  = Point2D(startCol_i, row);
	}
	for (int col = startCol_i + 1; col <= endCol_i; col++)
	{
		for (int row = 0; row < img->height; row++)
		{	
			weight_t weight10, weight20, weight30;
			weight_t value1, value2, value3;
			weight10 = weight20 = weight30 = TOP_VALUE;
			value1 = value2 = value3  = TOP_VALUE; 

			if(row > 0)
			{
				weight10 = graphWeight[row*width + width-1-col].weight_up; 
				value1 = weight10 + t_graphPath[(row-1)*width+(col-1)].weight;
			}

			weight20 = graphWeight[row*width + width-1-col].weight_hor; 
			value2 = weight20 + t_graphPath[(row+0)*width + (col-1)].weight;
			
			if (row < img->height-1)
			{
				weight30 = graphWeight[row*width + width-1-col].weight_down; 
				value3 = weight30 + t_graphPath[(row+1)*width + (col-1)].weight;
			}
			//
			weight_t min_v = min(min(value3, value2), value1);			
			if ((value3)<= (value2) && (value3)<= (value1))
			{				
				t_graphPath[(row)*width + (col)].previous = Point2D(col-1, row+1);
				t_graphPath[(row)*width + (col)].weight = value3;
				t_graphPath[(row)*width + (col)].start = t_graphPath[(row+1)*width + (col-1)].start;
			}
			else if ((value2)<= (value1) && (value2)<= (value3))
			{
				t_graphPath[(row)*width + (col)].previous = Point2D(col-1, row);
				t_graphPath[(row)*width + (col)].weight = value2;
				t_graphPath[(row)*width + (col)].start = t_graphPath[(row+0)*width + (col-1)].start;
			}
			else
			{
				t_graphPath[(row)*width + (col)].previous = Point2D(col-1, row-1);
				t_graphPath[(row)*width + (col)].weight = value1;
				t_graphPath[(row)*width + (col)].start = t_graphPath[(row-1)*width + (col-1)].start;
			}	
		}
	}
	for (int i = 0; i < img->height; i++)
		startRow_i.push_back(t_graphPath[i*width + endCol_i].start.y);
	return 0;
}

int staffLineFinder::findAllStablePaths(myIplImage* img, int startCol, int endCol, std::vector <vector<Point2D> > &stablePaths)
{
	stablePaths.clear();
	int width = img->width;
	vector<int> startRow_i;
#ifdef TWOTHREAD
	DWORD ThreadId;
	pstartRow_i = &startRow_i;
	t_startCol = startCol;
	t_endCol = endCol;
	t_img = img;
	HANDLE threadHandle = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)threadFunc, this, 0, &ThreadId);
#else
	int endCol_i = img->width-1-startCol;
	int startCol_i = img->width-1-endCol;
	for (int row = 0; row < img->height; row++)
	{			
		graphPath[row*width + startCol_i].weight = static_cast<weight_t>(0);
		graphPath[row*width + startCol_i].start  = Point2D(startCol_i, row);
	}
	for (int col = startCol_i + 1; col <= endCol_i; col++)
	{
		for (int row = 0; row < img->height; row++)
		{	
			weight_t weight10, weight20, weight30;
			weight_t value1, value2, value3;
			weight10 = weight20 = weight30 = TOP_VALUE;
			value1 = value2 = value3  = TOP_VALUE; 

			if(row > 0)
			{
				weight10 = graphWeight[row*width + width-1-col].weight_up; 
				value1 = weight10 + graphPath[(row-1)*width+(col-1)].weight;
			}

			weight20 = graphWeight[row*width + width-1-col].weight_hor; 
			value2 = weight20 + graphPath[(row+0)*width + (col-1)].weight;
			
			if (row < img->height-1)
			{
				weight30 = graphWeight[row*width + width-1-col].weight_down; 
				value3 = weight30 + graphPath[(row+1)*width + (col-1)].weight;
			}
			//
			weight_t min_v = min(min(value3, value2), value1);			
			if ((value3)<= (value2) && (value3)<= (value1))
			{				
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row+1);
				graphPath[(row)*width + (col)].weight = value3;
				graphPath[(row)*width + (col)].start = graphPath[(row+1)*width + (col-1)].start;
			}
			else if ((value2)<= (value1) && (value2)<= (value3))
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row);
				graphPath[(row)*width + (col)].weight = value2;
				graphPath[(row)*width + (col)].start = graphPath[(row+0)*width + (col-1)].start;
			}
			else
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row-1);
				graphPath[(row)*width + (col)].weight = value1;
				graphPath[(row)*width + (col)].start = graphPath[(row-1)*width + (col-1)].start;
			}	
		}
	}
	for (int i = 0; i < img->height; i++)
		startRow_i.push_back(graphPath[i*width + endCol_i].start.y); //*/
#endif
	//-------------------------------------------------------------------------
	for (int row = 0; row < img->height; row++)
	{			
		graphPath[row*width + startCol].weight = static_cast<weight_t>(0.0);
		graphPath[row*width + startCol].start = Point2D(startCol, row);
	}
	for (int c = startCol+1; c <= endCol; c++)
	{
		int col = c;
		for (int row = 0; row < img->height; row ++)
		{	
			weight_t weight10, weight20, weight30;
			weight_t value1, value2, value3;
			weight10 = weight20 = weight30 = TOP_VALUE;
			value1 = value2 = value3  = TOP_VALUE; 

			if(row > 0)
			{
				weight10 = graphWeight[(row-1)*width + col-1].weight_down; 
				value1 = weight10 + graphPath[(row-1)*width + (col-1)].weight;
			}
			weight20 = graphWeight[row*width + col-1].weight_hor; 
			value2 = weight20 + graphPath[(row+0)*width + (col-1)].weight;
			if (row < img->height-1)
			{
				weight30 = graphWeight[(row+1)*width + col-1].weight_up; 
				value3 = weight30 + graphPath[(row+1)*width + (col-1)].weight;
			}
			//
			if ((value3)<= (value2) && (value3)<= (value1))
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row+1);
				graphPath[(row)*width + (col)].weight = value3;
				graphPath[(row)*width + (col)].start = graphPath[(row+1)*width + (col-1)].start;
			}
			else if ((value2)<= (value1) && (value2)<= (value3))
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row);
				graphPath[(row)*width + (col)].weight = value2;
				graphPath[(row)*width + (col)].start = graphPath[(row+0)*width + (col-1)].start;
			}
			else 
			{
				graphPath[(row)*width + (col)].previous = Point2D(col-1, row-1);
				graphPath[(row)*width + (col)].weight = value1;
				graphPath[(row)*width + (col)].start = graphPath[(row-1)*width + (col-1)].start;
			}
		}
	}
#ifdef TWOTHREAD
	WaitForSingleObject(threadHandle, INFINITE);
#endif
	for (int i = 0; i < img->height; i ++)
	{
		int startR = graphPath[i*width + endCol].start.y;

		if (startRow_i[startR] == i)
		{	// STABLE PATH
			Point2D p (endCol, i);
			vector<Point2D> contour;
			contour.resize(endCol-startCol + 1);
			int pos = endCol-startCol;
			contour[pos] = p;
			pos--;
			while(p.x != startCol)
			{
				p = graphPath[p.y*width + p.x].previous;
				contour[pos] = p;
				pos--;
			}
			stablePaths.push_back(contour);
		}
	}
/*{
	IplImage* imgDebug = cvLoadImage (img_path.c_str(), 1);// output image is always color
	cvZero (imgDebug);
	for (int nvalid = 0; nvalid < stablePaths.size(); nvalid++)
	{
		vector <Point2D> staff = stablePaths[nvalid];
		unsigned char color1 = 127+rand()%128;
		unsigned char color2 = 127+rand()%128;
		unsigned char color3 = 127+rand()%128;
		for (int i = 0; i<staff.size(); i++)
		{
			int col = staff[i].x;
			int row = staff[i].y;
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 2] = color1;
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 1] = color2;
			imgDebug->imageData[row*imgDebug->widthStep + 3*col + 0] = color3;
		}
	}
	string out_path = img_path;
	out_path = out_path.substr(0, out_path.rfind ('.'));
	out_path += "_staffLines_vs0x.bmp";
	cvSaveImage(out_path.c_str(), imgDebug);
	cvReleaseImage(&imgDebug);

exit(0);}*/
	return 0;
}

void staffLineFinder::constructGraphWeights ()
{
#ifndef ICIP2008
	unsigned char WHITE = 255;
	// vertical runs
	for (int c = 0; c < img->width; c++)
	{
		int run = 0;
		unsigned char val = WHITE;
		for(int r = 0; r < img->height; r++)
		{
			unsigned char pel = img->imageData[r*img->widthStep + c];
			if (pel == val)
				run++;
			else
			{
				int len = run;
				for(int row = r-1; len >0; len--, row--)
				{
					verRun[row*img->width + c] = run;
				}				
				val = WHITE - val;
				run = 1;
			}
		}
	}
	//
	for (int c = 0; c < img->width; c++)
	{	
		for(int r = 0; r < img->height; r++)
		{
			unsigned char pel = img->imageData[r*img->widthStep + c];
			if (pel == (unsigned char) 255) // white distance are not currently used: comment this and the next line is otherwise
				continue;
			int row = r;
			unsigned char curr_pel = pel;
			while (row>0 && curr_pel == pel)
			{
				row--;
				curr_pel = img->imageData[row*img->widthStep + c];
			}
			int run1 = 1;
			while (row>0 && curr_pel != pel)
			{
				row--;
				curr_pel = img->imageData[row*img->widthStep + c];
				run1++;
			}
			
			row = r;
			curr_pel = pel;
			while (row<img->height-1 && curr_pel == pel)
			{
				row++;
				curr_pel = img->imageData[row*img->widthStep + c];
			}
			int run2 = 1;
			while (row<img->height-1 && curr_pel != pel)
			{
				row++;
				curr_pel = img->imageData[row*img->widthStep + c];
				run2++;
			}
			verDistance [r*img->width + c] = min(run1, run2);
		}
	}
#endif
	for (int r = 0; r < img->height; r++)
	{
		for (int c = 0; c <img->width-1; c++)
		{
			graphWeight[r*img->width + c].weight_hor = weightFunction(r,c, r, c+1, NEIGHBOUR4);
			if (r>0)
				graphWeight[r*img->width + c].weight_up = weightFunction(r,c, r-1, c+1, NEIGHBOUR8);
			else
				graphWeight[r*img->width + c].weight_up = TOP_VALUE;
			if (r<img->height-1)
				graphWeight[r*img->width + c].weight_down = weightFunction(r,c, r+1, c+1, NEIGHBOUR8);
			else
				graphWeight[r*img->width + c].weight_down = TOP_VALUE;
		}
	}
}
//==================================================
void staffLineFinder::shortStaffDetection (vector <vector<Point2D> > &validStaves)
{
	printf ("\t\tstaffDetection_vs01\n");
	
	int staffHeight= 0, staffDistance = 0;
	findStaffHeightAndDistance (img, staffHeight, staffDistance);

	staffLineHeight = staffHeight;
	staffSpaceDistance = staffDistance;

	constructGraphWeights ();
	myIplImage* imgErode  = myCloneImage(img);
	myVerticalErodeImage (imgErode);
	myIplImage* imgErodedCopy  = myCloneImage (imgErode);

	time_t start = time (0);
	vector<int> npaths;

	bool first_time = 1;
	double blackPerc;
	vector<Point2D> bestStaff;

	while(1)
	{
		vector<Point2D> shortestPath;
		int curr_n_paths = 0;

		findShortestPath(img, 0, img->width-1, shortestPath);

		if (first_time && shortestPath.size() > 0)
		{ 
			first_time = 0;
		
			bestStaff.clear();
			size_t bestSumOfValues = INT_MAX;

			size_t sumOfValues = sumOfValuesInVector (shortestPath, imgErode);
			if (sumOfValues < bestSumOfValues)
			{
				bestSumOfValues = sumOfValues;
			}
			bestStaff = shortestPath;
		
			double bestBlackPerc = 1- bestSumOfValues/(255.0*(bestStaff.size())); // IMPORTANT THAT 255.0 be real
			blackPerc = max(MIN_BLACK_PER, bestBlackPerc*0.8); 
			printf ("bestBlackPerc blackPerc %f %f\n", bestBlackPerc, blackPerc);
		}

		vector<Point2D> staff = shortestPath;
		
		if (tooMuchWhite(staff, imgErode, blackPerc))
		{
			break;
		}
		double dissimilarity = staffDissimilarity(bestStaff, staff);
		if (dissimilarity > 4*staffDistance)
		{
			printf ("\tToo Dissimilar\n");
			break;
		}
	
		validStaves.push_back(staff);
		curr_n_paths++;

		int path_half_width2 = max(staffHeight, staffDistance/2);
		for (size_t i = 0; i<staff.size(); i++)
		{
			int col = staff[i].x;
			int row = staff[i].y;
			// ERASE PATHS ALREADY SELECTED!
			for (int j = -path_half_width2-1; j <= path_half_width2+1; j++)
			{	
				if ( ((row+j)>img->height-1) || ((row+j)<0) )
					continue;
				img->imageData[(row+j)*img->widthStep + col] = 255;
				imgErode->imageData[(row+j)*imgErode->widthStep + col] = 255;
			}
		}
		for (size_t i = 0; i<staff.size(); i++)
		{
			int col = staff[i].x;
			int row = staff[i].y;
			// ERASE PATHS ALREADY SELECTED!
			for (int j = -path_half_width2-1; j <= path_half_width2+1; j++)
			{	
				if ( ((row+j)>img->height-1) || ((row+j)<0) )
					continue;
				if (col == img->width-1)
					continue;
				if (row+j > 0)
					graphWeight[(row+j)*img->width + col].weight_up = 12; //weightFunction(row+j, col, row+j-1, col+1, NEIGHBOUR8);//12
				else
					graphWeight[(row+j)*img->width + col].weight_up = TOP_VALUE;
				graphWeight[(row+j)*img->width + col].weight_hor = 8; //weightFunction(row+j, col, row+j, col+1, NEIGHBOUR4); //8
				if (row+j < img->height-1)
					graphWeight[(row+j)*img->width + col].weight_down = 12; //weightFunction(row+j, col, row+j+1, col+1, NEIGHBOUR8);//12
				else
					graphWeight[(row+j)*img->width + col].weight_down = TOP_VALUE;
			}
		}
		npaths.push_back(curr_n_paths);
		if (curr_n_paths == 0)
			break;
	}
	printf ("--MAIN CYCLE DURATION %d\n", time(0)-start);

	postProcessing(validStaves, staffDistance, imgErodedCopy);

	printf ("\n\n\t TOTAL time %d\n", time (0)-start);
	printf ("Npaths: ");
	int totalpaths = 0;
	for (int i = 0; i < npaths.size(); i++)
	{
		totalpaths += npaths[i];
        printf ("Npaths = %d ", npaths[i]);
	}
	printf ("TOTAL = %d\n", totalpaths);
	myReleaseImage(&imgErode);
	myReleaseImage(&imgErodedCopy);
}
//==================================================
void staffLineFinder::stableStaffDetection (vector <vector<Point2D> > &validStaves)
{
	printf ("\t\tstaffDetection_vs01\n");
	
	findStaffHeightAndDistance (img, staffLineHeight, staffSpaceDistance);
	constructGraphWeights ();
	//staffLineHeight    = max(staffLineHeight-1, 1);	
	//staffSpaceDistance = max(staffSpaceDistance-1, 1);

	myIplImage* imgErode  = myCloneImage(img);
	myVerticalErodeImage (imgErode);
	myIplImage* imgErodedCopy  = myCloneImage (imgErode);

	time_t start = time (0);
	vector<int> npaths;

	bool first_time = 1;
	double blackPerc;
	vector<Point2D> bestStaff;

	while(1)
	{
		vector <vector<Point2D> > stablePaths;
		int curr_n_paths = 0;

		findAllStablePaths(img, 0, img->width-1, stablePaths);

		if (first_time && stablePaths.size() > 0)
		{ 
			first_time = 0;
			
			//newFindStaffHeightAndDistance (img, stablePaths, staffLineHeight, staffSpaceDistance);

			bestStaff.clear();
			size_t bestSumOfValues = INT_MAX;
			size_t bestStaffIdx = 0;
			vector<size_t> allSumOfValues;
			for (size_t c = 0; c < stablePaths.size(); c++)
			{
				size_t sumOfValues = sumOfValuesInVector (stablePaths[c], imgErode);
				if (1- sumOfValues/(255.0*(stablePaths[c].size())) > MIN_BLACK_PER) // minimal acceptable black percentage
					allSumOfValues.push_back(sumOfValues);
				if (sumOfValues < bestSumOfValues)
				{
					bestSumOfValues = sumOfValues;
					bestStaffIdx = c;
				}
			}
			///bestStaff = stablePaths[bestStaffIdx];
			vector<size_t> copy_allSumOfValues = allSumOfValues; 
			sort(allSumOfValues.begin(),allSumOfValues.end());
			size_t medianSumOfValues = allSumOfValues[allSumOfValues.size()/2];
			int i;
			for (i = 0; i < copy_allSumOfValues.size(); i++)
				if (copy_allSumOfValues[i] == medianSumOfValues)
					break;
			bestStaff = stablePaths[i];
		
			double bestBlackPerc = 1- medianSumOfValues/(255.0*(bestStaff.size())); // IMPORTANT THAT 255.0 be real
			blackPerc = max(MIN_BLACK_PER, bestBlackPerc*0.8); 
			printf ("bestBlackPerc blackPerc %f %f\n", bestBlackPerc, blackPerc);
		}
		
		for (size_t i = 0; i < stablePaths.size(); i++)
		{
			vector<Point2D> staff = stablePaths[i];
			
			if (tooMuchWhite(staff, imgErode, blackPerc))
			{
				continue;
			}
			double dissimilarity = staffDissimilarity(bestStaff, staff);
			if (dissimilarity > 4*staffSpaceDistance)
			{
				printf ("\tToo Dissimilar\n");
				continue;
			}
		
			validStaves.push_back(staff);
			curr_n_paths++;

			int path_half_width2 = max(staffLineHeight, staffSpaceDistance/2);
			for (size_t i = 0; i<staff.size(); i++)
			{
				int col = staff[i].x;
				int row = staff[i].y;

				// ERASE PATHS ALREADY SELECTED!
				for (int j = -path_half_width2-1; j <= path_half_width2+1; j++)
				{	
					if ( ((row+j)>img->height-1) || ((row+j)<0) )
						continue;
					img->imageData[(row+j)*img->widthStep + col] = 255;
					imgErode->imageData[(row+j)*imgErode->widthStep + col] = 255;
					if ( ((row+j)>img->height-1) || ((row+j)<0) )
						continue;
					if (col == img->width-1)
						continue;
					if (row+j > 0)
						graphWeight[(row+j)*img->width + col].weight_up = 12;
					else
						graphWeight[(row+j)*img->width + col].weight_up = TOP_VALUE;
					graphWeight[(row+j)*img->width + col].weight_hor = 8;
					if (row+j < img->height-1)
						graphWeight[(row+j)*img->width + col].weight_down = 12;
					else
						graphWeight[(row+j)*img->width + col].weight_down = TOP_VALUE;
				} //
			}
		}
		npaths.push_back(curr_n_paths);
		if (curr_n_paths == 0)
			break;
	}
	printf ("--MAIN CYCLE DURATION %d\n", time(0)-start);

	postProcessing(validStaves, staffSpaceDistance, imgErodedCopy);

	printf ("\n\n\t TOTAL time %d\n", time (0)-start);
	printf ("Npaths: ");
	int totalpaths = 0;
	for (int i = 0; i < npaths.size(); i++)
	{
		totalpaths += npaths[i];
        printf ("Npaths = %d ", npaths[i]);
	}
	printf ("TOTAL = %d TOTAL STAFF LINES %d\n", totalpaths, validStaves.size());
	myReleaseImage(&imgErode);
	myReleaseImage(&imgErodedCopy);
}

double simplifiedAvgDistance(std::vector<Point2D>&staff1, std::vector<Point2D>& staff2)
{
	if ((staff1.size() == 0) || (staff2.size() == 0))
	{
		printf ("SIZE 0\n");
		return -1;
	}
	int simplifiedDistance = 0;
	int m = max (staff1[0].x, staff2[0].x);
	int M = min (staff1[staff1.size()-1].x, staff2[staff2.size()-1].x);
	int idx1 = 0, idx2 = 0;
	if (m > M)
	{
		printf ("Do not overlap\n");
		return -1;
	}
	while (staff1[idx1].x != m)
		idx1++;
	while (staff2[idx2].x != m)
		idx2++;
	while (staff1[idx1].x != M)
	{
		int dy = abs(staff1[idx1].y-staff2[idx2].y);
		simplifiedDistance += abs(dy);
		idx1++;
		idx2++;
	}
	return static_cast<double>(simplifiedDistance)/(M-m+1);
}

void staffLineFinder::postProcessing (vector <vector<Point2D> > &validStaves, int staffDistance, myIplImage* imgErodedCopy)
{
	if (validStaves.size() == 0)
		return;
	// UNCROSS LINES
	for (int c = validStaves[0][0].x; c <= validStaves[0][validStaves[0].size()-1].x; c++)
	{
		vector <Point2D> column;
		for (int nvalid = 0; nvalid < validStaves.size(); nvalid++)
		{
			column.push_back(validStaves[nvalid][c]);
		}
		sort(column.begin(), column.end());
		for (int nvalid = 0; nvalid < validStaves.size(); nvalid++)
		{
			validStaves[nvalid][c] = column[nvalid];
		}
	}

	vector<Point2D> medianStaff;
	computeMedianStaff(validStaves, medianStaff);

	int maxStaffDistance = max(2*staffDistance, 2*staffLineHeight);

	// ORGANIZE IN SETS
	vector <vector <vector<Point2D> > > setsOfValidStaves;
	int start = 0;
	for (size_t nvalid = 0; nvalid < validStaves.size(); nvalid++)
	{
		if (nvalid== validStaves.size()-1 || simplifiedAvgDistance (validStaves[nvalid], validStaves[nvalid+1]) > maxStaffDistance)			
		{
			vector <vector<Point2D> > singleSet;
			for(int i = start; i <= nvalid; i++)
			{
				singleSet.push_back(validStaves[i]);
			}			
			if (singleSet.size() > 2)//TODO: more complex rules to validate sets
			{
				setsOfValidStaves.push_back(singleSet);
				printf ("SET SIZE =%d\n", singleSet.size() );	
			}
			start = nvalid+1;
		}
	}

	// UNDOCUMENTED OPERATION
	for (int i = 0; i < setsOfValidStaves.size(); i++)
	{
		vector <vector<Point2D> > &setOfStaves = setsOfValidStaves[i];
		for (int nvalid = 0; nvalid < setOfStaves.size(); nvalid++)
		{
			for (int deltacolumn = 2, sgn = 1; deltacolumn <imgErodedCopy->width; deltacolumn++, sgn = (-1)*sgn)
			{
				int c = imgErodedCopy->width/2 + (deltacolumn>>1)*sgn; //starts on the middle of the img on both sides
				int y = setOfStaves[nvalid][c].y;
				int x = setOfStaves[nvalid][c].x;
				int my = medianStaff[c].y;
				int y0 = setOfStaves[nvalid][imgErodedCopy->width/2].y;
				int my0 = medianStaff[imgErodedCopy->width/2].y;
				double alpha = 0; 
				unsigned char pel = imgErodedCopy->imageData[y*imgErodedCopy->widthStep + x];
				if (pel != 0)
					alpha = pow(abs(c-imgErodedCopy->width/2)/double(imgErodedCopy->width/2), 1/4.0);

				int delta = static_cast<int>( (1-alpha)*(y-y0) + alpha*(my-my0) );

				y = y0+delta;
				int prev_y = setOfStaves[nvalid][c-sgn].y;
				if ((y-prev_y)>1) y = prev_y+1;
				if ((y-prev_y)<-1) y = prev_y-1;

				setOfStaves[nvalid][c].y = y;
			}
		}
	}

	// TRIM AND SMOOTH
	//for (int i = 0; i < setsOfValidStaves.size(); i++)
	vector <vector <vector<Point2D> > >::iterator set_it = setsOfValidStaves.begin();
	while (set_it != setsOfValidStaves.end())
	{
		vector <vector<Point2D> > &setOfStaves = *set_it; //setsOfValidStaves[i];
		// compute median staff in terms of colour
		vector<unsigned char> medianStaff;			
		for(int c = 0; c < imgErodedCopy->width; c++)
		{
			vector <unsigned char> medianValue;
			for(int i = 0; i <setOfStaves.size(); i++)
			{
				int x = setOfStaves[i][c].x;
				int y = setOfStaves[i][c].y;
				unsigned char pel = imgErodedCopy->imageData[y*imgErodedCopy->widthStep + x];
				medianValue.push_back(pel);
			}
			sort(medianValue.begin(), medianValue.end());
			medianStaff.push_back(medianValue[medianValue.size()/5]);
		}
		//1 find start and end
		int startx = 0, endx = imgErodedCopy->width-1;
		::trimPath (medianStaff, 2*staffDistance, startx, endx);
		if ( (endx-startx) < maxStaffDistance) // 
		{	// remove the whole set
			set_it = setsOfValidStaves.erase(set_it);
			continue;
		}
		//2 trim staffs from start to nvalid
		for(int i = 0; i < setOfStaves.size(); i++)
		{
			smoothStaffLine (setOfStaves[i], 2*staffDistance);
			continue; // DO TO ERASE PATH
			vector<Point2D>::iterator it = setOfStaves[i].begin();
			while (it->x != startx)
			{
				it++;
			}
			setOfStaves[i].erase(setOfStaves[i].begin(), it);

			it = setOfStaves[i].begin();
			while (it->x != endx)
				it++;
			setOfStaves[i].erase(it, setOfStaves[i].end() );
		}
		set_it ++;
	}
	validStaves.clear(); // the number of staff lines may be different
	int j = 0;
	for (int i = 0; i < setsOfValidStaves.size(); i++)
	{
		vector <vector<Point2D> > &setOfStaves = setsOfValidStaves[i];
		for (int s = 0; s<setOfStaves.size(); s++)
		{
			validStaves.push_back(setOfStaves[s]);
		}
	}
}

void staffLineFinder::smoothStaffLine (vector<Point2D> &staff, int halfwindowsize)
{	//TO IMPROVE: DO IT INPLACE
	if (staff.size() < halfwindowsize*2+1)
	{
		return;
	}
	// low pass filter
	vector <Point2D> cpStaff = staff;
	int accumY=0;
	for (int i = 0 ; i < halfwindowsize*2; i++)
	{
		accumY += cpStaff[i].y;
		staff[i/2].y = accumY/(i+1);
	}
	for (int i = halfwindowsize; i < staff.size()-halfwindowsize; i++)
	{
		accumY += cpStaff[i + halfwindowsize].y;
		staff[i].y = accumY/(halfwindowsize*2 + 1);
		accumY -= cpStaff[i - halfwindowsize].y;
	}
	accumY=0;
	for (int i = 0 ; i < halfwindowsize*2; i++)
	{
		accumY += cpStaff[staff.size()-1-i].y;
		staff[staff.size()-1-i/2].y = accumY/(i+1);
	}
}