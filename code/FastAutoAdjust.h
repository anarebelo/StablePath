#ifndef STAFFLINEFINDER_H
#define STAFFLINEFINDER_H

#include <vector>

//#define ICIP2008
//#define TRIMDIST

template <class T> 
T findMostRepresentedValueOnSortedVector(std::vector<T>& vec);

typedef struct Point2D
{
	int x;
	int y;
	Point2D () : x(0),y(0){}
	Point2D (int a, int b) : x(a),y(b){}
	Point2D (const Point2D& p):x(p.x),y(p.y){}
	bool operator<(Point2D &p){return y < p.y;}
}
Point2D;

typedef struct _myIplImage
{
	int  width;          
	int  height;        
	unsigned char *imageData;    
	int  widthStep;     
	unsigned char *imageDataOrigin; 
}
myIplImage;

myIplImage * myCreateImage (int width, int height, unsigned char initValue);
void myReleaseImage(myIplImage ** img);
void myErodeImage(myIplImage * img);
myIplImage * myCloneImage (myIplImage * imgOri);
myIplImage * myCloneImage (IplImage * imgOri);

template <class T>
static T get_min (T m1, T m2, T m3)
{
	return min(min(m1, m2), m3);
}
template <class T>
static T get_max (T m1, T m2, T m3)
{
	return max(max(m1, m2), m3);
}
class staffLineFinder
{	
public:
	typedef int weight_t;
	enum e_NEIGHBOUR {NEIGHBOUR4 = 0, NEIGHBOUR8};
	typedef enum e_NEIGHBOUR NEIGHBOUR;
private:

	struct NODE {
		Point2D previous;
		weight_t weight;
		Point2D start;
	};
	struct NODEGRAPH {
		weight_t weight_up;
		weight_t weight_hor;
		weight_t weight_down;
	};

	static const double MIN_BLACK_PER;
	static const weight_t TOP_VALUE;

	myIplImage* img;
	int* verRun;	  //vertical run of pixels of the same colour
	int* verDistance; //minimum distance of a pixel of the same colour NOT in the same run of pixels of the same colour
	NODE* graphPath;
	NODE* t_graphPath;
	NODEGRAPH * graphWeight;

	std::string img_path;
	int staffLineHeight;
	int staffSpaceDistance;
	time_t globalStart;

	void constructGraphWeights ();
	inline weight_t weightFunction (int row1, int col1, int row2, int col2, staffLineFinder::NEIGHBOUR neigh);
	int findAllStablePaths(myIplImage* img, int startCol, int endCol, std::vector <std::vector<Point2D> > &strongPaths);
	int findShortestPath(myIplImage* img, int startCol, int endCol, std::vector<Point2D> &shortestPath);
	static bool tooMuchWhite (std::vector<Point2D>& vec, myIplImage * img, double minBlackPerc);
	void postProcessing (std::vector <std::vector<Point2D> > &validStaves, int staffDistance, myIplImage* imgErodedCopy);
	void smoothStaffLine (std::vector<Point2D> &validStaves, int halfwindowsize);
	double staffDissimilarity(std::vector<Point2D>& staff1, std::vector<Point2D>& staff2);
	void computeMedianStaff (std::vector< std::vector <Point2D> > &staves, std::vector <Point2D> & medianStaff);

public:
	static void findStaffHeightAndDistance(myIplImage* img, int& staffHeight, int& staffDistance);
	static void newFindStaffHeightAndDistance (myIplImage * img, std::vector<std::vector<Point2D> >&stablePaths, int& staffHeight, int& staffDistance);

	staffLineFinder (const char * imgpath);
	~staffLineFinder ();
	void stableStaffDetection (std::vector <std::vector<Point2D> > &validStaves);
	void shortStaffDetection (std::vector <std::vector<Point2D> > &validStaves);
	int svc (myIplImage* img, int startCol, int endCol, std::vector<int>& startRow_i);
};


#endif // STAFFLINEFINDER_H
