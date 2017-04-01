#ifndef _HOUGH_
#define _HOUGH_

#include "canny.cpp"
#include <vector>
#include <algorithm> 
#include <iostream>

using namespace std;

#define PI 3.14159265359
#define threshold 100

struct Position {
    int rowindex;
    int colindex;
    int votenum;
};

struct Point {
    int x;
    int y;
    Point(int _x, int _y) {
        x = _x;
        y = _y;
    }
    Point() {
        x = -1;
        y = -1;
    }
};

class Hough {
public:
    CImg<unsigned char> edge;
    CImg<unsigned char> img;
    int edge_w;
    int edge_h;
    int acc_w;
    int acc_h;
    int center_h;
    vector<Point> vertex;
    int **acc;

    Hough(const char *filename);
    ~Hough();
    void detectline();
    vector<Position> vote();
    vector<Position> gethighestvote(vector<Position> v);
    void drawline(vector<Position> v);
    void drawpoint(vector<Position> v);
    
};

#endif