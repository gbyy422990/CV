#include "hough.cpp"

int main() {

    Hough h1("./Dataset/img1.bmp");
    h1.detectline();
    Hough h2("./Dataset/img2.bmp");
    h2.detectline();
    Hough h3("./Dataset/img3.bmp");
    h3.detectline();
    Hough h4("./Dataset/img4.bmp");
    h4.detectline();
    Hough h5("./Dataset/img5.bmp");
    h5.detectline();

    return 0;
}