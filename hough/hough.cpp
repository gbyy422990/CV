#include "hough.h"

bool cmp(Position &a, Position &b) {
    return a.votenum > b.votenum;
}

Hough::Hough(const char *filename) {
    Canny can(filename);
    can.performcanny();        //找出图像边缘
    edge = can.edge;        //边缘图像
    img = can.img;        //原图像
    
    edge_w = edge.width();
    edge_h = edge.height();

    acc_w = 180;
    acc_h = 2 * ceil(sqrt(2) * (edge_w > edge_h ? edge_w : edge_h));        //r最大取原图像的对角线的长度
    center_h = acc_h / 2;

    acc = new int*[acc_h];
    for (int i = 0; i < acc_h; i++) {
        acc[i] = new int [acc_w];
        for (int j = 0; j < acc_w; j++)
            acc[i][j] = 0;
    }
}

Hough::~Hough() {
    for (int i = 0; i < acc_h; i++)
        delete [] acc[i];
    delete [] acc;
}

void Hough::detectline() {
    edge.display();
    vector<Position> v = vote();
    v = gethighestvote(v);
    drawline(v);
    drawpoint(v);
    img.display();
}

vector<Position> Hough::vote() {
    cimg_forXY(edge, x, y) {        //对每个边缘点在另一个坐标系上进行投票
        if (edge(x, y) > 0) {        //该点是边缘点。
            for (int i = 0; i < 180; i++) {
                double angle = (double)i / 180 * PI;        //角度制转换为弧度值
                double dr = (double)x * cos(angle) + (double)y * sin(angle);
                int r = round(dr);
                acc[r + center_h][i]++;        //r可为负值，加上矩阵中心
            }
        }
    }

    vector<Position> v;

    for (int i = 0; i < acc_h; i++) {        //找出投票数局部最大的r， theta组合。
        for (int j = 0; j < acc_w; j++) {
            if (acc[i][j] > threshold) {
                int flag = 1;
                if (i > 0) {
                    if (j > 0) {
                        if (acc[i][j] < acc[i-1][j-1]) flag = 0;
                        if (acc[i][j] < acc[i][j-1]) flag = 0;
                    }
                    if (j < acc_w - 1) {
                        if (acc[i][j] < acc[i-1][j+1]) flag = 0;
                        if (acc[i][j] < acc[i][j+1]) flag = 0;
                    }
                    if (acc[i][j] < acc[i-1][j]) flag = 0;
                }
                if (i < acc_h - 1) {
                    if (j > 0) {
                        if (acc[i][j] < acc[i+1][j-1]) flag = 0;
                    }
                    if (j < acc_w - 1) {
                        if (acc[i][j] < acc[i+1][j+1]) flag = 0;
                    }
                    if (acc[i][j] < acc[i+1][j]) flag = 0;
                }

                if (flag == 1) {
                    Position po;
                    po.rowindex = i;    //r
                    po.colindex = j;    //theta
                    po.votenum = acc[i][j];
                    v.push_back(po);
                }
            }
        }
    }

    return v;
}

vector<Position> Hough::gethighestvote(vector<Position> v) {
    sort(v.begin(), v.end(), cmp);        //按照投票数对直线参数进行排序
    vector<Position>::iterator iter;
    vector<Position>::iterator iter1;
    for (iter = v.begin(); iter != v.end(); iter++) {
        for (iter1 = iter+1; iter1 != v.end(); iter1++) {
            if (abs(iter->rowindex - iter1->rowindex) < 60 && iter->votenum != 0 && iter1->votenum != 0) {
                iter->votenum = iter->votenum + iter1->votenum;
                iter1->votenum = 0;
            }
        }
    }

    sort(v.begin(), v.end(), cmp);
    while(v.size() > 4) {
        v.pop_back();
    }
    return v;
}

void Hough::drawline(vector<Position> v) {
    vector<Position>::iterator iter;
    vector<Position>::iterator iter1;
    for (iter = v.begin(); iter != v.end(); iter++) {        //将投票数最多的直线画在原图像上
        int x1, y1, x2, y2;
        x1 = y1 = x2 = y2 = 0;
        double angle = (double)(iter->colindex) / 180 * PI;
        double si = sin(angle);
        double co = cos(angle);
        if (iter->colindex >= 45 && iter->colindex <= 135) {        //在这个范围内sin值比较大，使用sin做分母误差较小
            x1 = 0;
            y1 = (iter->rowindex - center_h) / si;        //加上之前减去的值才是真正的r
            x2 = edge_w - 1;
            y2 = ((iter->rowindex - center_h) - (double)x2 * co) / si;
        } else {
            y1 = 0;
            x1 = (iter->rowindex - center_h) / co;
            y2 = edge_h - 1;
            x2 = ((iter->rowindex - center_h) - (double)y2 * si) / co;
        }
        const unsigned char color[] = {255, 0, 0};
        img.draw_line(x1, y1, x2, y2, color);
    }
}

void Hough::drawpoint(vector<Position> v) {
    vector<Position>::iterator iter;
    vector<Position>::iterator iter1;
    for (iter = v.begin(); iter != v.end(); iter++) {        //求出直线间的交点
        double angle0 = (double)(iter->colindex) / 180 * PI;
        double si0 = sin(angle0);
        double co0 = cos(angle0);
        int r0 = iter->rowindex - center_h;
        for (iter1 = v.begin(); iter1 != v.end(); iter1++) {
            if (iter == iter1) continue;
            double angle1 = (double)(iter1->colindex) / 180 * PI;
            double si1 = sin(angle1);
            double co1 = cos(angle1);
            int r1 = iter1->rowindex - center_h;

            int deta = iter->colindex - iter1->colindex;
            if (abs(deta) < 30 || abs(deta) > 150) continue;        //两直线夹角太小不是我们要求的纸的顶点
            double detaangle = (double)deta / 180 * PI;
            int x = (double(r1)*si0 - double(r0)*si1) / sin(detaangle);
            int y;
            if (iter->colindex >= 30 && iter->colindex <= 150) {
                y = ((double)r0 - double(x)*co0) / si0;
            } else {
                y = ((double)r1 - (double)(x)*co1) / si1; 
            }
            if (x < edge_w && x >= 0 && y < edge_h && y > 0) {
                int flag = 0;
                int deta = 5;
                for (int i = x - deta; i <= x + deta; i++)
                    for (int j = y - deta; j <= y + deta; j++)
                        if (img(i, j ,0) == 0 && img(i, j, 1) == 0 && img(i, j, 2) == 255) {
                            flag = 1;
                            break;
                        }
                    if (flag == 1) break;
                if (flag == 0) {
                    const unsigned char color[] = {0, 0, 255};
                    img.draw_point(x, y, color);
                    Point p(x, y);
                    vertex.push_back(p);
                }
            }
        }
    }
}