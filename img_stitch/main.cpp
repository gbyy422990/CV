extern "C" {
    #include <vl/generic.h>
    #include <vl/sift.h>
}

#include <cstdlib>
#include <ctime>
#include "CImg.h"
#include "kd_tree.cpp"
using namespace cimg_library;

struct match_pair {
	vector<float> orginal;
	vector<float> target;
};

void get_sift(CImg<unsigned char> img, int noctaves, int nlevels, int o_min, vector< vector<float> > &descris, vector< vector<int> > &areas) {

    int width = img.width();
	int height = img.height();

    CImg<unsigned char> gray(width, height, 1, 1);
    cimg_forXY(gray, x, y) {
        int newValue = (img(x,y,0) * 0.2126 + img(x,y,1) * 0.7152 + img(x,y,2) * 0.0722);
        gray(x, y) = newValue;
    }

	vl_sift_pix *vlimg = new vl_sift_pix[width * height];
    cimg_forXY(gray, x, y) {
    	int index = x + width * y;
    	vlimg[index] = (float)gray(x, y);
    }

	VlSiftFilt * filter = vl_sift_new(width, height, noctaves, nlevels, o_min);

	if (vl_sift_process_first_octave(filter, vlimg) != VL_ERR_EOF) {
		while (1) {
			vl_sift_detect(filter);
			VlSiftKeypoint *pkeypoint = filter->keys;

			for (int i = 0; i < filter->nkeys; i++) {
				VlSiftKeypoint one_point = *pkeypoint;
				pkeypoint++;

				vector<int> one_area;
				one_area.push_back(one_point.x);
				one_area.push_back(one_point.y);
				one_area.push_back(one_point.sigma/2);

				double angles[4];

				int angleCount = vl_sift_calc_keypoint_orientations(filter, angles, &one_point);
				one_area.push_back(angleCount);
				areas.push_back(one_area);

				for (int j = 0; j < angleCount; j++) {
					float *descriptor = new float[128];    //最后一维代表该描述子在areas中的下标
					vl_sift_calc_keypoint_descriptor(filter, descriptor, &one_point, angles[j]);

					vector<float> one_descri;
					for (int k = 0; k < 128; k++) {
						one_descri.push_back(descriptor[k]);
					}
					one_descri.push_back((float)i);
					descris.push_back(one_descri);
					delete [] descriptor;
					descriptor = NULL;
				}
			}

			if (vl_sift_process_next_octave(filter) == VL_ERR_EOF) {
				break;
			}
		}
	}
	vl_sift_delete(filter);
    delete [] vlimg;
    vlimg = NULL;

}

CImg<float> getParam(int ind1, int ind2, int ind3, int ind4, int is_reserve, vector<match_pair> Pairs, vector< vector<int> >left_areas, vector< vector<int> >right_areas) {
    
    int org_ind1, org_ind2, org_ind3, org_ind4, tar_ind1, tar_ind2, tar_ind3, tar_ind4;
    if (is_reserve) {
        tar_ind1 = Pairs[ind1].orginal[128];
        tar_ind2 = Pairs[ind2].orginal[128];
        tar_ind3 = Pairs[ind3].orginal[128];
        tar_ind4 = Pairs[ind4].orginal[128];
        org_ind1 = Pairs[ind1].target[128];
        org_ind2 = Pairs[ind2].target[128];
        org_ind3 = Pairs[ind3].target[128];
        org_ind4 = Pairs[ind4].target[128];
    } else {
        org_ind1 = Pairs[ind1].orginal[128];
        org_ind2 = Pairs[ind2].orginal[128];
        org_ind3 = Pairs[ind3].orginal[128];
        org_ind4 = Pairs[ind4].orginal[128];
        tar_ind1 = Pairs[ind1].target[128];
        tar_ind2 = Pairs[ind2].target[128];
        tar_ind3 = Pairs[ind3].target[128];
        tar_ind4 = Pairs[ind4].target[128];
    }

    int x[4];
    int y[4];
    int u[4];
    int v[4];

    x[0] = right_areas[org_ind1][0];
    y[0] = right_areas[org_ind1][1];
    x[1] = right_areas[org_ind2][0];
    y[1] = right_areas[org_ind2][1];
    x[2] = right_areas[org_ind3][0];
    y[2] = right_areas[org_ind3][1];
    x[3] = right_areas[org_ind4][0];
    y[3] = right_areas[org_ind4][1];

    u[0] = left_areas[tar_ind1][0];
    v[0] = left_areas[tar_ind1][1];
    u[1] = left_areas[tar_ind2][0];
    v[1] = left_areas[tar_ind2][1];
    u[2] = left_areas[tar_ind3][0];
    v[2] = left_areas[tar_ind3][1];
    u[3] = left_areas[tar_ind4][0];
    v[3] = left_areas[tar_ind4][1];

    CImg<float> matrix1 = CImg<float>(8, 8, 1, 1);
    CImg<float> matrix2 = CImg<float>(1, 8, 1, 1);
    CImg<float> param = CImg<float>(1, 8, 1, 1);

    for (int i = 0; i < 8; i++) {
        if (i % 2 == 0) {
            matrix1(0, i) = x[i/2];
            matrix1(1, i) = y[i/2];
            matrix1(2, i) = 1;
            matrix1(3, i) = 0;
            matrix1(4, i) = 0;
            matrix1(5, i) = 0;
            matrix1(6, i) = -x[i/2] * u[i/2];
            matrix1(7, i) = -y[i/2] * u[i/2];
        } else {
            matrix1(0, i) = 0;
            matrix1(1, i) = 0;
            matrix1(2, i) = 0;
            matrix1(3, i) = x[i/2];
            matrix1(4, i) = y[i/2];
            matrix1(5, i) = 1;
            matrix1(6, i) = -x[i/2] * v[i/2];
            matrix1(7, i) = -y[i/2] * v[i/2];
        }
    }
    
    matrix2(0, 0) = u[0];
    matrix2(0, 1) = v[0];
    matrix2(0, 2) = u[1];
    matrix2(0, 3) = v[1];
    matrix2(0, 4) = u[2];
    matrix2(0, 5) = v[2];
    matrix2(0, 6) = u[3];
    matrix2(0, 7) = v[3];

    param = matrix2.get_solve(matrix1);
    return param;

}

int main (int argc, const char * argv[]) {

	CImg<unsigned char>left_img(argv[1]);
	CImg<unsigned char>right_img(argv[2]);

    int width1 = left_img.width();
    int height1 = left_img.height();
    int width2 = right_img.width();
    int height2 = right_img.height();

    int width3 = width1 + width2;
    int height3 = height1 > height2 ? height1 : height2;
	CImg<unsigned char>img_join(width3, height3, 1, 3);
	CImg<unsigned char>img_ransac(width3, height3, 1, 3);

	cimg_forXY(img_join, x, y) {
		if (x < width1 && y < height1) {
			img_join(x, y, 0) = left_img(x, y, 0);
			img_join(x, y, 1) = left_img(x, y, 1);
			img_join(x, y, 2) = left_img(x, y, 2);

			img_ransac(x, y, 0) = left_img(x, y, 0);
			img_ransac(x, y, 1) = left_img(x, y, 1);
			img_ransac(x, y, 2) = left_img(x, y, 2);
		} else if (x < width1 && y >= height1) {
			img_join(x, y, 0) = 0;
			img_join(x, y, 1) = 0;
			img_join(x, y, 2) = 0;

			img_ransac(x, y, 0) = 0;
			img_ransac(x, y, 1) = 0;
			img_ransac(x, y, 2) = 0;
		} else if (x > width1 && y < height2) {
			img_join(x, y, 0) = right_img(x - width1, y, 0);
			img_join(x, y, 1) = right_img(x - width1, y, 1);
			img_join(x, y, 2) = right_img(x - width1, y, 2);

			img_ransac(x, y, 0) = right_img(x - width1, y, 0);
			img_ransac(x, y, 1) = right_img(x - width1, y, 1);
			img_ransac(x, y, 2) = right_img(x - width1, y, 2);
		} else {
			img_join(x, y, 0) = 0;
			img_join(x, y, 1) = 0;
			img_join(x, y, 2) = 0;

			img_ransac(x, y, 0) = 0;
			img_ransac(x, y, 1) = 0;
			img_ransac(x, y, 2) = 0;
		}
	}

	vector< vector<float> >left_descris;
	vector< vector<int> >left_areas;
	vector< vector<float> >right_descris;
	vector< vector<int> >right_areas;

    get_sift(left_img, 4, 2, 0, left_descris, left_areas);    //sift特征提取
    get_sift(right_img, 4, 2, 0, right_descris, right_areas);

    int tree_size = left_descris.size();
    kd_node *head = createKDNode(left_descris, 0, tree_size, NULL);    //kd树创建

    int search_size = right_descris.size();
    vector<float> n1;
    double min_dist, second_min_dist;
    vector<match_pair> Pairs;
    for (int i = 0; i < search_size; i++) {
		nearestSearch(head, right_descris[i], n1, min_dist, second_min_dist);    //kd树查询
		double dist_rate = min_dist / second_min_dist;
		if (dist_rate < 0.49) {
			match_pair one_pair;
			one_pair.orginal = right_descris[i];
			one_pair.target = n1;
			Pairs.push_back(one_pair);
			int src_index = (int)right_descris[i][128];
			int tar_index = (int)n1[128];
			int x1 = left_areas[tar_index][0];
			int y1 = left_areas[tar_index][1];
            int x2 = right_areas[src_index][0];
            int y2 = right_areas[src_index][1];
			if (x1 < width1 && y1 < height1 && x2 < width2 && y2 < height2) {
				const unsigned char color[] = {255, 0, 0};
	            img_join.draw_line(x1, y1, x2+width1, y2, color);
			}
		}
    }
    // img_join.display();
    
    int pair_size = Pairs.size();
    int iterator = log(1 - 0.999) / log(1 - pow((1 - 0.4), 4));
    srand((unsigned)time(NULL));
    vector<match_pair> best_set;
    double min_error = 10000000;
    CImg<float> best_param = CImg<float>(1, 8, 1, 1);

    for (int k = 0; k < iterator; k++) {    //RANSAC算法实现
    	int ind1 = rand() % pair_size;
    	int ind2 = rand() % pair_size;
    	while (ind2 == ind1) {
    		ind2 = rand() % pair_size;
    	}
    	int ind3 = rand() % pair_size;
    	while (ind3 == ind1 || ind3 == ind2) {
    		ind3 = rand() % pair_size;
    	}
    	int ind4 = rand() % pair_size;
    	while (ind4 == ind1 || ind4 == ind2 || ind4 == ind3) {
    		ind4 = rand() % pair_size;
    	}

        CImg<float> param = CImg<float>(1, 8, 1, 1);
	    param = getParam(ind1, ind2, ind3, ind4, 0, Pairs, left_areas, right_areas);

        vector<match_pair> current_set;
        double cur_error = 0;
        current_set.push_back(Pairs[ind1]);
        current_set.push_back(Pairs[ind2]);
        current_set.push_back(Pairs[ind3]);
        current_set.push_back(Pairs[ind4]);
        for (int i = 0; i < pair_size; i++) {
        	if (i != ind1 && i != ind2 && i != ind3 && i != ind4) {
        		int org_ind = Pairs[i].orginal[128];
        		int tar_ind = Pairs[i].target[128];

        		int org_x = right_areas[org_ind][0];
        		int org_y = right_areas[org_ind][1];
        		int tar_x = left_areas[tar_ind][0];
        		int tar_y = left_areas[tar_ind][1];

		        float hx = (param(0,0) * org_x + param(0,1) * org_y + param(0,2)) / (param(0,6) * org_x + param(0,7) * org_y + 1);
                float hy = (param(0,3) * org_x + param(0,4) * org_y + param(0,5)) / (param(0,6) * org_x + param(0,7) * org_y + 1);

                double error = sqrt((hx-tar_x)*(hx-tar_x) + (hy-tar_y)*(hy-tar_y));
                if (error < 3.0) {
                	current_set.push_back(Pairs[i]);
                    cur_error += error;
                }
        	}
        }
        int best_size = best_set.size();
        int cur_size = current_set.size();
        if (cur_size > best_size || (cur_size == best_size && cur_error < min_error)) {
        	best_set.clear();
        	best_set = current_set;
            min_error = cur_error;
        	cimg_forXY(best_param, x, y) {
        		best_param(x, y) = param(x, y);
        	}
        }
    }

    int best_size = best_set.size();
    for (int i = 0; i < best_size; i++) {
        int src_index = best_set[i].orginal[128];
        int tar_index = best_set[i].target[128];
        int x1 = left_areas[tar_index][0];
        int y1 = left_areas[tar_index][1];
        int x2 = right_areas[src_index][0];
        int y2 = right_areas[src_index][1];
        if (x1 < width1 && y1 < height1 && x2 < width2 && y2 < height2) {
            const unsigned char color[] = {255, 0, 0};
            img_ransac.draw_line(x1, y1, x2+width1, y2, color);
        }
    }
    // img_ransac.display();

    float rightup_x = (best_param(0,0) * (width2-1) + best_param(0,1) * 0 + best_param(0,2)) / (best_param(0,6) * (width2-1) + best_param(0,7) * 0 + 1);
    float rightdown_x = (best_param(0,0) * (width2-1) + best_param(0,1) * (height2-1) + best_param(0,2)) / (best_param(0,6) * (width2-1) + best_param(0,7) * (height2-1) + 1);

    int stitch_width = rightup_x > rightdown_x ? rightup_x : rightdown_x;
    int stitch_height = height1 < height2 ? height1 : height2;
    CImg<unsigned char> stitch_img(stitch_width, stitch_height, 1, 3);

    CImg<float> param = CImg<float>(1, 8, 1, 1);
    param = getParam(0, 1, 2, 3, 1, best_set, right_areas, left_areas);

    float leftup_x = best_param(0,2);
    float leftdown_x = (best_param(0,1) * (height2-1) + best_param(0,2)) / (best_param(0,7) * (height2-1) + 1);

    int left_bound = leftup_x < leftdown_x ? leftup_x : leftdown_x;
    int right_bound = width1;

    cimg_forXY(stitch_img, x, y) {
        stitch_img(x, y, 0) = 0;
        stitch_img(x, y, 1) = 0;
        stitch_img(x, y, 2) = 0;
    }

    // cimg_forXY(left_img, x, y) {    //拼接图像会有黑点
    //     if (y < stitch_height) {
    //         stitch_img(x, y, 0) = left_img(x, y, 0);
    //         stitch_img(x, y, 1) = left_img(x, y, 1);
    //         stitch_img(x, y, 2) = left_img(x, y, 2);
    //     }
    // }

    // cimg_forXY(right_img, x, y) {
    //     float u = (best_param(0,0) * (float)x + best_param(0,1) * (float)y + best_param(0,2)) / (best_param(0,6) * (float)x + best_param(0,7) * (float)y + 1);
    //     float v = (best_param(0,3) * (float)x + best_param(0,4) * (float)y + best_param(0,5)) / (best_param(0,6) * (float)x + best_param(0,7) * (float)y + 1);

    //     if (u >= 0 && u < stitch_width && v >= 0 && v < stitch_height) {
    //         stitch_img(u, v, 0) = right_img(x, y, 0);
    //         stitch_img(u, v, 1) = right_img(x, y, 1);
    //         stitch_img(u, v, 2) = right_img(x, y, 2);
    //     }
    // }

    // cimg_forXY(stitch_img, x, y) {    //会有明显拼接痕迹
    //     float u = (param(0,0) * (float)x + param(0,1) * (float)y + param(0,2)) / (param(0,6) * (float)x + param(0,7) * (float)y + 1);
    //     float v = (param(0,3) * (float)x + param(0,4) * (float)y + param(0,5)) / (param(0,6) * (float)x + param(0,7) * (float)y + 1);

    //     if (u >= 0 && u < width2 && v >= 0 && v < height2) {
    //         float u_ceil = ceil(u);
    //         float v_ceil = ceil(v);
    //         float u_floor = floor(u);
    //         float v_floor = floor(v);

    //         float a = u - u_floor;
    //         float b = v - v_floor;

    //         stitch_img(x, y, 0) = right_img(u_floor, v_floor, 0)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 0)*a*(1-b) + right_img(u_floor, v_ceil, 0)*(1-a)*b + right_img(u_ceil, v_ceil, 0)*a*b;
    //         stitch_img(x, y, 1) = right_img(u_floor, v_floor, 1)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 1)*a*(1-b) + right_img(u_floor, v_ceil, 1)*(1-a)*b + right_img(u_ceil, v_ceil, 1)*a*b;
    //         stitch_img(x, y, 2) = right_img(u_floor, v_floor, 2)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 2)*a*(1-b) + right_img(u_floor, v_ceil, 2)*(1-a)*b + right_img(u_ceil, v_ceil, 2)*a*b;
    //     }
    // }

    // cimg_forXY(stitch_img, x, y) {
    //     if (stitch_img(x, y, 0) == 0 && stitch_img(x, y, 1) == 0 && stitch_img(x, y, 2) == 0) {
    //         if (x < width1 && y < height1) {
    //             stitch_img(x, y, 0) = left_img(x, y, 0);
    //             stitch_img(x, y, 1) = left_img(x, y, 1);
    //             stitch_img(x, y, 2) = left_img(x, y, 2);
    //         }
    //     }
    // }

    cimg_forXY(stitch_img, x, y) {
        float u = (param(0,0) * (float)x + param(0,1) * (float)y + param(0,2)) / (param(0,6) * (float)x + param(0,7) * (float)y + 1);
        float v = (param(0,3) * (float)x + param(0,4) * (float)y + param(0,5)) / (param(0,6) * (float)x + param(0,7) * (float)y + 1);

        if (u >= 0 && u < width2 && v >= 0 && v < height2) {
            float u_ceil = ceil(u);
            float v_ceil = ceil(v);
            float u_floor = floor(u);
            float v_floor = floor(v);

            float a = u - u_floor;
            float b = v - v_floor;

            float stitch_r = right_img(u_floor, v_floor, 0)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 0)*a*(1-b) + right_img(u_floor, v_ceil, 0)*(1-a)*b + right_img(u_ceil, v_ceil, 0)*a*b;
            float stitch_g = right_img(u_floor, v_floor, 1)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 1)*a*(1-b) + right_img(u_floor, v_ceil, 1)*(1-a)*b + right_img(u_ceil, v_ceil, 1)*a*b;
            float stitch_b = right_img(u_floor, v_floor, 2)*(1-a)*(1-b) + right_img(u_ceil, v_floor, 2)*a*(1-b) + right_img(u_floor, v_ceil, 2)*(1-a)*b + right_img(u_ceil, v_ceil, 2)*a*b;            
        
            if (x >= left_bound && x <= right_bound) {
                int dl = x - left_bound;
                int dr = right_bound - x;
                stitch_img(x, y, 0) = (float)dr/(dr+dl) * left_img(x, y, 0) + (float)dl/(dr+dl) * stitch_r;
                stitch_img(x, y, 1) = (float)dr/(dr+dl) * left_img(x, y, 1) + (float)dl/(dr+dl) * stitch_g;
                stitch_img(x, y, 2) = (float)dr/(dr+dl) * left_img(x, y, 2) + (float)dl/(dr+dl) * stitch_b;
            } else {
                stitch_img(x, y, 0) = stitch_r;
                stitch_img(x, y, 1) = stitch_g;
                stitch_img(x, y, 2) = stitch_b;
            }
        }
    }

    cimg_forXY(stitch_img, x, y) {
        if (stitch_img(x, y, 0) == 0 && stitch_img(x, y, 1) == 0 && stitch_img(x, y, 2) == 0) {
            if (x < width1 && y < height1) {
                stitch_img(x, y, 0) = left_img(x, y, 0);
                stitch_img(x, y, 1) = left_img(x, y, 1);
                stitch_img(x, y, 2) = left_img(x, y, 2);
            }
        }
    }
    stitch_img.display();
    // stitch_img.save_bmp("img1.bmp");

	return 0;
}
