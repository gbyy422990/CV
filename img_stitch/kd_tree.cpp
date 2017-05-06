#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <cmath>

using namespace std;

int kd_split;

struct kd_node{
	vector<float> node_data;
	int split;
	int dim_value;
	kd_node *Left;
	kd_node *Right;
	kd_node *Parent;
};

bool cmp(vector<float> a, vector<float> b) {
	return a[kd_split] < b[kd_split];
}

kd_node* createKDNode(vector< vector<float> > &dataset, int begin, int end, kd_node* parent) {
	
	if (begin == end){
		return NULL;
	} else if (end == begin+1) {
		kd_node *newNode = new kd_node;
		vector< vector<float> >::iterator iter = dataset.begin();
		for (int i = 0; i < begin; i++)
			iter++;
		newNode->node_data = (*iter);
		newNode->split = 0;
		newNode->Left = NULL;
		newNode->Right = NULL;
		newNode->Parent = parent;
		return newNode;
	} else {
		double mean[128] = {0.0};
		double var[128] = {0.0};
		int count = 0;

		vector< vector<float> >::iterator beginiter;
		vector< vector<float> >::iterator enditer;
		vector< vector<float> >::iterator iter = dataset.begin();
		for (int i = 0; i < begin; i++)
			iter++;
		beginiter = iter;
		for (int i = begin; i < end; i++)
			iter++;
		enditer = iter;

		for (iter = beginiter; iter != enditer; iter++) {
			for (int i = 0; i < 128; i++) {
				mean[i] += (*iter)[i];
			}
			count++;
		}
		for (int i = 0; i < 128; i++) {
			mean[i] /= count;
		}
		for (iter = beginiter; iter != enditer; iter++) {
			for (int i = 0; i < 128; i++) {
				var[i] += ((*iter)[i] - mean[i]) * ((*iter)[i] - mean[i]);
			}
		}
		double maxvar = -1;
		for (int i = 0; i < 128; i++) {
			if (var[i] > maxvar) {
				maxvar = var[i];
				kd_split = i;
			}
		}


		sort(beginiter, enditer, cmp);

		int size = end - begin;
		int mid = size / 2;
		iter = dataset.begin();
		for (int i = 0; i < begin; i++)
			iter++;
		for (int i = 0; i < mid; i++)
			iter++;
		vector< vector<float> >::iterator miditer = iter;

		kd_node *newNode = new kd_node;
		newNode->node_data = *miditer;
		newNode->split = kd_split;
		newNode->Parent = parent;
		int midindex = mid + begin;
		newNode->Left = createKDNode(dataset, begin, midindex, newNode);
		newNode->Right = createKDNode(dataset, midindex+1, end, newNode);

		return newNode;
	}

}

double Dist(vector<float> v1, vector<float> v2) {
	int size = v1.size();
	double dist = 0.0;
	for (int i = 0; i < size - 1; i++) {
		dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	dist = sqrt(dist);
	return dist;
}

void nearestSearch(kd_node* root, vector<float> target, vector<float> &n1, double &min_dist, double &second_min_dist) {
	list<kd_node*> priority_queue;
	list<kd_node*>::iterator iter;
	vector<float> current_data;
	priority_queue.push_back(root);

    int timeout = 0;
    min_dist = Dist(target, root->node_data);
    second_min_dist = min_dist;
    n1 = root->node_data;
	while (!priority_queue.empty()) {
		if (timeout == 200) break;
		timeout++;
		kd_node* current_node = priority_queue.front();
		priority_queue.pop_front();
		while (current_node) {
			current_data = current_node->node_data;
			int spl = current_node->split;
			if (target[spl] < current_node->node_data[spl]) {
				kd_node* que_node = current_node->Right;
				if (que_node != NULL) {
					int tmp_sql = que_node->split;
					que_node->dim_value = abs(target[tmp_sql] - que_node->node_data[tmp_sql]);
	                int list_size = priority_queue.size();
	                if (list_size > 0) {
						for (iter = priority_queue.begin(); iter != priority_queue.end(); iter++) {
							if (que_node->dim_value < (*iter)->dim_value) {
								priority_queue.insert(iter, que_node);
								break;
							}
						}
						if (list_size == priority_queue.size()) {
							priority_queue.push_back(que_node);
						}
	                } else {
	                	priority_queue.push_back(que_node);
	                }
				}
				current_node = current_node->Left;
			} else {
				kd_node* que_node = current_node->Left;
				if (que_node != NULL) {
					int tmp_sql = que_node->split;
					que_node->dim_value = abs(target[tmp_sql] - que_node->node_data[tmp_sql]);
	                int list_size = priority_queue.size();
	                if (list_size > 0) {
						for (iter = priority_queue.begin(); iter != priority_queue.end(); iter++) {
							if (que_node->dim_value < (*iter)->dim_value) {
								priority_queue.insert(iter, que_node);
								break;
							}
						}
						if (list_size == priority_queue.size()) {
							priority_queue.push_back(que_node);
						}
	                } else {
	                	priority_queue.push_back(que_node);
	                }
				}
				current_node = current_node->Right;
			}
			double cur_dist = Dist(target, current_data);
			if (min_dist == second_min_dist && cur_dist >= min_dist) {
				second_min_dist = cur_dist;
			}
			if (cur_dist < min_dist) {
				n1 = current_data;
				min_dist = cur_dist;
			}
			if (cur_dist < second_min_dist && cur_dist >= min_dist) {
				second_min_dist = cur_dist;
			}
		}
	}
}
