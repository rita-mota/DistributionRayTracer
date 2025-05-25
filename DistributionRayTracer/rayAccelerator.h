#ifndef ACCELERATOR_H
#define ACCELERATOR_H

#include <stack>
#include <queue>
#include <cmath>
#include <algorithm>
#include "scene.h"

using namespace std;

class Grid
{
public:
	Grid(void);
	//~Grid(void);
	int getNumObjects();
	void addObject(Object* o);
	void setAABB(AABB& bbox_);
	Object* getObject(unsigned int index);
	void Build(vector<Object*>& objs);   // set up grid cells
	bool Traverse(Ray& ray, Object **hitobject, HitRecord& hitRec);  
	bool Traverse(Ray& ray);  //Traverse for shadow ray

private:
	vector<Object *> objects;
	vector<vector<Object*> > cells;

	int nx, ny, nz; // number of cells in the x, y, and z directions
	float m = 2.0f; // factor that allows to vary the number of cells

	//Setup function for Grid traversal
	bool Init_Traverse(Ray& ray, int& ix, int& iy, int& iz, double& dtx, double& dty, double& dtz, double& tx_next, double& ty_next, double& tz_next, 
		int& ix_step, int& iy_step, int& iz_step, int& ix_stop, int& iy_stop, int& iz_stop);

	AABB bbox;
};

/*********************************BVH*****************************************************************/
class BVH
{
	class Comparator {
	public:
		int dimension;

		bool operator() (Object* a, Object* b) {
			float ca = a->getCentroid().getAxisValue(dimension);
			float cb = b->getCentroid().getAxisValue(dimension);
			return ca < cb;
		}
	};

	class BVHNode {
	private:
		AABB bbox;
		bool leaf;
		unsigned int n_objs;
		unsigned int index;	// if leaf == false: index to left child node,
							// else if leaf == true: index to first Intersectable (Object *) in objects vector

	public:
		BVHNode(void);
		void setAABB(AABB& bbox_);
		void makeLeaf(unsigned int index_, unsigned int n_objs_);
		void makeNode(unsigned int left_index_);
		bool isLeaf() { return leaf; }
		unsigned int getIndex() { return index; }
		unsigned int getNObjs() { return n_objs; }
		AABB& getAABB() { return bbox; };
	};

private:
	int Threshold = 2;
	vector<Object*> objects;
	vector<BVH::BVHNode*> nodes;

	struct StackItem {
		BVHNode* ptr;
		float t;
		StackItem(BVHNode* _ptr, float _t) : ptr(_ptr), t(_t) { }
	};

	//stack<StackItem> hit_stack;  just declare it in the traverse procedure in order to be parallelized with OMP

public:
	BVH(void);
	int getNumObjects();
	
	void Build(vector<Object*>& objects);
	void build_recursive(int left_index, int right_index, BVHNode* node, int depth);
	bool Traverse(Ray& ray, Object** hit_obj, HitRecord& hitRec);
	bool Traverse(Ray& ray);
	void printNodes();
};
#endif