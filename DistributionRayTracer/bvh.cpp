#include "rayAccelerator.h"
#include "macros.h"

using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }


void BVH::Build(vector<Object *> &objs) {

		
			BVHNode *root = new BVHNode();

			Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			AABB world_bbox = AABB(min, max);

			for (Object* obj : objs) {
				AABB bbox = obj->GetBoundingBox();
				world_bbox.extend(bbox);
				objects.push_back(obj);
			}
			world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
			world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
			root->setAABB(world_bbox);
			nodes.push_back(root);
			build_recursive(0, objects.size(), root); // -> root node takes all the objects
		}

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {

	//PUT YOUR CODE HERE

		//right_index, left_index and split_index refer to the indices in the objects vector
	   // do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	    // node.index can have a index of objects vector or a index of nodes vector
			
	}

bool BVH::Traverse(Ray& ray, Object** hit_obj, HitRecord& hitRec) {
			float tmp;
			bool hit = false;
			stack<StackItem> hit_stack;
			HitRecord rec;   //rec.isHit initialized to false and rec.t initialized with FLT_MAX
			BVHNode* currentNode = nodes[0];
			//PUT YOUR CODE HERE

			Ray localRay = ray;
			hitRec = rec;

			if (!currentNode->getAABB().hit(localRay, tmp))
				return false;

			while (true) {
				if (!currentNode->isLeaf()) {
					int leftIndex = currentNode->getIndex();
					BVHNode* leftNode = nodes[leftIndex];
					BVHNode* rightNode = nodes[leftIndex+1];

					float tmpL;
					float tmpR;
					bool leftHit = leftNode->getAABB().hit(localRay, tmpL);
					bool rightHit = rightNode->getAABB().hit(localRay, tmpR);

					if (leftNode->getAABB().isInside(ray.origin))
						tmpL = 0;
					if (rightNode->getAABB().isInside(ray.origin))
						tmpR = 0;

					if (leftHit && rightHit) {
						if (tmpL <= tmpR) {
							currentNode = leftNode;
							//push right to stash
							hit_stack.push(StackItem(rightNode, tmpR));
						}
						else {
							currentNode = rightNode;
							//push left to stash
							hit_stack.push(StackItem(leftNode, tmpL));
						}
					}
					else {
						if (leftHit) {
							currentNode = leftNode;
						}
						if (rightHit) {
							currentNode = rightNode;
						}
					}
				}
				else {
					int nObj = currentNode->getNObjs();
					int objIndex = currentNode->getIndex();
					Object* obj;
					for (int i = 0; i < nObj; i++) {
						obj = objects[objIndex + i];
						rec = obj->hit(localRay);
						if (rec.isHit && rec.t < hitRec.t) {	
							hitRec = rec;
							hit_obj = &obj;
							hit = true;
						}
					}
				}

			while (!hit_stack.empty()) {
				StackItem stack = hit_stack.top();
				hit_stack.pop();
				if (stack.t < hitRec.t) {
					currentNode = stack.ptr;
					break;
				}
			}
			if (hit_stack.empty())
				return hit;

			}
		return hit;
	}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
			float tmp;
			stack<StackItem> hit_stack;
			HitRecord rec;

			double length = ray.direction.length(); //distance between light and intersection point
			ray.direction.normalize();

			Ray localRay = ray;
			BVHNode* currentNode = nodes[0];

			if (!currentNode->getAABB().hit(localRay, tmp))
				return false;

			while (true) {
				if (!currentNode->isLeaf()) {
					int leftIndex = currentNode->getIndex();
					BVHNode* leftNode = nodes[leftIndex];
					BVHNode* rightNode = nodes[leftIndex + 1];

					float tmpL;
					float tmpR;
					bool leftHit = leftNode->getAABB().hit(localRay, tmpL);
					bool rightHit = rightNode->getAABB().hit(localRay, tmpR);

					if (leftHit && rightHit) {
						if (tmpL <= tmpR) {
							currentNode = leftNode;
							//push right to stash
							hit_stack.push(StackItem(rightNode, tmpR));
						}
						else {
							currentNode = rightNode;
							//push left to stash
							hit_stack.push(StackItem(leftNode, tmpL));
						}
					}
					else {
						if (leftHit) {
							currentNode = leftNode;
						}
						if (rightHit) {
							currentNode = rightNode;
						}
					}
				}
				else {
					int nObj = currentNode->getNObjs();
					int objIndex = currentNode->getIndex();
					Object* obj;
					for (int i = 0; i < nObj; i++) {
						obj = objects[objIndex + i];
						rec = obj->hit(localRay);
						if (rec.isHit)
							return true;
					}
				}

				if (hit_stack.empty())
					return false;

				StackItem stack = hit_stack.top();
				currentNode = stack.ptr;
				hit_stack.pop();
			}

			return false;  //no primitive intersection		
	}		
