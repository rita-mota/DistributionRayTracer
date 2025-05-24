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

void BVH::build_recursive(int left_index, int right_index, BVHNode* node) {
    const int LEAF_THRESHOLD = 2; // If the number of objects is less than or equal to this, we make a leaf
    int n_objects = right_index - left_index; // Number of objects in this range

    // ====== 1. BASE CASE: Make a LEAF node if we have few enough objects ======
    if (n_objects <= LEAF_THRESHOLD) {
        node->makeLeaf(left_index, n_objects);
        return;
    }

    // ====== 2. Compute bounding box of object CENTROIDS to determine spatial distribution ======
    Vector min_c(FLT_MAX, FLT_MAX, FLT_MAX);     // Initialize min corner to +inf
    Vector max_c(-FLT_MAX, -FLT_MAX, -FLT_MAX);  // Initialize max corner to -inf

    for (int i = left_index; i < right_index; ++i) {
        Vector c = objects[i]->GetBoundingBox().centroid();
        // Expand bounding box for centroid bounds
        if (c.x < min_c.x) min_c.x = c.x;
        if (c.y < min_c.y) min_c.y = c.y;
        if (c.z < min_c.z) min_c.z = c.z;
        if (c.x > max_c.x) max_c.x = c.x;
        if (c.y > max_c.y) max_c.y = c.y;
        if (c.z > max_c.z) max_c.z = c.z;
    }

    // ====== 3. Choose SPLIT AXIS (x, y, or z) by largest extent in centroid bounds ======
    Vector extent = max_c - min_c;
    int axis = 0; // Start with x-axis
    if (extent.y > extent.x) axis = 1;
    if (extent.z > extent.getAxisValue(axis)) axis = 2;

    // If extent is too small (all centroids are effectively equal), create leaf
    if (extent.getAxisValue(axis) < 1e-5f) {
        node->makeLeaf(left_index, n_objects);
        return;
    }

    // ====== 4. Sort objects in-place by their centroid along the chosen axis ======
    auto comparator = [axis](Object* a, Object* b) {
        return a->GetBoundingBox().centroid().getAxisValue(axis) <
            b->GetBoundingBox().centroid().getAxisValue(axis);
    };
    std::sort(objects.begin() + left_index, objects.begin() + right_index, comparator);

    // ====== 5. Find the midpoint for splitting the objects ======
    int split_index = (left_index + right_index) / 2;
    if (split_index == left_index || split_index == right_index) {
        // Fallback to a balanced split if midpoint failed
        split_index = left_index + (n_objects / 2);
    }

    // ====== 6. Create child BVH nodes and store them ======
    BVHNode* left_node = new BVHNode();
    BVHNode* right_node = new BVHNode();

    int left_node_index = (int)nodes.size(); // Index of left node in 'nodes' vector
    nodes.push_back(left_node);
    nodes.push_back(right_node); // Right node implicitly follows

    node->makeNode(left_node_index); // Register this node as internal with children

    // ====== 7. Compute and assign AABBs (Axis-Aligned Bounding Boxes) for both children ======
    AABB left_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));
    AABB right_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));

    for (int i = left_index; i < split_index; ++i)
        left_bbox.extend(objects[i]->GetBoundingBox());

    for (int i = split_index; i < right_index; ++i)
        right_bbox.extend(objects[i]->GetBoundingBox());

    left_node->setAABB(left_bbox);
    right_node->setAABB(right_bbox);

    // ====== 8. RECURSIVELY build child nodes ======
    build_recursive(left_index, split_index, left_node);   // Left subtree
    build_recursive(split_index, right_index, right_node); // Right subtree
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
