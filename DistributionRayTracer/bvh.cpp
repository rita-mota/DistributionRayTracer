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
	build_recursive(0, objects.size(), root, 1); // -> root node takes all the objects
}

void BVH::printNodes() {
	for (size_t i = 0; i < nodes.size(); ++i) {
		BVHNode* node = nodes[i];
		AABB aabb = node->getAABB();

		std::cout << "Node[" << i << "] ";
		if (node->isLeaf()) {
			std::cout << "LEAF";
		}
		else {
			std::cout << "INNER";
		}
		std::cout << " | AABB min: " << aabb.min << " max: " << aabb.max << std::endl;
	}
}

void BVH::build_recursive(int left_index, int right_index, BVHNode* node, int depth) {
	const int LEAF_THRESHOLD = 2;                      // Minimum number of objects per leaf
	int n_objects = right_index - left_index;          // Number of objects in this node's range

	// === BASE CASE: Create a leaf node if there are few enough objects ===
	if (n_objects <= LEAF_THRESHOLD) {
		node->makeLeaf(left_index, n_objects);         // Store start index and count in the leaf
		return;                                        // Stop recursion
	}

	// === 1. Determine split axis: choose longest AABB extent ===
	AABB box = node->getAABB();                        // Bounding box of the current node
	int axis = 0;                                      // Default to X-axis
	if (box.max.y > box.max.x) axis = 1;               // Choose Y if taller than wide
	if (box.max.z > box.max.getAxisValue(axis)) axis = 2; // Choose Z if deeper than current max

	// === 2. Sort the object list by centroid position along chosen axis ===
	Comparator cmp{};                                  // Custom comparator for sorting
	cmp.dimension = axis;
	std::sort(objects.begin() + left_index, objects.begin() + right_index, cmp);

	// === 3. Compute the midpoint of the AABB to use as a spatial split ===
	float midpoint = (box.min.getAxisValue(axis) + box.max.getAxisValue(axis)) / 2.0f;

	// === 4. Binary search for the first object beyond the midpoint ===
	int low = left_index;
	int high = right_index - 1;
	int split_index = -1;

	while (low <= high) {
		int mid = (low + high) / 2;
		float val = objects[mid]->getCentroid().getAxisValue(axis);
		if (val > midpoint) {
			split_index = mid;
			high = mid - 1;
		}
		else {
			low = mid + 1;
		}
	}

	// === Fallback: if no valid split, use median split ===
	if (split_index <= left_index || split_index >= right_index) {
		split_index = (left_index + right_index) / 2;
	}

	// === 5. Allocate and prepare two new child BVH nodes ===
	BVHNode* left_node = new BVHNode();
	BVHNode* right_node = new BVHNode();
	int left_node_index = nodes.size();                // Index of left node in flat array

	node->makeNode(left_node_index);                   // Store index of first child in current node

	// === 6. Compute bounding boxes for each child node ===
	AABB left_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));
	AABB right_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));

	for (int i = left_index; i < split_index; ++i)
		left_bbox.extend(objects[i]->GetBoundingBox());  // Expand left AABB to fit objects

	for (int i = split_index; i < right_index; ++i)
		right_bbox.extend(objects[i]->GetBoundingBox()); // Expand right AABB to fit objects

	left_node->setAABB(left_bbox);                     // Assign AABBs to new nodes
	right_node->setAABB(right_bbox);
	nodes.push_back(left_node);                        // Add new nodes to global list
	nodes.push_back(right_node);

	// === 7. Recurse into left and right children ===
	build_recursive(left_index, split_index, left_node, depth + 1);     // Left subtree
	build_recursive(split_index, right_index, right_node, depth + 1);   // Right subtree
}




bool BVH::Traverse(Ray& ray, Object** hit_obj, HitRecord& hitRec) {
	float tmp;
	bool hit = false;
	stack<StackItem> hit_stack;
	HitRecord rec;   //rec.isHit initialized to false and rec.t initialized with FLT_MAX
	BVHNode* currentNode = nodes[0];

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
				if (tmpL < tmpR) {
					currentNode = leftNode;
					hit_stack.push(StackItem(rightNode, tmpR));
				}
				else {
					currentNode = rightNode;
					hit_stack.push(StackItem(leftNode, tmpL));
						
				}
				continue;
			}
			else {
				if (leftHit) {
					currentNode = leftNode;
					continue;
				}
				if (rightHit) {
					currentNode = rightNode;
					continue;
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
					*hit_obj = obj;
					hit = true;
				}
			}
		}

		bool hasBetter = false;

		while (!hit_stack.empty()) {
			StackItem stack = hit_stack.top();
			hit_stack.pop();
			if (stack.t < hitRec.t) {
				currentNode = stack.ptr;
				hasBetter = true;
				break;
			}
		}

		if (!hasBetter)
			break;
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
				continue;
			}
			else {
				if (leftHit) {
					currentNode = leftNode;
					continue;
				}
				if (rightHit) {
					currentNode = rightNode;
					continue;
							
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
