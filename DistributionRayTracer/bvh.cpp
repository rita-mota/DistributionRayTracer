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


void BVH::Build(vector<Object*>& objs) {

	BVHNode* root = new BVHNode();

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
	const int LEAF_THRESHOLD = 2; // Minimum number of objects per leaf
	const int BUCKET_COUNT = 12; // Number of bins for SAH
	const float TRAVERSAL_COST = 1.0f;
	const float INTERSECTION_COST = 1.0f;

	int n_objects = right_index - left_index; // Number of objects in this node's range

	// === BASE CASE: Create a leaf node if there are few enough objects ===
	if (n_objects <= LEAF_THRESHOLD) {
		node->makeLeaf(left_index, n_objects);
		return;
	}

	AABB box = node->getAABB();
	// Calculate parent surface area directly
	Vector extent = box.max - box.min;
	float parent_surface_area = 2.0f * (extent.x * extent.y + extent.x * extent.z + extent.y * extent.z);

	int best_axis = 0;
	float best_cost = FLT_MAX;
	int best_split = left_index;

	// Try all 3 axes
	for (int axis = 0; axis < 3; axis++) {
		// Sort objects along current axis
		Comparator cmp{};
		cmp.dimension = axis;
		std::sort(objects.begin() + left_index, objects.begin() + right_index, cmp);

		// Initialize buckets
		struct Bucket {
			int count = 0;
			Vector min_bounds = Vector(FLT_MAX, FLT_MAX, FLT_MAX);
			Vector max_bounds = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		};
		Bucket buckets[BUCKET_COUNT];

		// Determine range and scale factor for bucketing
		float min_bound = box.min.getAxisValue(axis);
		float max_bound = box.max.getAxisValue(axis);
		float scale = (max_bound - min_bound) > 0.0f ?
			BUCKET_COUNT / (max_bound - min_bound) : 0.0f;

		// Place objects into buckets
		for (int i = left_index; i < right_index; i++) {
			float centroid = objects[i]->getCentroid().getAxisValue(axis);
			int bucket_idx = std::min(BUCKET_COUNT - 1,
				(int)((centroid - min_bound) * scale));
			buckets[bucket_idx].count++;

			// Update bucket bounds
			AABB obj_box = objects[i]->GetBoundingBox();
			buckets[bucket_idx].min_bounds = Vector(
				std::min(buckets[bucket_idx].min_bounds.x, obj_box.min.x),
				std::min(buckets[bucket_idx].min_bounds.y, obj_box.min.y),
				std::min(buckets[bucket_idx].min_bounds.z, obj_box.min.z)
			);
			buckets[bucket_idx].max_bounds = Vector(
				std::max(buckets[bucket_idx].max_bounds.x, obj_box.max.x),
				std::max(buckets[bucket_idx].max_bounds.y, obj_box.max.y),
				std::max(buckets[bucket_idx].max_bounds.z, obj_box.max.z)
			);
		}

		// Evaluate all possible splits
		for (int i = 1; i < BUCKET_COUNT; i++) {
			Vector left_min(FLT_MAX, FLT_MAX, FLT_MAX);
			Vector left_max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			Vector right_min(FLT_MAX, FLT_MAX, FLT_MAX);
			Vector right_max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			int left_count = 0, right_count = 0;

			// Compute left bucket properties
			for (int j = 0; j < i; j++) {
				left_min = Vector(
					std::min(left_min.x, buckets[j].min_bounds.x),
					std::min(left_min.y, buckets[j].min_bounds.y),
					std::min(left_min.z, buckets[j].min_bounds.z)
				);
				left_max = Vector(
					std::max(left_max.x, buckets[j].max_bounds.x),
					std::max(left_max.y, buckets[j].max_bounds.y),
					std::max(left_max.z, buckets[j].max_bounds.z)
				);
				left_count += buckets[j].count;
			}

			// Compute right bucket properties
			for (int j = i; j < BUCKET_COUNT; j++) {
				right_min = Vector(
					std::min(right_min.x, buckets[j].min_bounds.x),
					std::min(right_min.y, buckets[j].min_bounds.y),
					std::min(right_min.z, buckets[j].min_bounds.z)
				);
				right_max = Vector(
					std::max(right_max.x, buckets[j].max_bounds.x),
					std::max(right_max.y, buckets[j].max_bounds.y),
					std::max(right_max.z, buckets[j].max_bounds.z)
				);
				right_count += buckets[j].count;
			}

			// Calculate surface areas directly
			Vector left_extent = left_max - left_min;
			float left_area = 2.0f * (left_extent.x * left_extent.y +
				left_extent.x * left_extent.z +
				left_extent.y * left_extent.z);

			Vector right_extent = right_max - right_min;
			float right_area = 2.0f * (right_extent.x * right_extent.y +
				right_extent.x * right_extent.z +
				right_extent.y * right_extent.z);

			// Compute SAH cost
			float cost = TRAVERSAL_COST +
				(left_count * left_area + right_count * right_area) *
				INTERSECTION_COST / parent_surface_area;

			if (cost < best_cost) {
				best_cost = cost;
				best_axis = axis;
				best_split = left_index;

				// Find the split index corresponding to this bucket
				int count = 0;
				for (int j = 0; j < i; j++) {
					count += buckets[j].count;
				}
				best_split = left_index + count;
			}
		}
	}

	// If no good split found, create leaf
	if (best_split <= left_index || best_split >= right_index || best_cost >= n_objects * INTERSECTION_COST) {
		node->makeLeaf(left_index, n_objects);
		return;
	}

	// Sort along best axis
	Comparator cmp{};
	cmp.dimension = best_axis;
	std::sort(objects.begin() + left_index, objects.begin() + right_index, cmp);

	BVHNode* left_node = new BVHNode();
	BVHNode* right_node = new BVHNode();
	int left_node_index = nodes.size();
	node->makeNode(left_node_index);

	AABB left_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));
	AABB right_bbox(Vector(FLT_MAX, FLT_MAX, FLT_MAX), Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX));

	for (int i = left_index; i < best_split; i++) {
		left_bbox.extend(objects[i]->GetBoundingBox());  // Expand left AABB to fit objects
	}
	for (int i = best_split; i < right_index; i++) {
		right_bbox.extend(objects[i]->GetBoundingBox()); // Expand right AABB to fit objects
	}

	left_node->setAABB(left_bbox);                     // Assign AABBs to new nodes
	right_node->setAABB(right_bbox);


	nodes.push_back(left_node);
	nodes.push_back(right_node);

	build_recursive(left_index, best_split, left_node, depth + 1);
	build_recursive(best_split, right_index, right_node, depth + 1);
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
				if (rec.isHit && rec.t <= length + EPSILON)
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