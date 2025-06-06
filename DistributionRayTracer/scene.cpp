#include <iostream>
#include <string>
#include <fstream>

#include "maths.h"
#include "scene.h"
#include "macros.h"


Triangle::Triangle(Vector& P0, Vector& P1, Vector& P2)
{
	points[0] = P0; points[1] = P1; points[2] = P2;

	// Calculate the Min and Max for bounding box
	Min = Vector(+FLT_MAX, +FLT_MAX, +FLT_MAX);
	Max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for (auto &point : points) {
		if (point.x < Min.x)
			Min.x = point.x;
		if (point.x > Max.x)
			Max.x = point.x;
		if (point.y < Min.y)
			Min.y = point.y;
		if (point.y > Max.y)
			Max.y = point.y;
		if (point.z < Min.z)
			Min.z = point.z;
		if (point.z > Max.z)
			Max.z = point.z;
	}
	// enlarge the bounding box a bit just in case...
	Min -= EPSILON;
	Max += EPSILON;
}

AABB Triangle::GetBoundingBox() {
	return(AABB(Min, Max));
}



// Ray/Triangle intersection test using Tomas Moller-Ben Trumbore algorithm.
HitRecord Triangle::hit(Ray& r) {
	HitRecord rec;
	rec.t = FLT_MAX;      // Default to no intersection
	rec.isHit = false;    // Will be set true if an intersection is found

	// === 1: Compute triangle normal ===
	Vector normal = (points[1] - points[0]) % (points[2] - points[1]); // cross product of two edges
	normal.normalize(); // make sure it's a unit vector for later shading

	// Triangle vertices
	Vector v0 = points[0];
	Vector v1 = points[1];
	Vector v2 = points[2];

	// === 2: M�ller�Trumbore intersection ===
	Vector edge1 = v1 - v0;  // Edge from v0 to v1
	Vector edge2 = v2 - v0;  // Edge from v0 to v2

	Vector h = r.direction % edge2;  // Cross product between ray direction and edge2
	float a = edge1 * h;             // Dot product: if close to 0, ray is parallel to triangle

	float f = 1.0f / a;
	Vector s = r.origin - v0;           // Vector from v0 to ray origin
	float u = f * (s * h);              // Barycentric coordinate u

	if (u < 0.0 || u > 1.0) return rec; // u out of bounds -> outside triangle

	Vector q = s % edge1;              // Cross product of s and edge1
	float v = f * (r.direction * q);   // Barycentric coordinate v

	if (v < 0.0 || u + v > 1.0) return rec; // v or (u+v) out of bounds -> outside triangle

	// === 3. Compute intersection distance t ===
	float t = f * (edge2 * q);         // Distance from ray origin to intersection point

	// If t is positive and large enough, it's a valid hit
	if (t > EPSILON) {
		rec.t = t;
		rec.isHit = true;

		// Compute and store the surface normal at the intersection point
		// (already normalized earlier)
		rec.normal = (edge1 % edge2).normalize();

		return rec;  // Return successful hit
	}

	return rec;  // If intersection was behind the ray origin, discard
}



Plane::Plane(Vector& a_PN, float a_D)
	: PN(a_PN), D(a_D)
{}

Plane::Plane(Vector& P0, Vector& P1, Vector& P2)
{
   float l;
  
   PN = (P1 - P0) % (P2 - P0); // counter-clockwise vectorial product.
   if ((l=PN.length()) == 0.0)
   {
     cerr << "DEGENERATED PLANE!\n";
   }
   else
   {
     PN.normalize();
     D  = -(PN * P0);
   }
}


// Ray/Plane intersection test.
HitRecord Plane::hit(Ray& r) {
	HitRecord rec;
	rec.t = FLT_MAX;     // Initialize to maximum distance (no intersection yet)
	rec.isHit = false;   // Initially assume no hit

	float PNxRd = PN * r.direction;  // Dot product between plane normal and ray direction

	// === 1. Check if Ray is Parallel to the Plane ===

	// If dot product is zero (or very close), ray is parallel to plane -> no hit
	if (fabs(PNxRd) < EPSILON)
		return rec;

	// === 2. Compute Intersection Point (Parameter t) ===

	// Plane equation: (PN � P) + D = 0
	// Ray equation: P = O + tD
	// Solve for t: t = -(PN � O + D) / (PN � D)
	float t = -((PN * r.origin) + D) / PNxRd;

	// === 3. Check if Intersection is in Front of the Ray ===

	if (t > 0) {
		rec.t = t;          // Store intersection distance
		rec.normal = PN;    // Plane normal is constant everywhere
		rec.isHit = true;   // Mark intersection as valid
		return rec;
	}

	// === 4. Otherwise, Return No Hit (Behind the Ray) ===
	return rec;
}


HitRecord Sphere::hit(Ray& r) {
	HitRecord rec;
	rec.t = FLT_MAX;        // Initialize with maximum distance (no hit yet)
	rec.isHit = false;      // Mark as no hit initially

	// === 1. Account for Motion Blur (if enabled) ===
	Vector moved_center = center; // Default: sphere center is static
	if (motion_blur_enabled) {
		velocity.y = 1.0f;  // Assume motion is in the Y direction (for this example)
		moved_center = center + velocity * r.time; // Move sphere according to ray time
	}

	// === 2. Ray-Sphere Intersection Setup ===
	Vector oc = r.origin - moved_center;         // Vector from sphere center to ray origin
	float a = r.direction * r.direction;         // a = dot(D, D)
	float b = 2.0f * (oc * r.direction);         // b = 2 * dot(oc, D)
	float c = (oc * oc) - radius * radius;       // c = dot(oc, oc) - R^2

	float discriminant = b * b - 4 * a * c;      // Discriminant of quadratic equation

	// === 3. Check for No Intersection ===
	if (discriminant < 0) return rec; // Ray misses the sphere

	// === 4. Compute Both Roots of the Quadratic Equation ===
	float sqrt_discriminant = sqrt(discriminant);
	float t1 = (-b - sqrt_discriminant) / (2.0f * a); // Closer intersection
	float t2 = (-b + sqrt_discriminant) / (2.0f * a); // Farther intersection

	// === 5. Choose the First Valid Intersection (In Front of Ray) ===
	if (t1 > EPSILON) {
		rec.t = t1;
	}
	else if (t2 > EPSILON) {
		rec.t = t2;
	}
	else {
		return rec; // Both intersections are behind the ray
	}

	// === 6. Record Hit Information ===
	rec.isHit = true;
	// Compute normal at the hit point: (P - Center) normalized
	rec.normal = (r.origin + r.direction * rec.t - moved_center).normalize();

	return rec;
}



AABB Sphere::GetBoundingBox() {
	Vector a_min = center - Vector(radius, radius, radius);
	Vector a_max = center + Vector(radius, radius, radius);

	return(AABB(a_min, a_max));
}

aaBox::aaBox(Vector& minPoint, Vector& maxPoint) //Axis aligned Box: another geometric object
{
	this->min = minPoint;
	this->max = maxPoint;
}

AABB aaBox::GetBoundingBox() {
	return(AABB(min, max));
}

HitRecord aaBox::hit(Ray& ray)
{
	HitRecord rec;
	rec.t = FLT_MAX;
	rec.isHit = false;

	// --- X slab intersection ---
	// Compute t-values for intersection with the planes perpendicular to the X-axis
	float tmin = (min.x - ray.origin.x) / ray.direction.x;
	float tmax = (max.x - ray.origin.x) / ray.direction.x;

	// Ensure tmin <= tmax
	if (tmin > tmax) std::swap(tmin, tmax);

	// --- Y slab intersection ---
	float tymin = (min.y - ray.origin.y) / ray.direction.y;
	float tymax = (max.y - ray.origin.y) / ray.direction.y;
	if (tymin > tymax) std::swap(tymin, tymax);

	// Check for overlap between X and Y intervals.
	// If the ray misses the box on any axis, there is no intersection.
	if ((tmin > tymax) || (tymin > tmax)) return rec;

	// Update intersection interval to reflect overlap with Y slab
	if (tymin > tmin) tmin = tymin;
	if (tymax < tmax) tmax = tymax;

	// --- Z slab intersection ---
	float tzmin = (min.z - ray.origin.z) / ray.direction.z;
	float tzmax = (max.z - ray.origin.z) / ray.direction.z;
	if (tzmin > tzmax) std::swap(tzmin, tzmax);

	// Final overlap test with Z slab
	if ((tmin > tzmax) || (tzmin > tmax)) return rec;

	// Update tmin and tmax with the Z slab intersection
	if (tzmin > tmin) tmin = tzmin;
	if (tzmax < tmax) tmax = tzmax;

	// If the closest intersection is ahead of the ray origin, it's valid
	if (tmin > EPSILON) {
		rec.t = tmin;
		rec.isHit = true;

		// Compute intersection point
		Vector hitPoint = ray.origin + ray.direction * tmin;

		// Determine which face was hit by comparing hit point with box boundaries
		Vector normal(0, 0, 0);
		if (fabs(hitPoint.x - min.x) < EPSILON) normal = Vector(-1, 0, 0); // Left face
		else if (fabs(hitPoint.x - max.x) < EPSILON) normal = Vector(1, 0, 0); // Right face
		else if (fabs(hitPoint.y - min.y) < EPSILON) normal = Vector(0, -1, 0); // Bottom face
		else if (fabs(hitPoint.y - max.y) < EPSILON) normal = Vector(0, 1, 0); // Top face
		else if (fabs(hitPoint.z - min.z) < EPSILON) normal = Vector(0, 0, -1); // Back face
		else if (fabs(hitPoint.z - max.z) < EPSILON) normal = Vector(0, 0, 1); // Front face

		rec.normal = normal;  // Store surface normal for shading
	}

	return rec; // Return hit info; unchanged if no valid intersection
}



Scene::Scene()
{}

Scene::~Scene()
{
	objects.erase(objects.begin()+1, objects.end()-1);
}

int Scene::getNumObjects()
{
	return objects.size();
}


void Scene::addObject(Object* o)
{
	objects.push_back(o);
}


Object* Scene::getObject(unsigned int index)
{
	if (index >= 0 && index < objects.size())
		return objects[index];
	return NULL;
}


int Scene::getNumLights()
{
	return lights.size();
}


void Scene::addLight(Light* l)
{
	lights.push_back(l);
}


Light* Scene::getLight(unsigned int index)
{
	if (index >= 0 && index < lights.size())
		return lights[index];
	return NULL;
}

void Scene::LoadSkybox(const char *sky_dir)
{
	char *filenames[6];
	char buffer[100];
	const char *maps[] = { "/right.jpg", "/left.jpg", "/top.jpg", "/bottom.jpg", "/front.jpg", "/back.jpg" };

	for (int i = 0; i < 6; i++) {
		strcpy_s(buffer, sizeof(buffer), sky_dir);
		strcat_s(buffer, sizeof(buffer), maps[i]);
		filenames[i] = (char *)malloc(sizeof(buffer));
		strcpy_s(filenames[i], sizeof(buffer), buffer);
	}
	
	
	ILuint ImageName;

	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);

	for (int i = 0; i < 6; i++) {
		ilGenImages(1, &ImageName);
		ilBindImage(ImageName);

		if (ilLoadImage(filenames[i]))  //Image loaded with lower left origin
			printf("Skybox face %d: Image sucessfully loaded.\n", i);
		else
			exit(0);

		ILint bpp = ilGetInteger(IL_IMAGE_BITS_PER_PIXEL);

		ILenum format = IL_RGB;
		printf("bpp=%d\n", bpp);
		if (bpp == 24)
			format = IL_RGB;
		else if (bpp == 32)
			format = IL_RGBA;

		ilConvertImage(format, IL_UNSIGNED_BYTE);

		int size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
		skybox_img[i].img = (ILubyte *)malloc(size);
		ILubyte *bytes = ilGetData();
		memcpy(skybox_img[i].img, bytes, size);
		skybox_img[i].resX = ilGetInteger(IL_IMAGE_WIDTH);
		skybox_img[i].resY = ilGetInteger(IL_IMAGE_HEIGHT);
		format == IL_RGB ? skybox_img[i].BPP = 3 : skybox_img[i].BPP = 4;
		ilDeleteImages(1, &ImageName);
	}
	ilDisable(IL_ORIGIN_SET);
}

Color Scene::GetSkyboxColor(Ray& r) {
	float t_intersec;
	Vector cubemap_coords; //To index the skybox

	float ma;
	CubeMap img_side;
	float sc, tc, s, t;
	unsigned int xp, yp, width, height, bytesperpixel;

	//skybox indexed by the ray direction
	cubemap_coords = r.direction;


	if (fabs(cubemap_coords.x) > fabs(cubemap_coords.y)) {
		ma = fabs(cubemap_coords.x);
		cubemap_coords.x >= 0 ? img_side = LEFT : img_side = RIGHT;    //left cubemap at X = +1 and right at X = -1
	}
	else {
		ma = fabs(cubemap_coords.y);
		cubemap_coords.y >= 0 ? img_side = TOP : img_side = BOTTOM; //top cubemap at Y = +1 and bottom at Y = -1
	}

	if (fabs(cubemap_coords.z) > ma) {
		ma = fabs(cubemap_coords.z);
		cubemap_coords.z >= 0 ? img_side = FRONT : img_side = BACK;   //front cubemap at Z = +1 and back at Z = -1
	}

	switch (img_side) {

	case 0:  //right
		sc = -cubemap_coords.z;
		tc = cubemap_coords.y;
		break;

	case 1:  //left
		sc = cubemap_coords.z;
		tc = cubemap_coords.y;
		break;

	case 2:  //top
		sc = -cubemap_coords.x;
		tc = -cubemap_coords.z;
		break;

	case 3: //bottom
		sc = -cubemap_coords.x;
		tc = cubemap_coords.z;
		break;

	case 4:  //front
		sc = -cubemap_coords.x;
		tc = cubemap_coords.y;
		break;

	case 5: //back
		sc = cubemap_coords.x;
		tc = cubemap_coords.y;
		break;
	}

	double invMa = 1 / ma;
	s = (sc * invMa + 1) / 2;
	t = (tc * invMa + 1) / 2;

	width = skybox_img[img_side].resX;
	height = skybox_img[img_side].resY;
	bytesperpixel = skybox_img[img_side].BPP;

	xp = int((width - 1) * s);
	xp < 0 ? 0 : (xp > (width - 1) ? width - 1 : xp);
	yp = int((height - 1) * t);
	yp < 0 ? 0 : (yp > (height - 1) ? height - 1 : yp);

	float red = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel]);
	float green = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel + 1]);
	float blue = u8tofloat(skybox_img[img_side].img[(yp*width + xp) * bytesperpixel + 2]);

	return(Color(red, green, blue));
}




////////////////////////////////////////////////////////////////////////////////
// P3F file parsing methods.
//
void next_token(ifstream& file, char *token, const char *name)
{
  file >> token;
  if (strcmp(token, name))
    cerr << "'" << name << "' expected.\n";
}


bool Scene::load_p3f(const char *name)
{
  const	int	lineSize = 1024;
  string	cmd;
  char		token	[256];
  ifstream	file(name, ios::in);
  Material *	material;

  material = NULL;
  this->SetSkyBoxFlg(false);  //init with no skybox

  if (file >> cmd)
  {
    while (true)
    {
      if (cmd == "accel") {  //Acceleration data structure
		string accel_type; // type of acceleration data structure
		file >> accel_type;
		if (accel_type == "none")
			this->SetAccelStruct(NONE);
		else if (accel_type == "grid")
			this->SetAccelStruct(GRID_ACC);
		else if (accel_type == "bvh")
			this->SetAccelStruct(BVH_ACC);
		else {
			printf("Unsupported acceleration type\n");
			break;
		}

	  }

	  else if (cmd == "spp")    //samples per pixel
	  {
		  unsigned int spp; // number of samples per pixel 

		  file >> spp;
		  this->SetSamplesPerPixel(spp);
	  }
	  else if (cmd == "mat")   //Material
      {
		  double Kd, Ks, Shine, T, ior;
		  Color cd, cs;

		  file >> cd >> Kd >> cs >> Ks >> Shine >> T >> ior;

		  material = new Material(cd, Kd, cs, Ks, Shine, T, ior);
      }

      else if (cmd == "s")    //Sphere
      {
	     Vector center;
    	 float radius;
         Sphere* sphere;

	    file >> center >> radius;
        sphere = new Sphere(center,radius);
	    if (material) sphere->SetMaterial(material);
        this->addObject( (Object*) sphere);
      }

	  else if (cmd == "box")    //axis aligned box
	  {
		  Vector minpoint, maxpoint;
		  aaBox	*box;

		  file >> minpoint >> maxpoint;
		  box = new aaBox(minpoint, maxpoint);
		  if (material) box->SetMaterial(material);
		  this->addObject((Object*)box);
	  }
	  else if (cmd == "p")  // Polygon: just accepts triangles for now
      {
		  Vector P0, P1, P2;
		  Triangle* triangle; 
		  unsigned total_vertices;
		  
		  file >> total_vertices;
		  if (total_vertices == 3)
		  {
			  file >> P0 >> P1 >> P2;
			  triangle = new Triangle(P0, P1, P2);
			  if (material) triangle->SetMaterial(material);
			  this->addObject( (Object*) triangle);
		  }
		  else
		  {
			  cerr << "Unsupported number of vertices.\n";
			  break;
		  }
      }

	  else if (cmd == "mesh") {
		  unsigned total_vertices, total_faces;
		  unsigned P0, P1, P2;
		  Triangle* triangle;
		  Vector* verticesArray, vertex;

		  file >> total_vertices >> total_faces;
		  verticesArray = (Vector*)malloc(total_vertices * sizeof(Vector));
		  for (int i = 0; i < total_vertices; i++) {
			  file >> vertex;
			  verticesArray[i] = vertex;
		  }

		  for (int i = 0; i < total_faces; i++) {
			  file >> P0 >> P1 >> P2;
			  if (P0 > 0) {
				  P0 -= 1;
				  P1 -= 1;
				  P2 -= 1;
			  }
			  else {
				  P0 += total_vertices;
				  P1 += total_vertices;
				  P2 += total_vertices;
			  }
			  triangle = new Triangle(verticesArray[P0], verticesArray[P1], verticesArray[P2]); //vertex index start at 1
			  if (material) triangle->SetMaterial(material);
			  this->addObject((Object*)triangle);
		  }
	  }
      
      else if (cmd == "npl")  //Plane in Hessian form
	  {
          Vector N;
          float D;
          Plane* plane;

          file >> N >> D;
          plane = new Plane(N, D);
	      if (material) plane->SetMaterial(material);
          this->addObject( (Object*) plane);
	  }
	  else if (cmd == "pl")  // General Plane
	  {
          Vector P0, P1, P2;
		  Plane* plane;

          file >> P0 >> P1 >> P2;
          plane = new Plane(P0, P1, P2);
	      if (material) plane->SetMaterial(material);
          this->addObject( (Object*) plane);
	  }

      else if (cmd == "light")  // Need to check light color since by default is white
      {
	    Vector pos, at;
        Color color;
		Vector v1, v2;
		unsigned int grid_res;

		string type;
		file >> type;

		if (type == "punctual") {
			file >> pos >> color;
			this->addLight(new Light(pos, color));
		}
		else if(type == "quad") {
			file >> pos >> color >> v1 >> v2 >> grid_res;
			this->addLight(new Light(pos, color, v1, v2, grid_res));
		}
		else
		{
			cerr << "Unsupported light type.\n";
			break;
		}
	    
      }
      else if (cmd == "camera")
      {
	    Vector up, from, at;
	    float fov, hither;
	    int xres, yres;
        Camera* camera;
		float focal_ratio; //ratio beteween the focal distance and the viewplane distance
		float aperture_ratio; // number of times to be multiplied by the size of a pixel

	    next_token (file, token, "eye");
	    file >> from;

	    next_token (file, token, "at");
	    file >> at;

	    next_token (file, token, "up");
	    file >> up;

	    next_token (file, token, "angle");
	    file >> fov;

	    next_token (file, token, "hither");
	    file >> hither;

	    next_token (file, token, "resolution");
	    file >> xres >> yres;

		next_token(file, token, "aperture");
		file >> aperture_ratio;

		next_token(file, token, "focal");
		file >> focal_ratio;
	    // Create Camera
		camera = new Camera( from, at, up, fov, hither, 1000.0*hither, xres, yres, aperture_ratio, focal_ratio);
        this->SetCamera(camera);
      }

	  else if (cmd == "bclr")   //Background color
	  {
	  Color bgcolor;
	  file >> bgcolor;
	  this->SetBackgroundColor(bgcolor);
	  }

	  else if (cmd == "env")
	  {
	  file >> token;

	  this->LoadSkybox(token);
	  this->SetSkyBoxFlg(true);
	  }

	  /*
      else if (cmd == "b")   //Background color
      {
		Color bgcolor;
		bool skybox_flg;

		file >> token;
		if (!strcmp(token, "clr")) {
			skybox_flg = false;
			this->SetSkyBoxFlg(skybox_flg);
			file >> bgcolor;
			this->SetBackgroundColor(bgcolor);
		}
		else if (!strcmp(token, "cubemap")) {
			skybox_flg = true;
			this->SetSkyBoxFlg(skybox_flg);
			file >> token;   
			this->LoadSkybox(token);
		}
		else {
			skybox_flg = false;
			this->SetSkyBoxFlg(skybox_flg);
			file.ignore(lineSize, '\n');
			this->SetBackgroundColor(Color(0, 0, 0));
		}

      }
	  */

      else if (cmd[0] == '#')
      {
	    file.ignore (lineSize, '\n');
      }
      else
      {
	    cerr << "unknown command '" << cmd << "'.\n";
	    break;
      }
      if (!(file >> cmd))
        break;
    }
  }

  file.close();
  return true;
};

void Scene::create_random_scene() {
	Camera* camera;
	Material* material;
	Sphere* sphere;

	set_rand_seed(time(NULL)* time(NULL)* time(NULL));
	material = NULL;
	this->SetSkyBoxFlg(false);  //init with no skybox

	this->SetBackgroundColor(Color(0.5, 0.7, 1.0));
	this->SetAccelStruct(NONE);
	//this->SetAccelStruct(GRID_ACC);
	this->SetSamplesPerPixel(0);
	camera = new Camera(Vector(-5.312192, 4.456562, 11.963158), Vector(0.0, 0.0, 0), Vector(0.0, 1.0, 0.0), 40.0, 0.01, 10000.0, 800, 600, 0, 1.5f);
	this->SetCamera(camera);

	this->addLight(new Light(Vector(7, 10, -5), Color(1.0, 1.0, 1.0)));
	this->addLight(new Light(Vector(-7, 10, -5), Color(1.0, 1.0, 1.0)));
	this->addLight(new Light(Vector(0, 10, 7), Color(1.0, 1.0, 1.0)));

	material = new Material(Color(0.5, 0.5, 0.5), 1.0, Color(0.0, 0.0, 0.0), 0.0, 10, 0, 1);


	sphere = new Sphere(Vector(0.0, -1000, 0.0), 1000.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	for(int a = -5; a < 5; a++)
		for (int b = -5; b < 5; b++) {

			double choose_mat = rand_double();

			Vector center = Vector(a + 0.9 * rand_double(), 0.2, b + 0.9 * rand_double());

			if ((center - Vector(4.0, 0.2, 0.0)).length() > 0.9) {
				if (choose_mat < 0.4) {  //diffuse
					material = new Material(Color(rand_double(), rand_double(), rand_double()), 1.0, Color(0.0, 0.0, 0.0), 0.0, 10, 0, 1);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}
				else if (choose_mat < 0.7) {   //metal
					material = new Material(Color(0.0, 0.0, 0.0), 0.0, Color(rand_double(0.5, 1), rand_double(0.5, 1), rand_double(0.5, 1)), 1.0, 220, 0, 1);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}
				else {   //glass 
					//material = new Material(Color(0.8, 0.3, 0.3), 0.0, Color(1.0, 1.0, 1.0), 0.7, 20, 1, 1.5);
					material = new Material(Color(rand_double(0.6, 1), rand_double(0.6, 1), rand_double(0.6, 1)), 0.0, Color(1.0, 1.0, 1.0), 0.7, 20, 1, 1.5);
					sphere = new Sphere(center, 0.2);
					if (material) sphere->SetMaterial(material);
					this->addObject((Object*)sphere);
				}

			}

		}

	material = new Material(Color(1.0, 1.0, 1.0), 0.0, Color(1.0, 1.0, 1.0), 0.7, 20, 1, 1.5);
	sphere = new Sphere(Vector(0.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	material = new Material(Color(0.4, 0.2, 0.1), 0.9, Color(1.0, 1.0, 1.0), 0.0, 10, 0, 1.0);
	sphere = new Sphere(Vector(-4.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);

	material = new Material(Color(0.4, 0.2, 0.1), 0.0, Color(0.7, 0.6, 0.5), 1.0, 220, 0, 1.0);
	sphere = new Sphere(Vector(4.0, 1.0, 0.0), 1.0);
	if (material) sphere->SetMaterial(material);
	this->addObject((Object*)sphere);
}