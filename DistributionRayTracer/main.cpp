 ///////////////////////////////////////////////////////////////////////
//
// P3D Course
// (c) 2025 by Jo�o Madeiras Pereira
//Distribution Ray Tracing P3F scenes and drawing points with Modern OpenGL
// It explores parallelism through OMP
//
///////////////////////////////////////////////////////////////////////
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <chrono>
#include <conio.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <IL/il.h>

#include "scene.h"
#include "rayAccelerator.h"
#include "maths.h"
#include "macros.h"

//Enable OpenGL drawing.  
bool drawModeEnabled = true;
bool P3F_scene = true; //choose between P3F scene or the built-in Peter Shirley scene
bool Progressive_flg = false;
bool motion_blur_enabled = true;

#define MAX_DEPTH 4  //number of bounces

#define CAPTION "Accel Distribution RT"
#define VERTEX_COORD_ATTRIB 0
#define COLOR_ATTRIB 1
#define MAX_SAMPLES 10000

// Frame counting and FPS computation
unsigned int FramesPerSecond = 0, FPS = 0;
unsigned long FrameCount = 1;

// Current Camera Position
float camX, camY, camZ;

//Original Camera position;
Vector Eye;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 0.0f;
float _beta = 0.0f;
float r = 4.0f;

// Color Gamma correction
double invGamma = 1.0f / 2.2f;

// Points defined by 2 attributes: positions which are stored in vertices array and colors which are stored in colors array
float *colors;
float *vertices;
int size_vertices;
int size_colors;

//Array of Pixels to be stored in a file by using DevIL library
uint8_t *img_Data;

GLfloat m[16];  //projection matrix initialized by ortho function

GLuint VaoId;
GLuint VboId[2];

GLuint VertexShaderId, FragmentShaderId, ProgramId;
GLint UniformId;

Scene* scene = NULL;
Grid* grid_ptr = NULL;
BVH* bvh_ptr = NULL;

int RES_X, RES_Y;


int WindowHandle = 0;

bool AA = false;
unsigned int spp; //samples per pixel
bool DOF = false;
bool SoftShadows = false;

accelerator Accel_Struct = NONE;

/////////////////////////////////////////////////////////////////////// ERRORS

bool isOpenGLError() {
	bool isError = false;
	GLenum errCode;
	const GLubyte *errString;
	while ((errCode = glGetError()) != GL_NO_ERROR) {
		isError = true;
		errString = gluErrorString(errCode);
		std::cerr << "OpenGL ERROR [" << errString << "]." << std::endl;
	}
	return isError;
}

void checkOpenGLError(std::string error)
{
	if(isOpenGLError()) {
		std::cerr << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

/////////////////////////////////////////////////////////////////////// SHADERs

const GLchar* VertexShader =
{
	"#version 430 core\n"

	"in vec2 in_Position;\n"
	"in vec3 in_Color;\n"
	"uniform mat4 Matrix;\n"
	"out vec4 color;\n"

	"void main(void)\n"
	"{\n"
	"	vec4 position = vec4(in_Position, 0.0, 1.0);\n"
	"	color = vec4(in_Color, 1.0);\n"
	"	gl_Position = Matrix * position;\n"

	"}\n"
};

const GLchar* FragmentShader =
{
	"#version 430 core\n"

	"in vec4 color;\n"
	"out vec4 out_Color;\n"

	"void main(void)\n"
	"{\n"
	"	out_Color = color;\n"
	"}\n"
};

void createShaderProgram()
{
	VertexShaderId = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(VertexShaderId, 1, &VertexShader, 0);
	glCompileShader(VertexShaderId);

	FragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(FragmentShaderId, 1, &FragmentShader, 0);
	glCompileShader(FragmentShaderId);

	ProgramId = glCreateProgram();
	glAttachShader(ProgramId, VertexShaderId);
	glAttachShader(ProgramId, FragmentShaderId);

	glBindAttribLocation(ProgramId, VERTEX_COORD_ATTRIB, "in_Position");
	glBindAttribLocation(ProgramId, COLOR_ATTRIB, "in_Color");
	
	glLinkProgram(ProgramId);
	UniformId = glGetUniformLocation(ProgramId, "Matrix");

	checkOpenGLError("ERROR: Could not create shaders.");
}

void destroyShaderProgram()
{
	glUseProgram(0);
	glDetachShader(ProgramId, VertexShaderId);
	glDetachShader(ProgramId, FragmentShaderId);

	glDeleteShader(FragmentShaderId);
	glDeleteShader(VertexShaderId);
	glDeleteProgram(ProgramId);

	checkOpenGLError("ERROR: Could not destroy shaders.");
}

/////////////////////////////////////////////////////////////////////// VAOs & VBOs


void createBufferObjects()
{
	glGenVertexArrays(1, &VaoId);
	glBindVertexArray(VaoId);
	glGenBuffers(2, VboId);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);

	/* S� se faz a aloca��o dos arrays glBufferData (NULL), e o envio dos pontos para a placa gr�fica
	� feito na drawPoints com GlBufferSubData em tempo de execu��o pois os arrays s�o GL_DYNAMIC_DRAW */
	glBufferData(GL_ARRAY_BUFFER, size_vertices, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glVertexAttribPointer(VERTEX_COORD_ATTRIB, 2, GL_FLOAT, 0, 0, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferData(GL_ARRAY_BUFFER, size_colors, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(COLOR_ATTRIB);
	glVertexAttribPointer(COLOR_ATTRIB, 3, GL_FLOAT, 0, 0, 0);
	
// unbind the VAO
	glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
//	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB); 
//	glDisableVertexAttribArray(COLOR_ATTRIB);
	checkOpenGLError("ERROR: Could not create VAOs and VBOs.");
}

void destroyBufferObjects()
{
	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glDisableVertexAttribArray(COLOR_ATTRIB);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glDeleteBuffers(1, VboId);
	glDeleteVertexArrays(1, &VaoId);
	checkOpenGLError("ERROR: Could not destroy VAOs and VBOs.");
}

void drawPoints()
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	glBindVertexArray(VaoId);
	glUseProgram(ProgramId);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_vertices, vertices);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_colors, colors);

	glUniformMatrix4fv(UniformId, 1, GL_FALSE, m);
	glDrawArrays(GL_POINTS, 0, RES_X*RES_Y);
	//glFinish();
	glUseProgram(0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	checkOpenGLError("ERROR: Could not draw scene.");
	
}

ILuint saveImgFile(const char *filename) {
	ILuint ImageId;

	ilEnable(IL_FILE_OVERWRITE);
	ilGenImages(1, &ImageId);
	ilBindImage(ImageId);

	ilTexImage(RES_X, RES_Y, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, img_Data /*Texture*/);
	ilSaveImage(filename);

	ilDisable(IL_FILE_OVERWRITE);
	ilDeleteImages(1, &ImageId);
	if (ilGetError() != IL_NO_ERROR)return ilGetError();

	return IL_NO_ERROR;
}

/////////////////////////////////////////////////////////////////////// CALLBACKS

void timer(int value)
{
	FramesPerSecond = FPS;
	if(!Progressive_flg) {
		std::ostringstream oss;
		oss << CAPTION << ": " << FramesPerSecond << " FPS @ (" << RES_X << "x" << RES_Y << ")";
		std::string s = oss.str();
		glutSetWindow(WindowHandle);
		glutSetWindowTitle(s.c_str());
	}
	FPS = 0;
	glutTimerFunc(1000, timer, 0);
}


///////////////////////////////////////////////////YOUR CODE HERE////////////////////////////////////////////////////////////////////////

float schlick(float cosTheta, float ior1, float ior2) {
	float r0 = (ior1 - ior2) / (ior1 + ior2);
	r0 = pow(r0, 2);
	return r0 + (1.0f - r0) * pow(1.0f - cosTheta, 5.0f);
}


Color rayTracing(Ray ray, int depth, float ior_1, Vector lightSample)  //index of refraction of medium 1 where the ray is travelling
{
	Color color_Acc; //Class constructor init the color with zero

	Object* hitObj = NULL; //nearest object
	HitRecord closestHit;  //isHit=false and t=FLT_MAX
	Vector hitPoint; //closest hit point
	Vector N;
	HitRecord auxRec;
	bool skybox_flg;
	int num_lights = scene->getNumLights();
	skybox_flg = scene->GetSkyBoxFlg();
	Accel_Struct = scene->GetAccelStruct();   //Type of acceleration data structure
	int num_objects = scene->getNumObjects();

	if (Accel_Struct == NONE) {  //no acceleration

		// Find closest hit
		Object* obj = NULL;

		for (int i = 0; i < num_objects; i++)
		{
			obj = scene->getObject(i);
			auxRec = obj->hit(ray);


			if (auxRec.isHit && auxRec.t < closestHit.t) {
				closestHit = auxRec;
				hitObj = obj;

			}
		}

		if (hitObj == NULL) {  // No intersected object
			if (skybox_flg)  //skybox cubemap overrides background color 
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = (scene->GetBackgroundColor());

			return color_Acc.clamp();
		}
	}

	else if (Accel_Struct == GRID_ACC) {  // regular Grid
		if (!grid_ptr->Traverse(ray, &hitObj, closestHit)) {
			if (skybox_flg)
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = (scene->GetBackgroundColor());
			return color_Acc.clamp();
		}
	}

	else if (Accel_Struct == BVH_ACC) { //BVH
		if (!bvh_ptr->Traverse(ray, &hitObj, closestHit)) {
			if (skybox_flg)
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = (scene->GetBackgroundColor());
			return color_Acc.clamp();
		}
	}

	hitPoint = ray.origin + ray.direction * closestHit.t;
	N = closestHit.normal.normalize();

	//CALCULATE THE COLOR OF THE PIXEL

	bool outside = (ray.direction * N) < 0.0f;

	if (!outside) N = -N;

	Material* mat = hitObj->GetMaterial();

	Color diff_color = mat->GetDiffColor();
	Color spec_color= mat->GetSpecColor();
	float kr = mat->GetReflection();
	float kd = mat->GetDiffuse();
	float ks = mat->GetSpecular();
	float shine = mat->GetShine();
	float ior2 = mat->GetRefrIndex();
	float trans = mat->GetTransmittance();
	Vector V = -ray.direction.normalize(); // view vector

	float offset = 1e-4f;

	// ========== SOFT SHADOWS ==========
	Vector lightPos;
	for (int j = 0; j < num_lights; j++) {

		Light* light = scene->getLight(j);
		float lightDist;

		if (light->type == QUAD) {
			// Use jittered sampling for quad lights
			lightPos = light->getAreaLightPoint(lightSample); //GET JITTERED LIGHT SAMPLE IN THE QUAD LIGHT PERIMETER!!
		}
		else {  // Point light
			lightPos = light->position;
		}


		Vector L = lightPos - hitPoint;
		Vector Ls = L;
		L.normalize();
		Vector H = (L + V).normalize();

		lightDist = Ls.length();


		float NdotL = std::max(N * L, 0.0f);
		float NdotH = std::max(N * H, 0.0f);

		// === Shadow ray ===

		Vector shadowDir;
		if (Accel_Struct == GRID_ACC) {
			shadowDir = Ls;
		}
		if (Accel_Struct == BVH_ACC) {
			shadowDir = Ls;
		}
		else {
			shadowDir = L;
		}

		Ray shadowRay(hitPoint + N * offset, shadowDir); // offset a bit to avoid acne
		bool inShadow = false;
		
		if (Accel_Struct == GRID_ACC) {
			inShadow = grid_ptr->Traverse(shadowRay);
		}
		else if (Accel_Struct == BVH_ACC) {
			inShadow = bvh_ptr->Traverse(shadowRay);
		}
		else {
			for (int o = 0; o < num_objects; o++) {
				if (scene->getObject(o) == hitObj) continue; // skip self
				HitRecord shadowHit = scene->getObject(o)->hit(shadowRay);
				if (shadowHit.isHit && shadowHit.t > offset && shadowHit.t < L.length()) {
					inShadow = true;
					break;
				}
			}
		}

		

		if (!inShadow) {

			Color diffuseTerm = diff_color * kd * NdotL;
			Color specularTerm = spec_color* ks * pow(NdotH, shine);

			color_Acc += (diffuseTerm + specularTerm);
		}
	}

	if (depth > MAX_DEPTH) {
		return color_Acc;
	}

	float kr_fresnel = kr;
	if (!outside) ior2 = 1.0;
	float eta = ior_1 / ior2;
	Vector Vt = N * (V * N) - V;
	float sin_i = Vt.length();
	Vector t = Vt/Vt.length();
	float sin_t = eta * sin_i;

	if (trans == 1 && sin_t < 1) {

		float sin_t2 = pow(sin_t, 2);
		float cos_t = sqrt(1 - sin_t2);
		Vector r_t = (t * sin_t + (-N) * cos_t).normalize();

		// === Schlick�s Approximation ===
		float cos_i = N * V;
		
		float cosTheta = cos_i;

		if (ior_1 > ior2) {
			cosTheta = cos_t;
		}
		else {
			cosTheta = cos_i;
		}

		float r0 = (ior_1 - ior2) / (ior_1 + ior2);
		r0 = pow(r0, 2);
		kr_fresnel = r0 + (1.0f - r0) * pow((1.0f - cosTheta), 5);

		Ray refractRay(hitPoint - (N * offset), r_t);
		Color refractColor = rayTracing(refractRay, depth + 1, ior2, lightPos).clamp();

		Color one = Color(1.0f, 1.0f, 1.0f);
		if (!outside) {
			refractColor = refractColor * ((one - diff_color) * (-closestHit.t)).exp_();
			
		}

		color_Acc += refractColor * (1 -kr_fresnel);

	}
	else if (trans > 0.0f && sin_t >= 1) {
		kr_fresnel = 1;
	}

	if (ks > 0) {
		// === Reflection ===
		Vector reflectDir = N * (V * N) * 2.0f - V;

		float roughness_param = 0.5;

		reflectDir = (reflectDir + rnd_unit_sphere() * roughness_param).normalize();
		Ray reflectRay(hitPoint + N * offset, reflectDir);

		Color reflectColor = rayTracing(reflectRay, depth + 1, ior_1, lightPos).clamp();
		float k_ref = kr_fresnel;

		if (reflectDir*N > 0)
			color_Acc += reflectColor * k_ref * spec_color;

	}
	
	return color_Acc.clamp();
}


// Render function by primary ray casting from the eye towards the scene's objects
void renderScene()	
{
	unsigned int counter = 0;
	set_rand_seed(time(NULL) * time(NULL)); // Use current time as seed for random generator 

	if (drawModeEnabled) {
		//glClear(GL_COLOR_BUFFER_BIT);
		scene->GetCamera()->SetEye(Vector(camX, camY, camZ));  //Camera motion
	}
	

	if (Progressive_flg){          ///////////////////////// ZONE A  - Progressive RayTracer/////////////////////////////////////////
		if (FrameCount < MAX_SAMPLES) {
#pragma omp parallel for collapse(2)
			for (int y = 0; y < RES_Y; y++) {
				for (int x = 0; x < RES_X; x++) {

					int index_pos = 0;
					int index_col = 0;
					Color color;
					Ray ray;
					Vector pixel_sample;  //viewport coordinates
					Vector light_sample = Vector(0.0f, 0.0f, 0.0f); // sample in Light coordinates
					float jitter_time =0;
					if (motion_blur_enabled == true) {
						jitter_time = (float)rand() / (float)RAND_MAX; //for motion blur
					}
					

					pixel_sample.x = x + rand_double();
					pixel_sample.y = y + rand_double();

					if (!DOF) ray = scene->GetCamera()->PrimaryRay(pixel_sample, jitter_time);
					else {        // sample_unit_disk() returns [-1 1] and aperture is the diameter of the lens

						Vector lens_sample = rnd_unit_disk() * scene->GetCamera()->GetAperture() / 2.0f;  // lens sample in Camera coordinates

						/////////PROGRAM THE FOLLOWING FUNCTION//////////////////////
						ray = scene->GetCamera()->PrimaryRay(lens_sample, pixel_sample, jitter_time);
					}
					/////////PROGRAM THE FOLLOWING FUNCTION//////////////////////
					color = rayTracing(ray, 1, 1.0, Vector(rand_float(), rand_float(), 0.0f));

					index_pos = 2 * (x + RES_X * y);
					vertices[index_pos] = (float)x;
					vertices[index_pos + 1] = (float)y;

					index_col = 3 * (x + RES_X * y);
					///////////UNDERSTAND AND EXPLAIN THE FOLLOWING CODE//////////////////////
					if (FrameCount == 1) {
						colors[index_col] = (float)color.r();
						colors[index_col + 1] = (float)color.g();
						colors[index_col + 2] = (float)color.b();
					}

					else {
						colors[index_col] = lerp(colors[index_col], (float)color.r(), 1.0 / FrameCount);
						index_col++;
						colors[index_col] = lerp(colors[index_col], (float)color.g(), 1.0 / FrameCount);
						index_col++;
						colors[index_col] = lerp(colors[index_col], (float)color.b(), 1.0 / FrameCount);
					}
				}
			}
		}
		drawPoints();
		if(FrameCount != MAX_SAMPLES)  FrameCount++;
		FPS++;  
		std::ostringstream oss;
		oss << CAPTION << ": " << FramesPerSecond << " FPS @ (" << RES_X << "x" << RES_Y << ") @ Samples number: " << FrameCount;
		std::string s = oss.str();
		glutSetWindow(WindowHandle);
		glutSetWindowTitle(s.c_str());
		glutSwapBuffers();
	}

	//////////////ZONE B - NOT Progressive RayTracer////////////////////////////////////
	else {				
#pragma omp parallel for collapse(2)
		for (int y = 0; y < RES_Y; y++) {
			for (int x = 0; x < RES_X; x++) {
				Color color;
				int index_pos = 0;
				int index_col = 0;
				Ray ray;
				Vector pixel_sample;  //viewport coordinates
				Vector light_sample = Vector(0.0f, 0.0f, 0.0f); // sample in Light coordinates
				float jitter_time = 0;
				if (motion_blur_enabled == true) {
					jitter_time = (float)rand() / (float)RAND_MAX; //for motion blur
				}

				////// ZONE B.1  -  Distribution Ray Tracer: pixel, area light and lens supersampling with jittering (or stratified)
				if(AA) {  
					int n = (int)sqrt(spp);

					// Create arrays for jittered samples
					Vector* r = new Vector[spp]; // pixel samples
					Vector* s = new Vector[spp]; // light samples

					// Generate jittered samples for pixels and lights
					for (int p = 0; p < spp; p++) {
						index_col = p / n;
						index_pos = p % n;

						// Jittered samples for pixel
						float epsilon_x = (float)rand() / (float)RAND_MAX;
						float epsilon_y = (float)rand() / (float)RAND_MAX;
						r[p] = Vector((index_pos + epsilon_x) / n, (index_col + epsilon_y) / n, 0.0f);

						// Jittered samples for light
						epsilon_x = (float)rand() / (float)RAND_MAX;
						epsilon_y = (float)rand() / (float)RAND_MAX;
						s[p] = Vector(epsilon_x, epsilon_y, 0.0f);
					}

					// Shuffle the light samples
					for (int i = spp - 1; i > 0; i--) {
						int j = rand() % (i + 1); // random integer between 0 and i
						// Swap elements i and j
						Vector temp = s[i];
						s[i] = s[j];
						s[j] = temp;
					}

#pragma omp parallel for
					for (int p = 0; p < spp; p++) {
						pixel_sample = Vector(x + r[p].x, y + r[p].y, 0.0f);

						if (!DOF) ray = scene->GetCamera()->PrimaryRay(pixel_sample, jitter_time);
						else {        // sample_unit_disk() returns [-1 1] and aperture is the diameter of the lens

							Vector lens_sample = rnd_unit_disk() * scene->GetCamera()->GetAperture() / 2.0f;  // lens sample in Camera coordinates

							/////////PROGRAM THE FOLLOWING FUNCTION//////////////////////
							ray = scene->GetCamera()->PrimaryRay(lens_sample, pixel_sample, jitter_time);
						}

						/////////PROGRAM THE FOLLOWING FUNCTION//////////////////////
						color += rayTracing(ray, 1, 1.0, s[p]);
					}
					color *= 1.0 / ((float)spp);

					// Clean up
					delete[] r;
					delete[] s;
				}

				//ZONE B.2  - Whitted ray tracer  (without antialiasing)
				else {	

					pixel_sample.x = x + 0.5f;
					pixel_sample.y = y + 0.5f;
					Ray ray1 = scene->GetCamera()->PrimaryRay(pixel_sample, jitter_time);

					// Regular sampling for soft shadows when no AA
					if (SoftShadows) {
						Light* light = scene->getLight(0);
						const int lightSamples = light->gridRes; // 4x4 grid usually
						Color tempColor;

						for (int s = 0; s < lightSamples; s++) {
							// Calculate regular grid coordinates
							int gridSize = (int)sqrt(lightSamples);
							float u = (s % gridSize + 0.5f) / gridSize;
							float v = (s / gridSize + 0.5f) / gridSize;
							Vector regularLightSample = Vector(u, v, 0.0f);

							tempColor += rayTracing(ray1, 1, 1.0, regularLightSample);
						}
						color = tempColor * (1.0f / lightSamples); // Average the samples
					}
					else {
						// Default center sample for hard shadows
						color = rayTracing(ray1, 1, 1.0, Vector(0.5f, 0.5f, 0.0f));
					}
				}

				if (drawModeEnabled) {
					index_pos = 2 * (x + RES_X * y);
					vertices[index_pos] = (float)x;
					vertices[index_pos + 1] = (float)y;

					index_col = 3 * (x + RES_X * y);
					colors[index_col] = (float)color.r();
					colors[index_col+1] = (float)color.g();
					colors[index_col+2] = (float)color.b();
				}
				else {
					img_Data[counter++] = u8fromfloat((float)color.r());
					img_Data[counter++] = u8fromfloat((float)color.g());
					img_Data[counter++] = u8fromfloat((float)color.b());
				}
			}
		}


		if (drawModeEnabled){
			FPS++;
			drawPoints();
			glutSwapBuffers();
		}
		else {
			printf("Rendering ended!\n"); 	
			if (saveImgFile("RT_Output.png") != IL_NO_ERROR) {
				printf("Error saving Image file\n");
				exit(0);
			}
			printf("Image file created\n");
		}
	}
}

// Callback function for glutCloseFunc
void cleanup()
{
	destroyShaderProgram();
	destroyBufferObjects();
}

void ortho(float left, float right, float bottom, float top, 
			float nearp, float farp)
{
	m[0 * 4 + 0] = 2 / (right - left);
	m[0 * 4 + 1] = 0.0;
	m[0 * 4 + 2] = 0.0;
	m[0 * 4 + 3] = 0.0;
	m[1 * 4 + 0] = 0.0;
	m[1 * 4 + 1] = 2 / (top - bottom);
	m[1 * 4 + 2] = 0.0;
	m[1 * 4 + 3] = 0.0;
	m[2 * 4 + 0] = 0.0;
	m[2 * 4 + 1] = 0.0;
	m[2 * 4 + 2] = -2 / (farp - nearp);
	m[2 * 4 + 3] = 0.0;
	m[3 * 4 + 0] = -(right + left) / (right - left);
	m[3 * 4 + 1] = -(top + bottom) / (top - bottom);
	m[3 * 4 + 2] = -(farp + nearp) / (farp - nearp);
	m[3 * 4 + 3] = 1.0;
}

void reshape(int w, int h)
{
    glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, w, h);
	ortho(0, (float)RES_X, 0, (float)RES_Y, -1.0, 1.0);
}

void processKeys(unsigned char key, int xx, int yy)
{
	switch (key) {

	case 27:
		glutLeaveMainLoop();
		break;

	case 'p':
		if (Progressive_flg) Progressive_flg = false;
		else { Progressive_flg = true; FrameCount = 1; }
		break;

	case 'r':
		camX = Eye.x;
		camY = Eye.y;
		camZ = Eye.z;
		r = Eye.length();
		_beta = asinf(camY / r) * 180.0f / 3.14f;
		alpha = atanf(camX / camZ) * 180.0f / 3.14f;
		FrameCount = 1;
		break;

	case 'c':
		printf("Camera Spherical Coordinates (%f, %f, %f)\n", r, _beta, alpha);
		printf("Camera Cartesian Coordinates (%f, %f, %f)\n", camX, camY, camZ);
		break;
	}
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		//FrameCount = 1;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha -= (xx - startX);
			_beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX = -xx + startX;
	deltaY = yy - startY;

	// left mouse button: move camera
	if (tracking == 1) {


		alphaAux = alpha + deltaX;
		betaAux = _beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = _beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	FrameCount = 1;
	camX = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camZ = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camY = rAux * sin(betaAux * 3.14f / 180.0f);

	//  uncomment this if not using an idle func
	//	glutPostRedisplay();
}

void mouseWheel(int wheel, int direction, int x, int y) {

	r += direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;
	FrameCount = 1;
	camX = r * sin(alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camZ = r * cos(alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camY = r * sin(_beta * 3.14f / 180.0f);

	//  uncomment this if not using an idle func
	//	glutPostRedisplay();
}


/////////////////////////////////////////////////////////////////////// SETUP

void setupCallbacks() 
{
	glutKeyboardFunc(processKeys);
	glutCloseFunc(cleanup);
	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc(mouseWheel);
	glutTimerFunc(1000, timer, 0);
	glutIdleFunc(renderScene);
}

void setupGLEW() {
	glewExperimental = GL_TRUE;
	GLenum result = glewInit() ; 
	if (result != GLEW_OK) { 
		std::cerr << "ERROR glewInit: " << glewGetString(result) << std::endl;
		exit(EXIT_FAILURE);
	} 
	GLenum err_code = glGetError();
	printf ("Vendor: %s\n", glGetString (GL_VENDOR));
	printf ("Renderer: %s\n", glGetString (GL_RENDERER));
	printf ("Version: %s\n", glGetString (GL_VERSION));
	printf ("GLSL: %s\n", glGetString (GL_SHADING_LANGUAGE_VERSION));
}

void setupGLUT(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitContextVersion(4, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutInitWindowPosition(100,250);
	glutInitWindowSize(RES_X, RES_Y);
	
	glDisable(GL_DEPTH_TEST);
	WindowHandle = glutCreateWindow(CAPTION);
	if(WindowHandle < 1) {
		std::cerr << "ERROR: Could not create a new rendering window." << std::endl;
		exit(EXIT_FAILURE);
	}

}


void init(int argc, char* argv[])
{
	// set the initial camera position on its spherical coordinates
	Eye  =  scene->GetCamera()->GetEye();
	camX = Eye.x;
	camY = Eye.y;
	camZ = Eye.z;
	r = Eye.length();
	_beta = asinf(camY/r) * 180.0f / 3.14f;
	alpha = atanf(camX / camZ) * 180.0f / 3.14f;

	setupGLUT(argc, argv);
	setupCallbacks();
	setupGLEW();
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	std::cerr << "CONTEXT: OpenGL v" << glGetString(GL_VERSION) << std::endl;
	createShaderProgram();
	createBufferObjects();
}

void init_scene(void)
{
	char scenes_dir[70] = "P3D_Scenes/";
	char input_user[50];
	char scene_name[70];


	scene = new Scene();

	if (P3F_scene) {  //Loading a P3F scene

		while (true) {
			cout << "Input the Scene Name: ";
			cin >> input_user;
			strcpy_s(scene_name, sizeof(scene_name), scenes_dir);
			strcat_s(scene_name, sizeof(scene_name), input_user);

			ifstream file(scene_name, ios::in);
			if (file.fail()) {
				printf("\nError opening P3F file.\n");
			}
			else
				break;
		}

		scene->load_p3f(scene_name);
		printf("Scene loaded.\n\n");
	}
	else {
		printf("Creating Peter Shirley Scene.\n\n");
		scene->create_random_scene();
	}

	RES_X = scene->GetCamera()->GetResX();
	RES_Y = scene->GetCamera()->GetResY();
	printf("\nResolutionX = %d  ResolutionY= %d.\n", RES_X, RES_Y);

	spp = scene->GetSamplesPerPixel();
	if (spp != 0) {
		AA = true; //Anti-aliasing
		printf("\nJittering with %d samples per pixel.\n", scene->GetSamplesPerPixel());
	}
	else AA = false;

	
	if ((scene->GetCamera()->GetAperture() != 0) && AA) {
		DOF = true; //Depth-Of-Field effect enabled only if AA enabled
		printf("Depth-Of-Field effect enabled\n");
	}
	else DOF = false;

	// Pixel buffer to be used in the Save Image function
	img_Data = (uint8_t*)malloc(3 * RES_X*RES_Y * sizeof(uint8_t));
	if (img_Data == NULL) exit(1);

	Accel_Struct = scene->GetAccelStruct();   //Type of acceleration data structure

	if (Accel_Struct == GRID_ACC) {
		grid_ptr = new Grid();
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		grid_ptr->Build(objs);
		printf("Grid built.\n\n");
	}
	else if (Accel_Struct == BVH_ACC) {
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();
		bvh_ptr = new BVH();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		bvh_ptr->Build(objs);
		//bvh_ptr->printNodes();
		printf("BVH built.\n\n");
	}
	else
		printf("No acceleration data structure.\n\n");

	unsigned int spp = scene->GetSamplesPerPixel();
	if (spp == 0)
		printf("Whitted Ray-Tracing\n");
	else
		printf("Distribution Ray-Tracing\n");
}

int main(int argc, char* argv[])
{
	//Initialization of DevIL 
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	int ch;
	if (!drawModeEnabled) {

		do {
			init_scene();
			
			auto timeStart = std::chrono::high_resolution_clock::now();
			renderScene();  //Just creating an image file
			auto timeEnd = std::chrono::high_resolution_clock::now();
			auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
			printf("\nDone: %.2f (sec)\n", passedTime / 1000);
			if (!P3F_scene) break;
			cout << "\nPress 'y' to render another image or another key to terminate!\n";
			delete(scene);
			if (Accel_Struct == GRID_ACC) delete(grid_ptr);
			else if (Accel_Struct == BVH_ACC) delete(bvh_ptr);
			free(img_Data);
			ch = _getch();
		} while((toupper(ch) == 'Y')) ;
	}

	else {   //Use OpenGL to draw image in the screen
		printf("OPENGL DRAWING MODE\n\n");
		init_scene();
		
	//	if (Accel_Struct == 1) grid_ptr->Build();  
		
		size_vertices = 2 * RES_X * RES_Y * sizeof(float);
		size_colors = 3 * RES_X * RES_Y * sizeof(float);
		vertices = (float*)malloc(size_vertices);
		if (vertices == NULL) exit(1);
		colors = (float*)malloc(size_colors);
		if (colors == NULL) exit(1);
		memset(colors, 0, size_colors);
		/* Setup GLUT and GLEW */
		init(argc, argv);
		glutMainLoop();
	}

	free(colors);
	free(vertices);
	printf("Programa terminado normalmente\n");
	exit(EXIT_SUCCESS);
}
///////////////////////////////////////////////////////////////////////