/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction

    vec3 eye_offset = cam.eye + cam.u * ls.x + cam.v * ls.y;
    float focal_ratio = cam.focusDist / cam.planeDist; 
    float px = ((pixel_sample.x / iResolution.x) - 0.5) * cam.width * focal_ratio;
    float py = ((pixel_sample.y / iResolution.y) - 0.5) * cam.height * focal_ratio;

    vec3 ray_direction = normalize(cam.u * (px - ls.x) + cam.v * (py - ls.y) - cam.n * cam.focusDist);
    
    return createRay(eye_offset, normalize(ray_direction), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIELECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for Dielectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDielectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIELECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};


float schlick(float cosine, float refIdx)
{
    float r0 = (1.0 - refIdx) / (1.0 + refIdx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    vec3 V = normalize(-rIn.d);
    vec3 outwardNormal = rec.normal;
    if(rec.material.type == MT_DIFFUSE)
    {
        //INSERT CODE HERE
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections
        //float VdotN = dot(V, rec.normal);
        //vec3 reflectDir = outwardNormal * VdotN * 2.0 - V;
        //reflectDir = normalize(reflectDir + randomInUnitSphere(gSeed) * rec.material.roughness);
        
        vec3 reflectDir = reflect(rIn.d, rec.normal); //calculate the reflected ray direction
        rScattered = createRay(rec.pos + rec.normal * epsilon, reflectDir);
        atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIELECTRIC) // fuzzy reflections and refractions
    {
        
        Ray rScattered;
        bool outside = dot(V, rec.normal) < 0.0; //check if the ray is outside or inside the object~
        vec3 refractDir = refract(rIn.d, rec.normal, outside ? rec.material.refIdx : 1.0 / rec.material.refIdx);
        float reflectProb;
        if(refractDir == vec3(0.0)) { //total internal reflection
            reflectProb = 1.0; //reflect always
        }else {
            //calculate the reflect probability using Schlick's approximation
            reflectProb = schlick(dot(V, rec.normal), rec.material.refIdx);
        }
        // Decide whether to reflect or refract
        if( hash1(gSeed) < reflectProb){ //Reflection
            // rScattered = calculate reflected ray
            rScattered = createRay(rec.pos + rec.normal * epsilon, reflect(rIn.d, rec.normal));
        }
        else { //Refraction
            if(!outside) {
                vec3 one = vec3(1.0, 1.0, 1.0); //white color
                atten = (one - rec.material.albedo * rec.t);
            } else {
                atten = vec3(1.0, 1.0, 1.0);
            }
            rScattered = createRay(rec.pos + rec.normal * epsilon, refract(rIn.d, rec.normal, outside? rec.material.refIdx : 1.0 / rec.material.refIdx));
        }
        return true;  
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    vec3 e1 = t.b - t.a;
    vec3 e2 = t.c - t.a;

    vec3 h = cross(r.d, e2);
    float a = dot(e1, h);
    float f = 1.0 / a;
    vec3 s = r.o - t.a;
    float u = f * dot(s, h);

    if(u < 0.0 || u > 1.0) return false;
    vec3 q = cross(s, e1);
    float v = f * dot(r.d, q);
    if(v < 0.0 || u + v > 1.0) return false;
    //calculate t
    float t_closest = f * dot(e2, q);

    vec3 normal = normalize(cross(e1, e2));

    //calculate a valid t and normal
    if(t_closest < tmax && t_closest > tmin)
    {
        rec.t = t_closest;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}


struct Quad {vec3 a; vec3 b; vec3 c; vec3 d; };

Quad createQuad(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{
    Quad q;
    q.a = v0; q.b = v1; q.c = v2; q.d = v3;
    return q;
}

bool hit_quad(Quad q, Ray r, float tmin, float tmax, out HitRecord rec)
{
    if(hit_triangle(createTriangle(q.a, q.b, q.c), r, tmin, rec.t, rec)) return true;
    else if(hit_triangle(createTriangle(q.a, q.c, q.d), r, tmin, rec.t, rec)) return true;
    else return false;  
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
	//Program it
    vec3 moving_center = mvsphere.center0 + (mvsphere.center1 - mvsphere.center0) * (time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0);
    return moving_center;
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //INSERT YOUR CODE HERE
    //calculate a valid t and normal
    bool hit = false;
    vec3 oc = r.o - s.center;
    float a = dot(r.d, r.d);
    float b = 2.0f *  dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;

    float discriminant = b * b - 4.0f * a * c;
    
    if(discriminant < 0.0) return false; //no intersection

    float sqrtD = sqrt(discriminant);
    float t1 = (-b - sqrtD) / (2.0f * a); //first root
    float t2 = (-b + sqrtD) / (2.0f * a); //second root
    float t;
    
    if(t1 > epsilon) {
        t = t1; //use the first root if it is positive
    } else if(t2 > epsilon) {
        t = t2; //use the second root if the first is negative
    } else {
        return false; //both roots are negative, no intersection
    }
     
    vec3 normal = normalize((pointOnRay(r, t) - s.center)); //calculate the normal at the intersection point
	
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        hit = true;
    }
        return hit;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    float B, C, delta;
    bool outside;
    float t;

     //Calculate the moving center
    vec3 center0 = center(s, r.t);
    vec3 center1 = center(s, r.t + epsilon);
    vec3 oc = r.o - center0;
    vec3 oc1 = r.o - center1;
    vec3 d = r.d; //ray direction
    vec3 d1 = (center1 - center0) / (s.time1 - s.time0); //velocity vector

    float a = dot(d, d) - dot(d, d1) * dot(d, d1); 
    float b = 2.0 * (dot(oc, d) - dot(oc, d1) * dot(d, d1));
    float c = dot(oc, oc) - dot(oc, d1) * dot(oc, d1) - s.radius * s.radius; 

    //calculate the discriminant  
    B = b / (2.0 * a);
    C = c - B * B;      
    delta = B * B - C;

    if(delta < 0.0) return false; //no intersection

    float sqrtD = sqrt(delta);
    
    
    //calculate a valid t
    if(B - sqrtD > epsilon) {
        t = B - sqrtD; //use the first root if it is positive
        outside = true;
    } else if(B + sqrtD > epsilon) {
        t = B + sqrtD; //use the second root if the first is negative
        outside = false;
    } else {
        return false; //both roots are negative, no intersection
    }

    //calculate the normal at the intersection point
    //if outside, normal is from center0, if inside, normal is from center1
    vec3 normal;
    if(outside) {
        normal = normalize((pointOnRay(r, t) - center0)); //calculate the normal at the intersection point
    } else {
        normal = normalize((pointOnRay(r, t) - center1)); //calculate the normal at the intersection point
    }
    

    //check if the intersection is within the valid range
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        return true;
    }
    else return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}