/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
*/

#include "./common.glsl"
#iChannel0 "self"
#iChannel1 "file://cubemaps/forest_{}.png" // Use a single wildcard for CubeMap
#iChannel1::Type "CubeMap"
#iKeyboard
 
#define SCENE 2
#define SCENE 2

bool hit_world(Ray r, float tmin, float tmax, inout HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;

    #if SCENE == 0       //Shirley Weekend scene

        if(hit_quad(createQuad(vec3(-10.0, -0.05, 10.0), vec3(10.0, -0.05, 10.0), vec3(10.0, -0.05, -10.0), vec3(-10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2));
        }

        if(hit_sphere(createSphere(vec3(-4.0, 1.0, 0.0), 1.0), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
            //rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
        }

        if(hit_sphere(createSphere(vec3(4.0, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            //rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
            rec.material = createMetalMaterial(vec3(0.562, 0.565, 0.578), 0.6);
            //rec.material = createPlasticMaterial(vec3(0.0, 0.5, 1.0), 0.0);
            //rec.material = createPlasticMaterial(vec3(1.0, 0.0, 1.0), 0.4);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), -0.5),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0, 0.9, 0.9), 1.5, 0.0);
        }
            
        int numxy = 5;
        
        for(int x = -numxy; x < numxy; ++x)
        {
            for(int y = -numxy; y < numxy; ++y)
            {
                float fx = float(x);
                float fy = float(y);
                float seed = fx + fy / 1000.0;
                vec3 rand1 = hash3(seed);
                vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
                float chooseMaterial = rand1.z;
                if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
                {
                    if(chooseMaterial < 0.3)
                    {
                        vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                        // diffuse
                        if(hit_movingSphere(createMovingSphere(center, center1, 0.2, 0.0, 1.0),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.5)
                    {
                        // diffuse
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.7)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                        }
                    }
                    else if(chooseMaterial < 0.9)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                        }
                    }
                    else
                    {
                        // glass (Dielectric)
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDielectricMaterial(hash3(seed), 1.33, 0.0);
                        }
                    }
                }
            }
        }
    #elif SCENE == 1 //from https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/

        // diffuse floor
        
        vec3 A = vec3(-25.0f, -12.5f, 10.0f);
        vec3 B = vec3( 25.0f, -12.5f, 10.0f);
        vec3 C = vec3( 25.0f, -12.5f, -5.0f);
        vec3 D = vec3(-25.0f, -12.5f, -5.0f);

        if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.7));
        }

        //stripped background
        {
            vec3 A = vec3(-25.0f, -10.5f, -5.0f);
            vec3 B = vec3( 25.0f, -10.5f, -5.0f);
            vec3 C = vec3( 25.0f, -1.5f, -5.0f);
            vec3 D = vec3(-25.0f, -1.5f, -5.0f);
        
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                float shade = floor(mod(rec.pos.x, 1.0f) * 2.0f);
                rec.material = createDiffuseMaterial(vec3(shade));
            }
        }

        // ceiling piece above light
        
        {
            vec3 A = vec3(-7.5f, 12.5f, 5.0f);
            vec3 B = vec3( 7.5f, 12.5f, 5.0f);
            vec3 C = vec3( 7.5f, 12.5f, -5.0f);
            vec3 D = vec3(-7.5f, 12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }
        }    
       
        // light
        
        {
            vec3 A = vec3(-5.0f, 12.3f,  2.5f);
            vec3 B = vec3( 5.0f, 12.3f,  2.5f);
            vec3 C = vec3( 5.0f, 12.3f,  -2.5f);
            vec3 D = vec3(-5.0f, 12.3f,  -2.5f);

             if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.9f) * 20.0f;
            }
        }
 
        const int c_numSpheres = 7;
        for (int sphereIndex = 0; sphereIndex < c_numSpheres; ++sphereIndex)
        {
            vec3 center = vec3(-18.0 + 6.0 * float(sphereIndex), -8.0, 0.0);
            if(hit_sphere(createSphere(center, 2.8),r,tmin,rec.t,rec))
            {
                hit = true;
                float r = float(sphereIndex) / float(c_numSpheres-1) * 0.1f;
                rec.material = createDielectricMaterial(vec3(0.0, 0.5, 1.0), 1.1, r);
            }
        }

    #elif SCENE == 2
        // diffuse floor
        
        vec3 A = vec3(-25.0f, -12.5f, 10.0f);
        vec3 B = vec3( 25.0f, -12.5f, 10.0f);
        vec3 C = vec3( 25.0f, -12.5f, -5.0f);
        vec3 D = vec3(-25.0f, -12.5f, -5.0f);

        if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.7));
        }

        //stripped background
        {
            vec3 A = vec3(-25.0f, -10.5f, -5.0f);
            vec3 B = vec3( 25.0f, -10.5f, -5.0f);
            vec3 C = vec3( 25.0f, -1.5f, -5.0f);
            vec3 D = vec3(-25.0f, -1.5f, -5.0f);
        
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                float shade = floor(mod(rec.pos.x, 1.0f) * 2.0f);
                rec.material = createDiffuseMaterial(vec3(shade));
            }
        }

        // ceiling piece above light
        
        {
            vec3 A = vec3(-7.5f, 12.5f, 5.0f);
            vec3 B = vec3( 7.5f, 12.5f, 5.0f);
            vec3 C = vec3( 7.5f, 12.5f, -5.0f);
            vec3 D = vec3(-7.5f, 12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }
        }    
       
        // light
        
        {
            vec3 A = vec3(-5.0f, 12.3f,  2.5f);
            vec3 B = vec3( 5.0f, 12.3f,  2.5f);
            vec3 C = vec3( 5.0f, 12.3f,  -2.5f);
            vec3 D = vec3(-5.0f, 12.3f,  -2.5f);

             if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.9f) * 20.0f;
            }
        }
 
        const int c_numSpheres = 7;
        for (int sphereIndex = 0; sphereIndex < c_numSpheres; ++sphereIndex)
        {
            vec3 center = vec3(-18.0 + 6.0 * float(sphereIndex), -8.0, 0.0);
            if(hit_sphere(createSphere(center, 2.8),r,tmin,rec.t,rec))
            {
                hit = true;
                float r = float(sphereIndex) / float(c_numSpheres-1) * 0.1f;
                rec.material = createPlasticMaterial(vec3(0.9, 0.25, 0.25), r);
            }
        }
    #elif SCENE == 3

       	// back wall
        {
            vec3 A = vec3(-12.6f, -12.6f, 25.0f);
            vec3 B = vec3( 12.6f, -12.6f, 25.0f);
            vec3 C = vec3( 12.6f,  12.6f, 25.0f);
            vec3 D = vec3(-12.6f,  12.6f, 25.0f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7f, 0.7f, 0.7f));
            }
        }    
        
        // floor
        {
            vec3 A = vec3(-12.6f, -12.45f, 25.0f);
            vec3 B = vec3( 12.6f, -12.45f, 25.0f);
            vec3 C = vec3( 12.6f, -12.45f, 15.0f);
            vec3 D = vec3(-12.6f, -12.45f, 15.0f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7f, 0.7f, 0.7f));
            }
        }
        // ceiling
        {
            vec3 A = vec3(-12.6f, 12.5f, 25.0f);
            vec3 B = vec3( 12.6f, 12.5f, 25.0f);
            vec3 C = vec3( 12.6f, 12.5f, 15.0f);
            vec3 D = vec3(-12.6f, 12.5f, 15.0f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7f, 0.7f, 0.7f));
            }       
        }    
        
        // left wall
        {
            vec3 A = vec3(-12.5f, -12.6f, 25.0f);
            vec3 B = vec3(-12.5f, -12.6f, 15.0f);
            vec3 C = vec3(-12.5f,  12.6f, 15.0f);
            vec3 D = vec3(-12.5f,  12.6f, 25.0f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7f, 0.1f, 0.1f));
            }      
        }
        
        // right wall 
        {
            vec3 A = vec3( 12.5f, -12.6f, 25.0f);
            vec3 B = vec3( 12.5f, -12.6f, 15.0f);
            vec3 C = vec3( 12.5f,  12.6f, 15.0f);
            vec3 D = vec3( 12.5f,  12.6f, 25.0f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.1f, 0.7f, 0.1f));
            }      
        }    
        
        // light
        {
            vec3 A = vec3(-5.0f, 12.4f,  22.5f);
            vec3 B = vec3( 5.0f, 12.4f,  22.5f);
            vec3 C = vec3( 5.0f, 12.4f,  17.5f);
            vec3 D = vec3(-5.0f, 12.4f,  17.5f);
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.7f) * 20.0f;
            }
    
        }
        if(hit_sphere(createSphere(vec3(-9.0f, -9.5f, 20.0f), 3.0f),r,tmin,rec.t,rec))
        {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.9f, 0.9f, 0.75f));
        }

        if(hit_sphere(createSphere(vec3(0.0f, -9.5f, 20.0f), 3.0f),r,tmin,rec.t,rec))
        {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.9f, 0.75f, 0.9f));
        }

        if(hit_sphere(createSphere(vec3(9.0f, -9.5f, 20.0f), 3.0f),r,tmin,rec.t,rec))
        {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.75f, 0.9f, 0.9f));
        }
   
    #elif SCENE == 4  // Pool scene
        // Pool floor (concrete)
        if (hit_quad(createQuad(vec3(-10.0, -5.0, 10.0), vec3(10.0, -5.0, 10.0), vec3(10.0, -5.0, -10.0), vec3(-10.0, -5.0, -10.0)), r, tmin, rec.t, rec)) {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.7, 0.7, 0.7));
        }

        // Water surface (with slight waves)
        if (hit_quad(createQuad(vec3(-10.0, 0.0, 10.0), vec3(10.0, 0.0, 10.0), vec3(10.0, 0.0, -10.0), vec3(-10.0, 0.0, -10.0)), r, tmin, rec.t, rec)) {
            hit = true;
            // rec.normal += 0.1 * hash3(gSeed); // Fake waves
            // rec.normal = normalize(rec.normal);
            rec.material = createDielectricMaterial(vec3(0.8, 0.9, 1.0), 1.33, 0.1);
        }

        // Red ball (half-submerged)
        if (hit_sphere(createSphere(vec3(0.0, -1.0, 0.0), 1.0), r, tmin, rec.t, rec)) {
            hit = true;
            rec.material = createPlasticMaterial(vec3(1.0, 0.0, 0.0), 0.2);
        }
        
    #endif

    return hit;
}

vec3 directlighting(quadLight l, Ray r, HitRecord rec){
    vec3 diffCol, specCol;
    vec3 colorOut = vec3(0.0, 0.0, 0.0);
    float shininess;
    HitRecord dummy;
    vec3 N = normalize(rec.normal);

    // light position in world coordinates
    vec3 lightpos = l.pos + l.e1 * hash2(gSeed).x + l.e2 * hash2(gSeed).y;

    // 1. calculate the direction to the light source
    vec3 lightDir = normalize(l.pos - rec.pos);

    if (dot(N, lightDir) > 0.0){


        // Diffuse
        vec3 diffCol = rec.material.albedo * max(dot(N, lightDir), 0.0);

        // Specular
        vec3 viewDir = normalize(-r.d);
        vec3 H = normalize(lightDir + viewDir);
        float shininess = 8.0 / (pow(rec.material.roughness, 4.0) + epsilon) - 2.0;
        vec3 specCol = rec.material.specColor * pow(max(dot(N, H), 0.0), shininess);

        vec3 F0 = rec.material.specColor;

        vec3 ks = fresnelSchlick(max(dot(N, -viewDir), 0.0), F0);
        vec3 kd = vec3(1.0) - ks;

        if(rec.material.type == MT_METAL || rec.material.type == MT_PLASTIC)
            specCol = BRDF_GGX(N, -viewDir, lightDir, F0, rec.material.roughness);
        if (rec.material.type == MT_PLASTIC)
            diffCol = kd * rec.material.albedo / pi;
        // 4. combine contributions and apply attenuation
        colorOut += (diffCol + specCol) * l.color * max(dot(N, lightDir), 0.0);
    }
    
	return colorOut; 
}

vec3 directlighting(pointLight l, Ray r, HitRecord rec){
    vec3 diffCol, specCol;
    vec3 colorOut = vec3(0.0, 0.0, 0.0);
    float shininess;
    HitRecord dummy;
    vec3 N = normalize(rec.normal);

    // light position in world coordinates
    vec3 lightpos = l.pos;

    // 1. calculate the direction to the light source
    vec3 lightDir = normalize(l.pos - rec.pos);

    // Direction and distance to light
    float distToLight = length(lightDir);

    // HARD SHADOWS: Cast a ray to the light
    Ray shadowRay = createRay(rec.pos + N * 0.001, lightDir);
    HitRecord shadowHit;
    if (hit_world(shadowRay, 0.001, distToLight - 0.001, shadowHit)) {
        return vec3(0.0); // In shadow: no light contribution
    }

    if (dot(N, lightDir) > 0.0){
        // 2. calculate the diffuse color contribution
        diffCol = rec.material.albedo * max(dot(N, lightDir), 0.0);
        
        // 3. calculate the specular color contribution
        vec3 viewDir = normalize(r.d);
        vec3 H = normalize(lightDir - viewDir); // half vector
        shininess = 8.0 / (pow(rec.material.roughness, 4.0)+epsilon) - 2.0;
        specCol = rec.material.specColor * pow(max(dot(N, H), 0.0), shininess);

        vec3 F0 = rec.material.specColor;

        vec3 ks = fresnelSchlick(max(dot(N, -viewDir), 0.0), F0);
        vec3 kd = vec3(1.0) - ks;

        if(rec.material.type == MT_METAL || rec.material.type == MT_PLASTIC)
            specCol = BRDF_GGX(N, -viewDir, lightDir, F0, rec.material.roughness);
        if (rec.material.type == MT_PLASTIC)
            diffCol = kd * rec.material.albedo / pi;
        // 4. combine contributions and apply attenuation
        colorOut += (diffCol + specCol) * l.color * max(dot(N, lightDir), 0.0);
    }
    
	return colorOut; 
}


#define MAX_BOUNCES 10

vec3 rayColor(Ray r)
{
    HitRecord rec;
    vec3 col = vec3(0.0);
    vec3 throughput = vec3(1.0f, 1.0f, 1.0f);
    for(int i = 0; i < MAX_BOUNCES; ++i)
    {
        if(hit_world(r, 0.001, 10000.0, rec))
        {

            if(rec.material.emissive != vec3(0.0))
            {
                // if the material is emissive, add its color to the output
                col +=  rec.material.emissive * throughput;
            }
            if(SCENE == 0){  

                pointLight l1 = createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0, 1.0, 1.0));
                pointLight l2 = createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0, 1.0, 1.0));
                pointLight l3 = createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0, 1.0, 1.0));

                col += directlighting(l1, r, rec) * throughput;
                col += directlighting(l2, r, rec) * throughput;
                col += directlighting(l3, r, rec) * throughput;

            } else if(SCENE == 1){

                quadLight l4 = createQuadLight(vec3( 5.0f, 12.3f,  2.5f), vec3(1.0, 1.0, 1.0), vec3( -5.0f, 12.3f,  2.5f), vec3( 5.0f, 12.3f,  -2.5f));
                col += directlighting(l4, r, rec) * throughput;     

            }else if(SCENE == 2){

                quadLight l4 = createQuadLight(vec3( 5.0f, 12.3f,  2.5f), vec3(1.0, 1.0, 1.0), vec3( -5.0f, 12.3f,  2.5f), vec3( 5.0f, 12.3f,  -2.5f));
                col += directlighting(l4, r, rec) * throughput;
            
            }else if(SCENE == 3){
                quadLight l4 = createQuadLight(vec3( 5.0f, 12.3f,  2.5f), vec3(1.0, 1.0, 1.0), vec3( -5.0f, 12.3f,  2.5f), vec3( 5.0f, 12.3f,  -2.5f));
                col += directlighting(l4, r, rec) * throughput; 
            } else if(SCENE == 4){
                // Pool lights (optional)
                pointLight light1 = createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0, 1.0, 1.0) * 10.0);
                //pointLight light2 = createPointLight(vec3(-5.0, 4.0, 5.0), vec3(1.0, 1.0, 1.0) * 10.0);
                col += directlighting(light1, r, rec) * throughput;
                //col += directlighting(light2, r, rec) * throughput;
            }
            

            //calculate secondary ray and update throughput

            Ray scatterRay;
            vec3 atten = vec3(1.0f, 1.0f, 1.0f);

            if(scatter(r, rec, atten, scatterRay))
            {   
                r = scatterRay;
                throughput *= atten;

            }
            else
            {
                col += throughput * rec.material.emissive;
                break; // no more scattering, exit loop
            }

            // Russian Roulette
            // As the throughput gets smaller, the ray is more likely to get terminated early.
            // Survivors have their value boosted to make up for fewer samples being in the average.
            
            float p = max(throughput.r, max(throughput.g, throughput.b));
            if (hash1(gSeed) > p)
                break;

            // Add the energy we 'lose' by randomly terminating paths
            throughput *= 1.0f / p;            
            
        }
        else
        {
            float t = 0.8 * (r.d.y + 1.0);
            if(SCENE == 0){
                col += throughput * mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            } else if(SCENE == 1 || SCENE == 2 || SCENE == 3){
                // Background color for other scenes
                col += throughput * SRGBToLinear(texture(iChannel1, r.d).rgb);
            } else if(SCENE == 4){
                // Pool scene background color
                col += throughput * mix(vec3(0.2, 0.4, 0.6), vec3(0.1, 0.2, 0.3), t);
            }
            break;
        }
    }
    return col; 
}

// Constants for controlling zoom and camera range
const float zoomSpeed = 0.2f; // Adjust the zoom speed (how fast the zoom happens)

// Mouse camera control parameters
const float c_minCameraAngle = 0.01f;
const float c_maxCameraAngle = (pi - 0.01f);



void GetCameraVectors(out vec3 cameraPos, out vec3 cameraFwd, out vec3 cameraUp, out vec3 cameraRight)
{
    float minZoom = 5.0f; // Minimum camera distance
    float maxZoom; // Maximum camera distance
    vec3 c_cameraAt = vec3(iMouse.xy.x, iMouse.xy.y, 0.0f); // Default camera target position
    if(SCENE == 0){
        c_cameraAt = vec3(0.0f, 0.0f, 1.0f);
        maxZoom = 9.0f; // Set maximum zoom for Shirley Weekend scene
    } else if(SCENE == 1 || SCENE == 2 || SCENE == 3 || SCENE == 4) {
        minZoom = 10.0f;
        c_cameraAt = vec3(0.0f, -3.0f, 10.0f);
        maxZoom = 40.0f; // Set maximum zoom for the other scene
    } 
    // Get mouse scroll input to simulate zoom (iMouse.z controls zooming)
    float scroll = iMouse.z; // Positive for zooming in, negative for zooming out

    // Dynamic camera distance adjustment based on scroll input
    float cameraDistance = 5.0f; // Default camera distance
    cameraDistance = clamp(cameraDistance - scroll * zoomSpeed, minZoom, maxZoom); // Clamp to ensure it stays within the set range
    // if the mouse is at (0,0) it hasn't been moved yet, so use a default camera setup
    vec2 mouse = iMouse.xy;
    if (dot(mouse, vec2(1.0f, 1.0f)) == 0.0f)
    {
        cameraPos = vec3(0.0f, 0.0f, -cameraDistance); // Use adjusted cameraDistance for zoom
        cameraFwd = vec3(0.0f, 0.0f, 1.0f);
        cameraUp = vec3(0.0f, 1.0f, 0.0f);
        cameraRight = vec3(1.0f, 0.0f, 0.0f);
        return;
    }
     
    // Calculate the camera position using mouse movement (orbiting)
    float angleX = -mouse.x * 16.0f / float(iResolution.x);
    float angleY = mix(c_minCameraAngle, c_maxCameraAngle, mouse.y / float(iResolution.y));
     
    cameraPos.x = sin(angleX) * sin(angleY) * cameraDistance;  // Use adjusted cameraDistance for zoom
    cameraPos.y = -cos(angleY) * cameraDistance;
    cameraPos.z = cos(angleX) * sin(angleY) * cameraDistance;
     
    cameraPos += c_cameraAt; // Offset by target position
    
    // Calculate forward, right, and up vectors
    cameraFwd = normalize(c_cameraAt - cameraPos);
    cameraRight = normalize(cross(vec3(0.0f, 1.0f, 0.0f), cameraFwd));
    cameraUp = normalize(cross(cameraFwd, cameraRight));   
}



#define MAX_SAMPLES 10000.0


void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

        // calculate subpixel camera jitter for anti aliasing
 
    // get the camera vectors
    vec3 cameraPos, cameraFwd, cameraUp, cameraRight;
    GetCameraVectors(cameraPos, cameraFwd, cameraUp, cameraRight);    

    vec3 camPos = cameraPos; // small jitter for anti-aliasing
    vec3 camTarget = cameraPos + cameraFwd;
    float fovy = 60.0;
    float aperture = 0.0;
    float distToFocus = 1.0;
    float time0 = 0.0;
    float time1 = 1.0;
    Camera cam = createCamera(
        camPos,
        camTarget,
        cameraUp,
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

//usa-se o 4 canal de cor para guardar o numero de samples e não o iFrame pois quando se mexe o rato faz-se reset

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));

    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.0;
    color = mix(prevLinear, color, 1.0/w);
    gl_FragColor = vec4(toGamma(color), w);
}
