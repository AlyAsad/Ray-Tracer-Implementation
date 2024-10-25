#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>




//global vars and function declarations
int depth;
std::vector<parser::Vec3f> triNormals;
std::vector<std::vector<parser::Vec3f>> meshNormals;

parser::Vec3f applyShading(parser::Scene &scene, parser::Vec3f &e, parser::Vec3f &d, int &mat_id, parser::Vec3f &point, parser::Vec3f &normal);


//dot product helper function
float dot(parser::Vec3f &a, parser::Vec3f &b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

//cross product helper function
parser::Vec3f cross(parser::Vec3f &a, parser::Vec3f &b) {
    parser::Vec3f u;
    u.x = (a.y * b.z) - (a.z * b.y);
    u.y = (a.z * b.x) - (a.x * b.z);
    u.z = (a.x * b.y) - (a.y * b.x);
    return u;
}

//helper to find point of intersection
parser::Vec3f pointOfIntersection(parser::Vec3f &e, parser::Vec3f &d, float &t) {
    parser::Vec3f point;
    point.x = e.x + (t*d.x);
    point.y = e.y + (t*d.y);
    point.z = e.z + (t*d.z);
    return point;
}

//helper to find normal of sphere with point of intersection
parser::Vec3f sphereNormal(parser::Vec3f &p, parser::Sphere &sphere, std::vector<parser::Vec3f> &vertex_data) {
    parser::Vec3f norm;
    norm.x = (p.x - vertex_data[sphere.center_vertex_id -1].x)/sphere.radius;
    norm.y = (p.y - vertex_data[sphere.center_vertex_id -1].y)/sphere.radius;
    norm.z = (p.z - vertex_data[sphere.center_vertex_id -1].z)/sphere.radius;
    return norm;
}


//Barycentric Coordinates
float intersection(parser::Vec3f &a, parser::Vec3f &b, parser::Vec3f &c, parser::Vec3f &d, parser::Vec3f &e)
{
	float m11 = a.x-b.x, m12 = a.x - c.x, m13 = d.x;
	float m21 = a.y-b.y, m22 = a.y - c.y, m23 = d.y;
	float m31 = a.z-b.z, m32 = a.z - c.z, m33 = d.z;
	
	
	float detA = m11*(m22*m33 - m23*m32) - m21*(m12*m33 - m13*m32) + m31*(m12*m23 - m13*m22);
	
	if (detA == 0) return -1;
	
	float cx = a.x-e.x, cy = a.y-e.y, cz = a.z-e.z;
	float beta = (cx*(m22*m33 - m23*m32) - cy*(m12*m33 - m13*m32) + cz*(m12*m23 - m13*m22)) / detA;	 

	if(beta < 0 || beta > 1) return -1;
	
	float gamma = (m11*(cy*m33 - m23*cz) - m21*(cx*m33 - m13*cz) + m31*(cx*m23 - m13*cy)) / detA;
	if(gamma < 0 || beta + gamma > 1) return -1;

	return (m11*(m22*cz - cy*m32) - m21*(m12*cz - cx*m32) + m31*(m12*cy - cx*m22)) / detA;
}



//helper to find if intersecting any object for shadow
bool inShadow(parser::Scene &scene, parser::Vec3f &e, parser::Vec3f &d) {
    float curr_t;
        // copy pasted code from below to check if impacting any object
    for (auto sphere : scene.spheres) { //checking intersection with spheres first
        parser::Vec3f c = scene.vertex_data[sphere.center_vertex_id -1];
        c.x = e.x - c.x;
        c.y = e.y - c.y;
        c.z = e.z - c.z;
        float dd = dot(d,d), dc = dot(d,c), cc = dot(c,c);
        float disc = pow(dc,2) - (dd * (cc - pow(sphere.radius,2)));
        if (disc == 0) { //one intersection
            curr_t = (-dc)/dd;
            if (curr_t > 0 && curr_t < 1) return true;
        }
        
        else if (disc > 0) { //two intersections
            curr_t = (-dc + sqrt(disc))/dd;
            if (curr_t > 0 && curr_t < 1) return true;
            curr_t = (-dc - sqrt(disc))/dd;
            if (curr_t > 0 && curr_t < 1) return true;
        }
    }
    
    long int size = scene.triangles.size();
    for (int i = 0; i < size; i++) { //checking intersection with triangles second
        
        if (dot(triNormals[i],d) == 0) continue; //checking if parallel(=0)
        
        parser::Triangle currTri = scene.triangles[i];
        
        curr_t = intersection(scene.vertex_data[currTri.indices.v0_id - 1], scene.vertex_data[currTri.indices.v1_id - 1], scene.vertex_data[currTri.indices.v2_id - 1], d, e);
        //finally found t, comparing with lowest t and if lesser replacing lowest one
        if (curr_t > 0 && curr_t < 1) return true;
    }
    
    size = scene.meshes.size();
    for (int i = 0; i < size; i++) { //checking intersection with meshes final
        long int size2 = scene.meshes[i].faces.size();
        
    for (int j = 0; j < size2; j++) { // iterate over triangles of each mesh
        
        if (dot(meshNormals[i][j],d) == 0) continue; //checking if parallel(=0)
    
        parser::Face currFace = (scene.meshes[i]).faces[j];
        
       curr_t = intersection(scene.vertex_data[currFace.v0_id - 1], scene.vertex_data[currFace.v1_id - 1], scene.vertex_data[currFace.v2_id - 1], d, e);
        //finally found t, comparing with lowest t and if lesser replacing lowest one
        if (curr_t > 0 && curr_t < 1) return true;
            
    }
    }
    
    return false;
}







//recursive function
parser::Vec3f computeColor(parser::Scene &scene, parser::Vec3f &e, parser::Vec3f &d) {
    
    //if max depth exceeded
    if (depth > scene.max_recursion_depth) return {0,0,0};
    
    //find closest hit point
    
    float low_t = INFINITY, curr_t;
    int mat_id;
    parser::Vec3f point, normal;
    
    for (auto sphere : scene.spheres) { //checking intersection with spheres first
        parser::Vec3f c = scene.vertex_data[sphere.center_vertex_id -1];
        c.x = e.x - c.x;
        c.y = e.y - c.y;
        c.z = e.z - c.z;
        float dd = dot(d,d), dc = dot(d,c), cc = dot(c,c);
        float disc = pow(dc,2) - (dd * (cc - pow(sphere.radius,2)));
        if (disc == 0) { //one intersection
            curr_t = (-dc)/dd;
            if (curr_t > 0 && curr_t < low_t) {
                low_t = curr_t;
                mat_id = sphere.material_id;
                point = pointOfIntersection(e, d, curr_t);
                normal = sphereNormal(point, sphere, scene.vertex_data);
            }
        }
        
        else if (disc > 0) { //two intersections
            curr_t = (-dc + sqrt(disc))/dd;
            if (curr_t > 0 && curr_t < low_t) {
                low_t = curr_t;
                mat_id = sphere.material_id;
                point = pointOfIntersection(e, d, curr_t);
                normal = sphereNormal(point, sphere, scene.vertex_data);
            }
            curr_t = (-dc - sqrt(disc))/dd;
            if (curr_t > 0 && curr_t < low_t) {
                low_t = curr_t;
                mat_id = sphere.material_id;
                point = pointOfIntersection(e, d, curr_t);
                normal = sphereNormal(point, sphere, scene.vertex_data);
            }
        }
    }
    
    long int size = scene.triangles.size();
    for (int i = 0; i < size; i++) { //checking intersection with triangles second
        
        if (dot(triNormals[i],d) >= 0) continue; //parallel(=0) or facing away(>0): Back-face culling
        
        parser::Triangle currTri = scene.triangles[i];
        
        curr_t = intersection(scene.vertex_data[currTri.indices.v0_id - 1], scene.vertex_data[currTri.indices.v1_id - 1], scene.vertex_data[currTri.indices.v2_id - 1], d, e);
        //found t, comparing with lowest t and if lesser replacing lowest one
        if (curr_t > 0 && curr_t < low_t) {
            low_t = curr_t;
            mat_id = currTri.material_id;
            point = pointOfIntersection(e, d, curr_t);
            normal = triNormals[i];
        }
    }
    
    size = scene.meshes.size();
    for (int i = 0; i < size; i++) { //checking intersection with meshes final
        long int size2 = scene.meshes[i].faces.size();
        
    for (int j = 0; j < size2; j++) { // iterate over triangles of each mesh
        
        if (dot(meshNormals[i][j],d) >= 0) continue; //parallel(=0) or facing away(>0): Back-face culling
    
        parser::Face currFace = (scene.meshes[i]).faces[j];
        
        curr_t = intersection(scene.vertex_data[currFace.v0_id - 1], scene.vertex_data[currFace.v1_id - 1], scene.vertex_data[currFace.v2_id - 1], d, e);
        //finally found t, comparing with lowest t and if lesser replacing lowest one
        if (curr_t > 0 && curr_t < low_t) {
            low_t = curr_t;
            mat_id = scene.meshes[i].material_id;
            point = pointOfIntersection(e, d, curr_t);
            normal = meshNormals[i][j];
        }
            
    }
    }
    
    //if found hit point
    if (low_t != INFINITY) return applyShading(scene, e, d, mat_id, point, normal);
    
    
    //if no hit and primary ray
    else if (depth == 0) {
    parser::Vec3f out;
    out.x = scene.background_color.x;
    out.y = scene.background_color.y;
    out.z = scene.background_color.z;
    return out;
    }
    
    //if reflected from mirror and no intersection
    else return {0,0,0};
}














//shading recursive function
parser::Vec3f applyShading(parser::Scene &scene, parser::Vec3f &e, parser::Vec3f &d, int &mat_id, parser::Vec3f &point, parser::Vec3f &normal) {

    //definitions
    parser::Vec3f color;
    parser::Material mat = scene.materials[mat_id-1];
    
    parser::Vec3f oPoint; //offset point
    oPoint.x = point.x + (scene.shadow_ray_epsilon * normal.x);
    oPoint.y = point.y + (scene.shadow_ray_epsilon * normal.y);
    oPoint.z = point.z + (scene.shadow_ray_epsilon * normal.z);
    
    
    //adding ambient light
    color.x = mat.ambient.x * scene.ambient_light.x;
    color.y = mat.ambient.y * scene.ambient_light.y;
    color.z = mat.ambient.z * scene.ambient_light.z;
    
    
    //checking for mirror and adding mirror light
    if (mat.is_mirror) {
        parser::Vec3f refd;
        refd.x = -(2 * normal.x * dot(normal, d)) + d.x;
        refd.y = -(2 * normal.y * dot(normal, d)) + d.y;
        refd.z = -(2 * normal.z * dot(normal, d)) + d.z;
        depth++;
        parser::Vec3f colorMir = computeColor(scene, oPoint, refd);
        color = {color.x + (mat.mirror.x*colorMir.x), color.y + (mat.mirror.y*colorMir.y), color.z + (mat.mirror.z*colorMir.z)};
    }
    
    
    //diffuse and specular shading using point lights
    
    for (auto light : scene.point_lights) {
        
        parser::Vec3f w_i = {light.position.x - oPoint.x, light.position.y - oPoint.y, light.position.z - oPoint.z};
        
        //checking if light parallel or hitting from behind
       float dotP = dot(w_i, normal);
        if (dotP <= 0) continue;
        
        if (inShadow(scene, oPoint, w_i)) continue; //if in shadow dont calculate diffuse and specular
        
        float distSq = dot(w_i,w_i); //normalizing w_i
        float dist = sqrt(distSq);
        w_i.x /= dist;
        w_i.y /= dist;
        w_i.z /= dist;
        
        //calculating diffuse
        dist = dotP/dist;
        color.x += (mat.diffuse.x * dist * light.intensity.x)/distSq;
        color.y += (mat.diffuse.y * dist * light.intensity.y)/distSq;
        color.z += (mat.diffuse.z * dist * light.intensity.z)/distSq;
        
        
        //calculating specular
        parser::Vec3f h = {w_i.x - d.x, w_i.y - d.y, w_i.z - d.z};
        dist = sqrt(dot(h,h));
        h = {h.x/dist, h.y/dist, h.z/dist};
        dist = pow(std::max((float)0, dot(h, normal)), mat.phong_exponent);
        
        color.x += (mat.specular.x * dist * light.intensity.x)/distSq;
        color.y += (mat.specular.y * dist * light.intensity.y)/distSq;
        color.z += (mat.specular.z * dist * light.intensity.z)/distSq;
        
    }
    
    return color;
}















int main(int argc, char* argv[])
{

    
    //loading scene from Xml file
    parser::Scene scene;
    scene.loadFromXml(argv[1]);
    float tempf;
    
    //prelim finding normals of triangles and meshes
    for (auto tri : scene.triangles) {
        parser::Vec3f temp1, temp2;
        temp1.x = scene.vertex_data[tri.indices.v2_id -1].x - scene.vertex_data[tri.indices.v1_id -1].x;
        temp1.y = scene.vertex_data[tri.indices.v2_id -1].y - scene.vertex_data[tri.indices.v1_id -1].y;
        temp1.z = scene.vertex_data[tri.indices.v2_id -1].z - scene.vertex_data[tri.indices.v1_id -1].z;
        temp2.x = scene.vertex_data[tri.indices.v0_id -1].x - scene.vertex_data[tri.indices.v1_id -1].x;
        temp2.y = scene.vertex_data[tri.indices.v0_id -1].y - scene.vertex_data[tri.indices.v1_id -1].y;
        temp2.z = scene.vertex_data[tri.indices.v0_id -1].z - scene.vertex_data[tri.indices.v1_id -1].z;
        temp1 = cross(temp1,temp2);
        tempf = sqrt(dot(temp1,temp1));
        if (tempf != 0) {
        temp1.x /= tempf;
        temp1.y /= tempf;
        temp1.z /= tempf;
        }
        triNormals.push_back(temp1);
    }
    
    for (auto mesh : scene.meshes) {
        std::vector<parser::Vec3f> temp;
        for (auto face : mesh.faces) {
            parser::Vec3f temp1, temp2;
            temp1.x = scene.vertex_data[face.v2_id -1].x - scene.vertex_data[face.v1_id -1].x;
            temp1.y = scene.vertex_data[face.v2_id -1].y - scene.vertex_data[face.v1_id -1].y;
            temp1.z = scene.vertex_data[face.v2_id -1].z - scene.vertex_data[face.v1_id -1].z;
            temp2.x = scene.vertex_data[face.v0_id -1].x - scene.vertex_data[face.v1_id -1].x;
            temp2.y = scene.vertex_data[face.v0_id -1].y - scene.vertex_data[face.v1_id -1].y;
            temp2.z = scene.vertex_data[face.v0_id -1].z - scene.vertex_data[face.v1_id -1].z;
            temp1 = cross(temp1,temp2);
            tempf = sqrt(dot(temp1,temp1));
            if (tempf != 0) {
            temp1.x /= tempf;
            temp1.y /= tempf;
            temp1.z /= tempf;
            }
            temp.push_back(temp1);
        }
        meshNormals.push_back(temp);
    }
    
    

    //loop to iterate over each camera
    for (auto cam : scene.cameras) {
    
    //cross product to find u vector
    parser::Vec3f u = cross(cam.gaze, cam.up);
    
    //prelim finding q for ray eq
    parser::Vec3f q;
    q.x = (cam.position.x + (cam.near_distance * cam.gaze.x)) + (cam.near_plane.x * u.x) + (cam.near_plane.w * cam.up.x);
    q.y = (cam.position.y + (cam.near_distance * cam.gaze.y)) + (cam.near_plane.x * u.y) + (cam.near_plane.w * cam.up.y);
    q.z = (cam.position.z + (cam.near_distance * cam.gaze.z)) + (cam.near_plane.x * u.z) + (cam.near_plane.w * cam.up.z);
    
    int tot_pixels = cam.image_width * cam.image_height;
    unsigned char* image = new unsigned char [tot_pixels * 3];
    int imageI = 0;
    
    //loop over each pixel
    for (int height = 0; height < cam.image_height; height++) {
    for (int width = 0; width < cam.image_width; width++) {
        
    //calculating ray eq scaling parameter and normalizing
    float s_u = (width + 0.5) * ((cam.near_plane.y - cam.near_plane.x)/cam.image_width);
    float s_v = (height + 0.5) * ((cam.near_plane.w - cam.near_plane.z)/cam.image_height);
    parser::Vec3f d;
    d.x = (q.x + (s_u*u.x) - (s_v*cam.up.x)) - cam.position.x;
    d.y = (q.y + (s_u*u.y) - (s_v*cam.up.y)) - cam.position.y;
    d.z = (q.z + (s_u*u.z) - (s_v*cam.up.z)) - cam.position.z;
    float dDot = sqrt(dot(d,d));
    d.x /= dDot;
    d.y /= dDot;
    d.z /= dDot;
    
    //calculating the color for the ray (ez)
    depth = 0;
    parser::Vec3f color = computeColor(scene, cam.position, d);
    if (color.x < 0) color.x = 0; if (color.x > 255) color.x = 255;
    if (color.y < 0) color.y = 0; if (color.y > 255) color.y = 255;
    if (color.z < 0) color.z = 0; if (color.z > 255) color.z = 255;
    
    if ((color.x - (int)color.x) > 0.5) color.x = (int)color.x + 1;
    if ((color.y - (int)color.y) > 0.5) color.y = (int)color.y + 1;
    if ((color.z - (int)color.z) > 0.5) color.z = (int)color.z + 1;
    image[imageI++] = color.x;
    image[imageI++] = color.y;
    image[imageI++] = color.z;
    }
    } // end loop over each pixel
    
    //writing the final image to ppm
    
    write_ppm(cam.image_name.c_str(), image, cam.image_width, cam.image_height);
    free(image);
    
    } // end loop to iterate over each camera
}





