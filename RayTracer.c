/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

#include "utils.h"
#define HEIGHT 1000
#define WIDTH 1000
#define max(A,B) ((A)<(B)?(B):(A))
// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;
int object_list_size = 0;



void buildScene(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Let's add a plane
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);	// Note the plane is highly-reflective (rs=rg=.75) so we
						// should see some reflections if all is done properly.
						// Colour is close to cyan, and currently the plane is
						// completely opaque (alpha=1). The refraction index is
						// meaningless since alpha=1
 Scale(o,6,6,1);				// Do a few transforms...
 RotateZ(o,PI/1.20);
 RotateX(o,PI/2.25);
 Translate(o,0,-3,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
						// and store the inverse
						// transform for this object!
 insertObject
 (o,&object_list);			// Insert into object list
 object_list_size += 1;

 // Let's add a couple spheres
 o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);
 Scale(o,.75,.5,1.5);
 RotateY(o,PI/2);
 Translate(o,-1.45,1.1,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 object_list_size += 1;

 o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
 Scale(o,.5,2.0,1.0);
 RotateZ(o,PI/1.5);
 Translate(o,1.75,1.25,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 object_list_size += 1;

 // Insert a single point light source.
 p.px=0;
 p.py=15.5;
 p.pz=-5.5;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
 // This function implements the shading model as described in lecture. It takes
 // - A pointer to the first object intersected by the ray (to get the colour properties)
 // - The coordinates of the intersection point (in world coordinates)
 // - The normal at the point
 // - The ray (needed to determine the reflection direction to use for the global component, as well as for
 //   the Phong specular component)
 // - The current racursion depth
 // - The (a,b) texture coordinates (meaningless unless texture is enabled)
 //
 // Returns:
 // - The colour for this ray (using the col pointer)
 //

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////

 struct point3D *Lv;
 struct point3D *Rv;
 struct point3D *Vv;
 Lv = (struct point3D*)malloc(sizeof(struct point3D));
 Rv = (struct point3D*)malloc(sizeof(struct point3D));
 Vv = (struct point3D*)malloc(sizeof(struct point3D));

 memcpy(&light_list->p0, Lv, sizeof(struct point3D)); // index through light list later!!!
 subVectors(p,Lv);
 Lv->px = (-1)*(Lv->px);
 Lv->py = (-1)*(Lv->py);
 Lv->pz = (-1)*(Lv->pz);
 Lv->pw = (-1)*(Lv->pw);
 normalize(Lv);
 normalize(n);

 double comp = dot(Lv,n);
 Rv->px = 2*comp*(n->px);
 Rv->py = 2*comp*(n->py);
 Rv->pz = 2*comp*(n->pz);
 Rv->pw = 2*comp*(n->pw);
 subVectors(Lv,Rv);
 Rv->px = (-1)*(Rv->px);
 Rv->py = (-1)*(Rv->py);
 Rv->pz = (-1)*(Rv->pz);
 Rv->pw = (-1)*(Rv->pw);
 normalize(Rv);
 Vv->px = (ray->p0.px)-(p->px); 
 Vv->py = (ray->p0.py)-(p->py);
 Vv->pz = (ray->p0.pz)-(p->pz);
 Vv->pw = (ray->p0.pw)-(p->pw);
 normalize(Vv);

 double ambR, ambG, ambB, difR, difG, difB, spcR, spcG, spcB;

 ambR = (obj->alb.ra)*(R);
 ambG = (obj->alb.ra)*(G);
 ambB = (obj->alb.ra)*(B);

 difR = (obj->alb.rd)*(max(0.0, comp))*(R);
 difG = (obj->alb.rd)*(max(0.0, comp))*(G);
 difB = (obj->alb.rd)*(max(0.0, comp))*(B);

 spcR = (obj->alb.rs)*(max(0.0, pow(dot(Vv,Rv),(obj->shinyness))))*(R);
 spcG = (obj->alb.rs)*(max(0.0, pow(dot(Vv,Rv),(obj->shinyness))))*(G);
 spcB = (obj->alb.rs)*(max(0.0, pow(dot(Vv,Rv),(obj->shinyness))))*(B);

 // Be sure to update 'col' with the final colour computed here!
 col->R = (obj->alb.rg)*(ambR + difR + spcR);
 col->G = (obj->alb.rg)*(ambG + difG + spcG);
 col->B = (obj->alb.rg)*(ambB + difB + spcB);
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
  double min_dist = 99999;
  double curr_len = 0;
  double lam;
  struct point3D *tp = (struct point3D *) malloc(sizeof(struct point3D));
  struct point3D *tn = (struct point3D *) malloc(sizeof(struct point3D));
  struct object3D *curr_obj;
  curr_obj = object_list;
  *lambda = 0;	
  	while(curr_obj != NULL){ //comment out loop, just get just plane to render
  		if(curr_obj!=Os){ 
			curr_obj->intersect(curr_obj, ray, &lam, tp, tn, a, b);
			if(lam>0 && obj!=NULL){
				curr_len = lam;
				if ((min_dist > curr_len) && (curr_len > 0)){
					min_dist = curr_len;
					*obj = curr_obj;
					*lambda = lam;
					memcpy(p, tp, sizeof(struct point3D));
					memcpy(n, tn, sizeof(struct point3D));
				}
			}else if(lam>0 && obj==NULL{
				*lambda = lam;
			}
	 	}
	 	curr_obj = curr_obj->next;
	}

}


void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //

	 double lambda, shadowlam;		// Lambda at intersection
	 double a,b;		// Texture coordinates
	 struct object3D *obj;
	 struct object3D *obj2;	// Pointer to object at intersection
	 struct point3D p, tp;	// Intersection point
	 struct point3D n, tn;	// Normal at intersection
	 struct colourRGB I;	// Colour returned by shading function
	
	 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
	 {
	  col->R=-1;
	  col->G=-1;
	  col->B=-1;
	  return;
	 } 
	 else
	 {
		findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
				
//		if (obj != NULL){
		if (lambda > 0){ //intersection == true, color it

			struct ray3D* lightray;
			struct point3D *dir = (struct point3D*)malloc(sizeof(struct point3D));
		 	struct pointLS *curr_light = light_list;
		 	//while (curr_light!=NULL)
		 	//{
				memcpy(dir, &(curr_light->p0), sizeof(struct point3D));
				subVectors(&p, dir);
				dir->pw = 0;
				lightray = newRay(&p, dir);
				normalize(&(lightray->d));
				findFirstHit(lightray, &shadowlam, obj, NULL, &p, &n, &a, &b);
				//printf("p %f %f %f \n", p.px,p.py,p.pz);
				if(shadowlam <=0){
					//printf("intersected!!!");
					//rtShade(obj, &p, &n, ray, depth, a, b, &I);
					
					col->R = obj->col.R;
					col->G = obj->col.G;
					col->B = obj->col.B;
					//printf("RGB %f, %f, %f\n", col->R, col->G, col->B);
					
				}else{ //don't color the shadow
					col->R = 0;
					col->G = 0;
					col->B = 0;
					printf("it never gets here\n");
				}
				//curr_light = curr_light->next;
			//}

			/*if(depth>0){
				//create reflected ray
				rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
			}*/
			
		
		}
		else //darkness
		{
			col->R =0.0;
			col->G =0.0;
			col->B =0.0;
		}


	 }	
 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////
}

int main(int argc, char *argv[])
{
	 // Main function for the raytracer. Parses input parameters,
	 // sets up the initial blank image, and calls the functions
	 // that set up the scene and do the raytracing.
	 struct image *im;	// Will hold the raytraced image
	 struct view *cam;	// Camera and view for this scene
	 int sx;		// Size of the raytraced image
	 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
	 char output_name[1024];	// Name of the output file for the raytraced .ppm image
	 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
	 struct point3D g;
	 struct point3D up;
	 double du, dv;			// Increase along u and v directions for pixel coordinates
	 struct point3D* pc;
	 struct point3D* d;		// Point structures to keep the coordinates of a pixel and
					// the direction or a ray
	 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
	 struct colourRGB col;		// Return colour for raytraced pixels
	 struct colourRGB background;   // Background colour
	 int i,j;			// Counters for pixel coordinates
	 unsigned char *rgbIm;
		int total = 0, total2 = 0;
	 if (argc<5)
	 {
	  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
	  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
	  fprintf(stderr,"   size = Image size (both along x and y)\n");
	  fprintf(stderr,"   rec_depth = Recursion depth\n");
	  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
	  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
	  exit(0);
	 }
	 sx=atoi(argv[1]);
	 MAX_DEPTH=atoi(argv[2]);
	 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
	 strcpy(&output_name[0],argv[4]);

	 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
	 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
	 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
	 else fprintf(stderr,"Antialising is on\n");
	 fprintf(stderr,"Output file name: %s\n",output_name);

	 object_list=NULL;
	 light_list=NULL;

	 // Allocate memory for the new image
	 im=newImage(sx, sx);
	 if (!im)
	 {
	  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
	  exit(0);
	 }
	 else rgbIm=(unsigned char *)im->rgbdata;

	 ///////////////////////////////////////////////////
	 // TO DO: You will need to implement several of the
	 //        functions below. For Assignment 3, you can use
	 //        the simple scene already provided. But
	 //        for Assignment 4 you need to create your own
	 //        *interesting* scene.
	 ///////////////////////////////////////////////////
	 buildScene();		// Create a scene. This defines all the
				// objects in the world of the raytracer

	 //////////////////////////////////////////
	 // TO DO: For Assignment 3 you can use the setup
	 //        already provided here. For Assignment 4
	 //        you may want to move the camera
	 //        and change the view parameters
	 //        to suit your scene.
	 //////////////////////////////////////////

	 // Mind the homogeneous coordinate w of all vectors below. DO NOT
	 // forget to set it to 1, or you'll get junk out of the
	 // geometric transformations later on.

	 // Camera center is at (0,0,-1)
	 e.px=0;
	 e.py=0;
	 e.pz=-3;
	 e.pw=1;

	 // To define the gaze vector, we choose a point 'pc' in the scene that
	 // the camera is looking at, and do the vector subtraction pc-e.
	 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
	 g.px=0;
	 g.py=0;
	 g.pz=1;
	 g.pw=0;

	 // Define the 'up' vector to be the Y axis
	 up.px=0;
	 up.py=1;
	 up.pz=0;
	 up.pw=0;

	 // Set up view with given the above vectors, a 4x4 window,
	 // and a focal length of -1 (why? where is the image plane?)
	 // Note that the top-left corner of the window is at (-2, 2)
	 // in camera coordinates.
	 cam=setupView(&e, &g, &up, -3, -1, 1, 2);

	 if (cam==NULL)
	 {
	  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
	  cleanup(object_list,light_list);
	  deleteImage(im);
	  exit(0);
	 }

	 // Set up background colour here
	 background.R=0;
	 background.G=0;
	 background.B=0;

	 // Do the raytracing
	 //////////////////////////////////////////////////////
	 // TO DO: You will need code here to do the raytracing
	 //        for each pixel in the image. Refer to the
	 //        lecture notes, in particular, to the
	 //        raytracing pseudocode, for details on what
	 //        to do here. Make sure you undersand the
	 //        overall procedure of raytracing for a single
	 //        pixel.
	 //////////////////////////////////////////////////////
	 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
	 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
					// and dv is negative since y increases downward in pixel
					// coordinates and upward in camera coordinates.

	 fprintf(stderr,"View parameters:\n");
	 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
	 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
	 printmatrix(cam->C2W);
	 fprintf(stderr,"World to camera conversion matrix\n");
	 printmatrix(cam->W2C);
	 fprintf(stderr,"\n");
	 #pragma omp parallel for
	 //fprintf(stderr,"Rendering row: ");
	 for (j=0;j<sx;j++)		// For each of the pixels in the image
	 {
	 // fprintf(stderr,"%d/%d, ",j,sx);
	  for (i=0;i<sx;i++)
	  {
		///////////////////////////////////////////////////////////////////
		// TO DO - complete the code that should be in this loop to do the
		//         raytracing!
		///////////////////////////////////////////////////////////////////
		
		//from 418notes =: (wl+(i+0.5)Δu,wt+(j+0.5)Δv,f)
		d = (struct point3D*)malloc(sizeof(point3D)); //dir vector
		pc = (struct point3D*)malloc(sizeof(point3D)); //image plane
		struct point3D *origin = newPoint(0,0,0); //camera origin
		//image plane
		pc->px = cam->wl+ (i+0.5)*du;
		pc->py = cam->wt + (j+0.5)*dv;
		pc->pz = -1;
		pc->pw = 1;

		//bring image plane and origin to world
		matVecMult(cam->C2W, pc);
		matVecMult(cam->C2W, origin);
		//create ray in world space
		memcpy(d, pc, sizeof(struct point3D));
		subVectors(origin, d);
		d->pw = 0;
		ray = newRay(origin, d);
				
		if(length(&(ray->d))<0) //not supposed to happen
			printf("%f,%f, %f\n", ray->d.px,ray->d.py, ray->d.pz);
		
		rayTrace(ray, MAX_DEPTH, &col, NULL);

		//coloring the pixel with u and v values as per davids suggestion for debugging
		//if(col.R!=0 && col.G!=0 && col.B!=0){
			((unsigned char*)im->rgbdata)[(j*sx + i)*3]   = (unsigned char) min(col.R*255, 255);
			((unsigned char*)im->rgbdata)[(j*sx + i)*3+1] = (unsigned char) min(col.G*255, 255);
			((unsigned char*)im->rgbdata)[(j*sx + i)*3+2] = (unsigned char) min(col.B*255, 255);
		//}

	  } // end for i
	 } // end for j
	printf("pixel %d pixel 2 %d\n",total, total2);
	 fprintf(stderr,"\nDone!\n");
	 //printmatrix(object_list[0].T);
	 //printmatrix(object_list[0].Tinv);

	 // Output rendered image
	 imageOutput(im,output_name);

	 // Exit section. Clean up and return.
	  cleanup(object_list,light_list);		// Object and light lists
	  deleteImage(im);				// Rendered image
	 free(cam);					// camera view
	 exit(0);
}
