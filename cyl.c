void cylinderIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
  
  struct ray3D *transformed_ray=(struct ray3D *)calloc(1,sizeof(struct ray3D));
  rayTransform(ray,transformed_ray,cylinder);
  double a1 = transformed_ray->d.px * transformed_ray->d.px + transformed_ray->d.pz * transformed_ray->d.pz;
  double b1 = 2 * transformed_ray->p0.px * transformed_ray->d.px + 2 * transformed_ray -> p0.pz * transformed_ray -> d.pz;
  double c1 = transformed_ray->p0.px * transformed_ray->p0.px + transformed_ray->p0.pz * transformed_ray->p0.pz - 1;
  double root = b1*b1 - 4*a1*c1;
  *a = 20;
  if (root < 0)
  {
    free(transformed_ray);
    *lambda = -1;
   // printf("doesnt \n");
    return;
  }
 
    float sqb24ac = sqrt(root);
    float t0 = (-b1 + sqb24ac)/(2*a1);
    float t1 = (-b1 - sqb24ac)/(2*a1);
    if (t0 > t1){float tmp = t0; t0 = t1; t1 = tmp;}
    float y0 = transformed_ray->p0.py + t0 * transformed_ray->d.py;
    float y1 = transformed_ray->p0.py + t1 * transformed_ray->d.py;
    
    if (y0 < -1)
    {
     //printf("DAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \n");
      if (y1 < -1)
      {
	*lambda = -1;
	free(transformed_ray);
	//printf("doesnt \n");
	return;
      }
      else
      {
	float th = t0 + (t1 - t0)*(y0 + 1)/(y0 - y1);
	if (th <= 0)
	{
	  free(transformed_ray);
	  *lambda = -1;
	  return;
	}
        struct point3D *normal = newPoint(0, -1, 0);
        normal->pw = 0;
	p->px = transformed_ray->p0.px + th*transformed_ray->d.px;
	p->py = transformed_ray->p0.py + th*transformed_ray->d.py;
	p->pz = transformed_ray->p0.pz + th*transformed_ray->d.pz;
	p ->pw = 1;
	*lambda = th;
        normalTransform(normal,n,cylinder);
        free(normal);

        matVecMult(cylinder-> T, p); //Do you need this
	free(transformed_ray);
	return;
	
	
	
	//ovde dodaj rezultat
      }
    }
    else if ((y0>= -1)&&(y0 <= 1))
    {
     
      if (t0 <= 0)
	{
	  free(transformed_ray);
	  *lambda = -1;
	  return;
	}
	
	p->px = transformed_ray->p0.px + t0*transformed_ray->d.px;
	p->py = transformed_ray->p0.py + t0*transformed_ray->d.py;
	p->pz = transformed_ray->p0.pz + t0*transformed_ray->d.pz;
	p ->pw = 1;
   *lambda = t0;
        struct point3D *normal = newPoint(p->px, 0, p->pz);
        normal->pw = 0;
        normalTransform(normal,n,cylinder);
        free(normal);

        matVecMult(cylinder-> T, p); //Do you need this
	free(transformed_ray);
	return;
    }
    else if (y0 > 1)
    {
     
      if (y1 > 1)
      {
	free(transformed_ray);
	*lambda = -1;
	return;
      }
      else
      {
	float th = t0 + (t1 - t0)*(y0 -1)/(y0 - y1);
	if (th <= 0) 
	{
	  free(transformed_ray);
	  *lambda = -1;
	  return;
	}
	p->px = transformed_ray->p0.px + th*transformed_ray->d.px;
	p->py = transformed_ray->p0.py + th*transformed_ray->d.py;
	p->pz = transformed_ray->p0.pz + th*transformed_ray->d.pz;
	p ->pw = 1;
	*lambda = th;
        struct point3D *normal = newPoint(0, 1, 0);
        normal->pw = 0;
        normalTransform(normal,n,cylinder);
        free(normal);
        matVecMult(cylinder-> T, p); //Do you need this
	free(transformed_ray);
	return;
	
      }
      
    }
    // a ako nista return false
    free(transformed_ray);  
    *lambda = -1;
    return;
}
Chat Conversation End
