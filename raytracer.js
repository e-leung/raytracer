let tmin = 0.001;
let tmax = 1e20;

function RayTracer( canvasID ) {
  // setup the canvas
  this.canvas = document.getElementById(canvasID);

  // setup the background style: current options are 'daylight' or 'white'
  this.sky = 'daylight';

  // initialize the objects and lights
  this.objects = new Array();
  this.lights  = new Array();
}

function Camera( eye , center , up , fov , aspect ){
  this.eye          = eye;
  this.center       = center;
  this.up           = up;
  this.fov          = fov;
  this.aspect       = aspect;
  this.orthonormalB = mat3.create();

  // calculate gaze = lookat - eye
  let gaze = vec3.create();
  vec3.subtract( gaze , this.center , this.eye );

  // w = -gaze / ||gaze||
  let w = vec3.create();
  vec3.normalize( w , gaze );
  vec3.scale( w , w , -1.0 );

 // compute u = up x w
  let u = vec3.create();
  vec3.cross(u, this.up, w);
  vec3.normalize(u, u);

  // v = w x u
  let v = vec3.create();
  vec3.cross( v , w , u );

  //  B = b0  b3  b6  
  //      b1  b4  b7   
  //      b2  b5  b8  
  for (let d = 0; d < 3; d++) {
    this.orthonormalB[d]    = u[d];
    this.orthonormalB[3+d]  = v[d];
    this.orthonormalB[6+d]  = w[d];
  }
}

Camera.prototype.makeRay = function ( px , py ){
  let d       = vec3.distance( this.center, this.eye );
  let height  = 2.0 * d * Math.tan( this.fov/2.0 );
  let width   = this.aspect * height;

  let q = vec3.fromValues(
    -width/2.0  + width  * px,
    -height/2.0 + height * py,
    -d
  );

  let direction = vec3.create();

  vec3.transformMat3( direction, q, this.orthonormalB );
  vec3.normalize(direction, direction);

  let result = new Ray(this.eye, direction);
  return result;
}

Sphere.prototype.intersect = function( ray , tmin , tmax ) {
  const r  = ray.direction; // ray direction
  vec3.normalize(r, r);

  const x0 = ray.origin;    // ray origin
  const R  = this.radius;   // sphere radius
  const c  = this.center;   // sphere center

  let o = vec3.create();
  vec3.subtract( o , x0 , c ); // compute o = x0 - c
  let B = vec3.dot(r,o);       // B = dot product of r and o
  let C = vec3.dot(o,o) - R*R; // C = ||x0 - c||^2 - R^2

  let discriminant = B*B - C;
  if (discriminant < 0) {
    return 1e10;
  }
  let t1 = -B - Math.sqrt(discriminant);
  let t2 = -B + Math.sqrt(discriminant);

  let minimum = Math.min(t1, t2);
  if ((minimum < tmax) && (minimum > tmin)) {
   return minimum;
  }  
}

RayTracer.prototype.hit = function( ray , tmin , tmax ){
  let mySpheres = [];

  for (let k = 0; k < this.objects.length; k++) {
    let hit = this.objects[k].intersect(ray, tmin, tmax);
    mySpheres.push(hit);
  }

  let lowestT = 1e10;
  let lowestI = 1e10;

  //find the index of the lowest t value
  for (let i = 0; i < mySpheres.length; i++) {
    if (mySpheres[i] < lowestT) {
      lowestT = mySpheres[i];
      lowestI = i;
    }
  }
  
  if (lowestI == 1e10) {
    return undefined;
  }
  
  else {
    // hit_info = [ object , intersection , normal , t ]
    let intersection = vec3.create();
    //x(t) = x0 + rt
    vec3.scaleAndAdd(intersection, ray.origin, ray.direction, lowestT); 
    let normal = vec3.create();
    vec3.subtract(normal, intersection, this.objects[lowestI].center);
    vec3.normalize(normal, normal);
    
    return new Hit_info (this.objects[lowestI] , intersection , normal , lowestT);
  }
}

RayTracer.prototype.color = function( ray , depth ) {
  //base case
  if(depth == 0){
    return this.background(ray);
  }
  
  let myHit_info = this.hit( ray  , 0.001 , 1e20 );

  // use the background color if there is no intersection
  if (myHit_info == undefined){
    return this.background(ray);
  }

  if(this.lights.length == 0){
    return myHit_info.object.color;
  }
  
  // initialize color from ambient light
  let l_a = vec3.fromValues(1.0, 1.0, 1.0);

  let k_a = vec3.create();
  vec3.scale(k_a, myHit_info.object.material.color, 0.4)

  let ambient = vec3.create()
  vec3.multiply(ambient, k_a, l_a);

  let color = vec3.create()
  vec3.add(color, color, ambient);  

  // compute the diffuse and specular components for each light in the scene
  for (let i = 0; i < this.lights.length; i++){

    let light = this.lights[i];
    // determine if we are in the shadow of another object
    let direction_to_light = vec3.create();
    vec3.subtract(direction_to_light, light.location, myHit_info.intersection);
    
    let shadow_ray = new Ray( myHit_info.intersection , direction_to_light );

    let blocking_object = this.hit( shadow_ray , 0.0001 , 1e20 );
    if (blocking_object){
      continue; 
    }

    let phong = myHit_info.object.material.shade( ray , light , myHit_info );
    vec3.scale(phong, phong, 0.7);
    vec3.add(color, color, phong);
  }

  let scattered_ray = myHit_info.object.material.scatter( ray , myHit_info );  
  if (scattered_ray == undefined) {
    return color;
  }
  let color_scatter = this.color( scattered_ray , depth-1 );
  if (myHit_info.object.material.type == "refractive"){
    return color_scatter;
  }
  vec3.scale(color_scatter, color_scatter, 0.3);
  vec3.add( color , color , color_scatter );

  return color;
}

Material.prototype.scatter = function ( ray , hit_info ){
  
  if (this.type == "refractive"){
    return this.refract(ray, hit_info);
  }
  else if (this.type == "reflective"){
    return this.reflect(ray, hit_info);
  } 
  else {
    return undefined;
  }
}

Material.prototype.reflect = function( ray , hit_info ){
  let r = vec3.create();
  vec3.scale( r , hit_info.normal , -2.0 * vec3.dot( ray.direction , hit_info.normal ) );

  vec3.add( r , r , ray.direction );
  return new Ray( hit_info.intersection , r );
}

Material.prototype.refract = function( ray , hit_info ){
  let eta = 1.5;
  let n = hit_info.normal.slice();
  let eta1_over_eta2 = 1.5;

  let dt = vec3.dot(ray.direction, n);
  if (dt > 0){
    vec3.scale(n, n, -1);
    //eta1_over_eta2 = eta/1.0;
  }

  else {
    eta1_over_eta2 = 1.0/eta1_over_eta2;
  }

  dt = vec3.dot(ray.direction, n);
  
  let discriminant = 1.0 - eta1_over_eta2 * eta1_over_eta2 * ( 1.0 - dt*dt );
  if (discriminant < 0.0) {
    return this.reflect(ray, hit_info);
  }

  let r1 = vec3.create();
  vec3.scaleAndAdd( r1 , ray.direction , n , -dt );
  vec3.scale( r1 , r1 , eta1_over_eta2 );

  let r2 = vec3.create();
  vec3.scale( r2 , n , Math.sqrt(discriminant) );

  let r = vec3.create();
  vec3.subtract(r,r1,r2);

  return new Ray (hit_info.intersection , r);
}

Material.prototype.shade = function ( ray , light , hit_info ){
  const x = hit_info.intersection; // intersection point
  const n = hit_info.normal;       // outward normal of surface
  const y = light.location;        // light location
  const L = light.color; 
  let kd = hit_info.object.material.k_d; 
  let ks = hit_info.object.material.k_s; 
  // compute the vector from the surface to the light
  let l = vec3.create();
  vec3.subtract( l , y , x );
  vec3.normalize(l , l);

  vec3.scale(ray.direction, ray.direction, -1);
  let h = vec3.create();
  vec3.add(h, ray.direction, l);
  vec3.normalize(h, h); 
  vec3.scale(ray.direction, ray.direction, -1);
  // fraction of light reflected
  let cos_theta = Math.max(0.0,vec3.dot(n,l));

  // evaluate the diffuse term
  let cd = vec3.create();
  vec3.multiply( cd , kd , L );
  vec3.scale( cd , cd , cos_theta );

  if(this.shine == undefined){
    return cd;
  }
  
  let power = Math.max(0.0 , Math.pow(vec3.dot(h, n),this.shine));
  let cs = vec3.create();
  vec3.multiply( cs , ks , L );
  vec3.scale( cs , cs , power );
  let sum = vec3.create();
  vec3.add(sum, cd, cs);

  return sum;
}

RayTracer.prototype.draw = function(eye, center, up, fov, aspect, alias) {
  // get the canvas and the image data we will write to
  let context = this.canvas.getContext('2d');
  let image = context.createImageData(this.canvas.width,this.canvas.height);

  // numbers of pixels in x- and y- directions
  const nx = image.width;
  const ny = image.height;
  let myCam = new Camera(eye, center, up, fov, aspect);
  
  // loop through the canvas pixels
  for (let j = 0; j < ny; j++) {   //ny
    for (let i = 0; i < nx; i++) { //nx

      // if alias is true, find the average of 20 subpixel samples
      if (alias == true){
        let accumColor = vec3.create();
        for (let k = 0; k < 21; k++){
          let px = (i + Math.random()) / nx;       
          let py = (ny - j - Math.random()) / ny;
          let color = vec3.create();
          let ray = myCam.makeRay( px , py );

          color = this.color( ray , 5 );
          vec3.add(accumColor, color, accumColor);
        }
        vec3.scale(accumColor, accumColor, 1/20);
        this.setPixel( image , i , j , accumColor[0] , accumColor[1] , accumColor[2] );
      }
      else {
        // compute pixel coordinates in [0,1] x [0,1]
        let px = (i + 0.5) / nx;      // sample at pixel center
        let py = (ny - j - 0.5) / ny; // canvas has y pointing down, but image plane has y going up

        let color = vec3.create();

        let ray = myCam.makeRay( px , py );
        color = this.color( ray , 5 );
        this.setPixel( image , i , j , color[0] , color[1] , color[2] );
      }
    }
  }
  context.putImageData(image,0,0);
}

RayTracer.prototype.background = function( ray ) {
  if (this.sky === 'white') {
    // a white sky
    return vec3.fromValues(1,1,1);
  }
  else if (this.sky === 'daylight') {
    // a light blue sky :)
    let t = 0.5*ray.direction[1] + 0.2; // uses the y-values of ray.direction
    if (ray.direction == undefined) t = 0.2; // remove this if you have a different name for ray.direction
    let color = vec3.create();
    vec3.lerp( color , vec3.fromValues(.5,.7,1.)  , vec3.fromValues(1,1,1) , t );
    return color;
  }
  else
    alert('unknown sky ',this.sky);
}

RayTracer.prototype.setPixel = function( image , x , y , r , g , b ) {
  let offset = (image.width * y + x) * 4;
  image.data[offset  ] = 255*Math.min(r,1.0);
  image.data[offset+1] = 255*Math.min(g,1.0);
  image.data[offset+2] = 255*Math.min(b,1.0);
  image.data[offset+3] = 255; // alpha: transparent [0-255] opaque
}

function Ray ( origin , direction ){
  this.origin    = origin;
  this.direction = direction;
}

function Hit_info( object , intersection , normal , t ){
  this.object       = object;
  this.intersection = intersection;
  this.normal       = normal;
  this.t            = t;
}

function Sphere(params) {
  // represents a sphere object
  this.center   = params['center'];   // center of the sphere (vec3)
  this.radius   = params['radius'];   // radius of sphere (float)
  this.material = params['material']; // material used to shade the sphere (see 'Material' below)
  this.name     = params['name'] || 'sphere'; // a name to identify the sphere (useful for debugging) (string)
  this.color    = params['color'];    // color of the sphere
}

function Light(params) {
  // describes a point light source, storing the location of the light
  // as well as ambient, diffuse and specular components of the light
  this.location = params.location; // location of 
  this.color    = params.color || vec3.fromValues(1,1,1); // default to white (vec3)
  this.l_a       = vec3.fromValues(1.0, 1.0, 1.0);
  this.l_s       = vec3.fromValues(1.0, 1.0, 1.0);
}

function Material( params ) {
  // represents a generic material class
  this.type  = params.type; // diffuse, reflective, refractive (string)
  this.shine = params.shine; // phong exponent (float)
  this.color = params.color || vec3.fromValues(0.5,0.5,0.5); // default to gray color (vec3)
  let k_a     = vec3.create();
  vec3.scale(k_a, this.color, 0.4);
  this.k_a    = k_a;
  this.k_d    = this.color;
  this.k_s    = vec3.fromValues(1, 1, 1);
}