# Hello Scalismo!

The goal in this tutorial is to present the most important data structures, as well as the visualization capabilities of Scalismo.

## Initializing the system

Before we start, we need to initialize Scalismo by calling:

```scala mdoc
scalismo.initialize()
implicit val rng = scalismo.utils.Random(42)
```

The call to ```scalismo.initialize``` loads all the dependencies to native C++ libraries (such as e.g. [vtk](https://www.vtk.org) or [hdf5](https://www.hdf-group.org)). The second call, tells scalismo, which source
of randomness it should be using and at the same time seeds the random number generator appropriately.

Later on we would like to visualize the objects we create. This is done using [Scalismo-ui](https://github.com/unibas-gravis/scalismo-ui) - the visualization library accompanying scalismo. 
We can load an instance of the GUI, which we name here simply ```ui``` as follows:

```scala mdoc
import scalismo.ui.api.ScalismoUI

val ui = ScalismoUI()
```


## Meshes (surface data)

The first fundamental data structure we discuss is the Triangle Mesh.

Meshes can be read from a File using the method ```readMesh``` from the ```MeshIO```:

```scala mdoc
import scalismo.io.MeshIO

val mesh = MeshIO.readMesh(new java.io.File("datasets/Paola.stl")).get
``` 

To visualize any object in Scalismo, we can use the ```show``` method of the ```ui``` object. We will see later 
on that we often want to organize different visualization of an object in a group. We start directly with this practice and 
first create a new group, to which we then add a visualization:

```scala mdoc
val paolaGroup = ui.createGroup("paola")
val meshView = ui.show(paolaGroup, mesh, "Paola")
```

Now that the mesh is displayed in the "Scalismo Viewer's 3D view", you can interact with it as follows: 

* to rotate: maintain the left mouse button clicked and drag  
* to shift/translate: maintain the middle mouse button clicked and drag
* to scale: maintain the right mouse button clicked and drag up or down
 
###### Note: if you are a Mac user, please find out how to emulate these events using your mouse or trackpad

Note also that you can use the *RC*, *X*, *Y* and *Z* buttons in the 3D view to recenter the camera on the displayed object.

#### Anatomy of a Triangle mesh
A 3D triangle mesh in scalismo consists of a ```pointSet```, which maintains a collection of 3D points and a list of triangle cells. We can access individual points using their PointId. 
Here we show how we can access the first point in the mesh:

```scala mdoc
import scalismo.common.PointId

println("first point " + mesh.pointSet.point(PointId(0)))
``` 

Similarly, we can access the first triangles as follows:

```scala mdoc
import scalismo.mesh.TriangleId;

println("first cell " + mesh.triangulation.triangle(TriangleId(0)))
``` 

The first cell is a triangle between the first, second and third points of the mesh.
Notice here that the cell indicates the identifiers of the points (their index in the point sequence)
instead of the geometric position of the points.

Instead of visualizing the mesh, we can also display the points forming the mesh. 

```scala mdoc
ui.show(paolaGroup, mesh.pointSet, "pointCloud")
``` 

This should add a new point cloud element to the scene with the name "pointCloud".

###### Note: depending on your computer, visualizing the full point cloud may slow down the visualization performance.

Note that to clean up the 3D scene, you can delete the objects either from the user interface (by right-clicking on the object's name), or programmatically as such :

```tut:silent
remove("pointCloud")
``` 

## Points and Vectors

Different operations exist on 3D points in Scalismo: 

The difference between two points yields a vector:
```tut:silent
val v : Vector[_3D] = Point(4f,5f,6f) - Point(1f,2f,3f) 
```
The sum of a point with a vector yields a new point:
```tut:silent
val p : Point[_3D] = mesh.point(PointId(0)) + v 
```
Points can be converted to vectors:
```tut:silent
val v2 : Vector[_3D]= p.toVector
```
and vice versa:
```tut:silent
val p2 : Point[_3D] = v.toPoint 
```

##### Exercise: Given a list of points, return their center (solution does not have to be a one-liner). Hint: this could be a good occasion to try out the *map* and *reduce* functions presented in the cheat sheet.

```tut:silent
val pointList : List[Point[_3D]] = List(Point(4f,5f,6f), Point(1f,2f,3f), Point(14f,15f,16f), Point(7f,8f,9f), Point(
10f,11f,12f))
```

```tut:silent:fail
val center : Point[_3D] = ???
```

```tut:book
val vectors = pointList.map{p : Point[_3D] => p.toVector}  // use map to turn points into vectors
val vectorSum = vectors.reduce{ (v1,v2) => v1+v2} // reduce the collection of vectors by summing it
val centerV: Vector[_3D] = vectorSum * (1f / pointList.length ) // divide the sum by the number of points  
val center = centerV.toPoint
```

## Scalar Images

A discrete scalar image (e.g. gray level image) in Scalismo is simply a function from a discrete domain of points to a scalar value. 

Let's read and display a 3D image (MRI of a human):

```tut:silent
val image = ImageIO.read3DScalarImage[Short](new File("datasets/PaolaMRI.vtk")).get
show(image, "mri")
```
###### Note: depending on your view on the scene, it could appear as if the image is not displayed. In this case, make sure to rotate the scene and change the position of the slices as indicated below.

To visualize the different image slices in the viewer, select "Scene" (the upper node in the scene tree graph) and use the X,Y,Z sliders.

You can also change the way of visualizing the 3D scene under the

*View -> Perspective* menu.

### Scalar Image domain

Let's inspect the domain of the image :

```tut:silent
val origin : Point[_3D] = image.domain.origin
val spacing : Vector[_3D] = image.domain.spacing
val size : IntVector[_3D] = image.domain.size  
```

The discrete image domain is a 3-dimensional regular grid of points originating at point (92.5485, -121.926, 135.267), 
with regular spacing of 1.5 mm in each dimension and containing 171, 171, 139 grid slots in the x, y and z directions respectively. 

To better see this, let's display the first 172 points of the image domain 

```tut:silent
val imagePoints : Iterator[Point[_3D]] = image.domain.points.take(172)
show(imagePoints.toSeq, "imagePoints")
```

### Scalar image values

The other important part of a discrete image are the values associated with the domain points

```tut:silent
val values : Iterator[Short] = image.values
```
This is an iterator of scalar values of type Short as encoded in the read image.

Let's check the first value, that is the value associated with the origin : 

```tut:silent
image.values.next 
```

Given that the origin point has index (0,0,0), this same value can be obtained by accessing the image at this index :
```tut:silent
image(IntVector(0,0,0)) 
```

Naturally, the number of scalar values should be equal to the number of point domains 

```tut:silent
image.values.size == image.domain.points.size
```

Notice that you can check the intensity value at a particular point position in the image, by maintaining the Ctrl key pressed and hovering over the image. The intensity value will then be displayed in the lower left corner of the Scalismo viewer window.

##### Creating scalar images

Given that discrete scalar images are simply a mapping between points and values, we can very easily create such images programmatically.

Here we create a new image defined on the same domain of points with artificially created values:  we have a value of 1000 at every 50 th point, otherwise 0 : 

```tut:silent
val values = (0 until image.domain.numberOfPoints).map{ i:Int => if (i % 50 == 0) 1000 else 0}
val image2 = DiscreteScalarImage(image.domain, ScalarArray(values.toArray))
show(image2,"pattern")
```

##### Exercise: create and display a thresholded version of the MRI image, where all intensity values above 300 are replaced with 0 (make sure to use 0.toShort to have the correct type)

```tut:silent:fail
val threshValues : Array[Short] = image.values.map{v : Short => ???}.toArray
 
val thresholded : DiscreteScalarImage[_3D,Short] = ???  
```

```tut:book
val threshValues : Array[Short] = image.values.map{v :Short => if(v<= 300) v else 0.toShort}.toArray
val thresholded : DiscreteScalarImage[_3D,Short] = DiscreteScalarImage[_3D, Short](image.domain, ScalarArray(threshValues))
show(thresholded, "thresh")
```

Note that the same result can be achieved much easier, by directly *mapping* the thresholding function to the image as such :

```tut:silent
val thresholded2 : DiscreteScalarImage[_3D,Short] = image.map{v :Short => if(v<=300) v else 0.toShort}
```


## Statistical Mesh Models

Now let's load and show a statistical shape model of faces.

```tut:silent
val faceModel = StatismoIO.readStatismoMeshModel(new File("datasets/bfm.h5")).get
show(faceModel, "faceModel")
```
(suggestion: previously created objects are not needed in the following and can be removed from the scene.
You can do this by right-clicking on the *Static Objects* element and selecting "Remove all"
)

### Sampling in the UI
##### Exercise: Sample random instances of faces by using the graphical tools in the scene pane : click on the "Instance 1" tree node and then the "Random" button


##### Exercise: click a landmark on a position of the face model, e.g. chin or eye corner.. (use the toggle button "LM" in the toolbar to activate landmark clicking). Now continue sampling from the model. What happens to the selected point?

###### Note: in case landmark clicking did not work for you, you can also add landmarks by hovering over the position and pressing the *X* key, once the LM toggle button activated. Otherwise, you can add a landmark at the tip of the nose programmatically by executing:  

```tut:silent
addLandmarksTo(Seq(Landmark("A", faceModel.mean.point(PointId(8156)))), "faceModel")
```


As you can see, a new instance of the face model is displayed each time along with the corresponding landmark point. Notice how the position of the landmark point changes in space while it keeps the same "meaning" on the face (eye corner, tip of nose ..)

##### Retrieving Landmarks from the scene 

Once you clicked a Landmark point on an object in the scene, the landmark information can be retrieved using the getLandmarksOf method : 

```tut:book
addLandmarksTo(IndexedSeq(Landmark("B", faceModel.mean.point(PointId(86)))),"faceModel")
```

```tut:silent
val landmarks : Seq[Landmark[_3D]] = getLandmarksOf("faceModel").get
val landmarkId : String = landmarks(0).id
val landmarkPosition : Point[_3D] = landmarks(0).point
```

### Sampling programmatically

Let's now sample from the model programmatically. Execute the code below several times (remember that you can select it and press *Shift+Enter* in the code pane)

```tut:silent
remove("sampledFace")
val sampledFace : TriangleMesh = faceModel.sample
show(sampledFace, "sampledFace")
```

```tut:book
addLandmarksTo(IndexedSeq(Landmark("B", landmarkPosition + Vector(2f,2f,2f))),"faceModel")
```

##### Exercise: sample a 100 faces from the model and trace the position of point with id 610 in each sample. Show all the positions taken by this point as a point cloud. What does the point distribution look like? How could one possibly represent it? 

```tut:silent:fail
val pc610 :IndexedSeq[Point[_3D]] = ???
```

```tut:book
val pc610 :IndexedSeq[Point[_3D]] = (0 until 100).map(i => faceModel.sample.point(PointId(610)))
```

```tut:silent
show(pc610, "pc610")
```


##### Exercise: compute the mean of the point positions pc610 above and compare it to the position of the same point id on the mean mesh of the face model (either by displaying both, or computing the distance between the two). What do you notice?


Finally, let's sample several faces from the  model and compare the number of points in each sample

```tut:silent
(0 until 10).foreach{ i :Int => println(faceModel.sample.numberOfPoints) }
```

##### Why is the number of points always the same? 
We will find out the answer in the coming week.

<br> <br>
```tut:invisible
gui.close()
```

```scala mdoc
ui.close()
```