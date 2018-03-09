

# Gaussian processes, sampling and marginalization


The goal in this tutorial is to experiment with the notions of Gaussian process sampling and marginalization, 
and also learn how to evaluate the probability of shapes according to a given model.

As in the previous tutorials, we start by importing some commonly used objects and initializing the system. 

```scala mdoc
import scalismo.geometry.{_, Vector=>SpatialVector}
import scalismo.common._
import scalismo.ui.api._

scalismo.initialize()
implicit val rng = scalismo.utils.Random(42)

val ui = ScalismoUI()
```



#### Discrete and Continuous Gaussian processes

We have seen in the last tutorial that a Point Distribution Model (PDM) can be interpreted as a Discrete Gaussian process, modelling deformation fields, 
defined on a reference mesh. To continue our exploration of Gaussian processes, we therefore start by loading (and visualizing) an existing PDM 


```scala mdoc
val model = StatismoIO.readStatismoMeshModel(new File("datasets/bfm.h5")).get

val modelGroup = ui.createGroup("modelGroup")
val ssmView = ui.show(modelGroup, model, "model")
```

We can access the Gaussian Process and sample from it 

```scala mdoc
val gp = model.gp
val sampleDF : DiscreteField[_3D,UnstructuredPointsDomain[_3D], SpatialVector[_3D]] = model.gp.sample

val sampleGroup = ui.createGroup("sample")  
ui.show(sampleGroup, sampleDF, "discreteSample")
```

As you can see here the sampled vector field is a **discrete function**, defined over a **discrete set of points**. 
This is due to the fact that our Gaussian Process is stored in a file and was therefore discretized over the points of the reference face.

As seen in the previous tutorial, we could now interpolate the sample ```sampleDf``` to obtain a continuous version of the deformation field. 
A more convenient approach is, however, to interpolate directly the Gaussian process: 
 
```scala mdoc
val interpolator = NearestNeighborInterpolator3D()
val contGP = model.gp.interpolate(interpolator)
```

Now, when sampling from the continuous GP, we obtain a vector-valued function, which is defined on the entire 3D Space: 

```scala mdoc
val contSample: Field[_3D,SpatialVector[_3D]] = contGP.sample
```

*Attention: While the interpolated Gaussian process is now defined on the entire 3D Space, the interpolation really only makes sense close to the mesh points*.


## From continuous to discrete: marginalization

While having a continuous GP, which models a distribution over continuous vector fields is conceptually satisfying, in practice, we are always interested in 
the distribution on a finite set of points. But the advantage of having a continuous GP, is now, that we can now get a sample for an *any* finite set of points. 
To illustrate this, we could, for example obtain a sample, which is defined on all the points of the original reference mesh.

```scala mdoc
val fullSample = contGP.sampleAtPoints(model.referenceMesh)
val fullSampleView = ui.show(sampleGroup, fullSample, "fullSample")
```

But any other discrete set of points, represented by an object of type ```DiscreteDomain``` is possible. We can, for example
use sample only a single point: 

```scala mdoc
fullSampleView.remove()
val singlePointDomain : DiscreteDomain[_3D] = UnstructuredPointsDomain(IndexedSeq(model.referenceMesh.point(PointId(8156))))
val singlePointSample = contGP.sampleAtPoints(singlePointDomain) 
ui.show(sampleGroup, singlePointSample, "singlePointSample")
```
(This should show a vector at the tip of the nose, which, could also be behind the face)

The marginalization property of a Gaussian process makes it even possible to obtain a distribution at any finite set of points. We can 
obtain this distribution, by calling the method ```marginal``` on the Gaussian process instance. 

```scala mdoc
val rightEyePt: Point[_3D] = model.referenceMesh.point(PointId(4281))
val leftEyePt: Point[_3D] = model.referenceMesh.point(PointId(11937))
val dom = UnstructuredPointsDomain(IndexedSeq(rightEyePt,leftEyePt))
val marginal : DiscreteGaussianProcess[_3D, UnstructuredPointsDomain[_3D], SpatialVector[_3D]] = contGP.marginal(dom)
```

As we see, the result of marginalization is again a discrete Gaussian process. 
We can now again obtain a sample for this GAussian process
```scala mdoc
val sample : DiscreteVectorField[_3D, _3D] = marginal.sample 
show(sample, "marginal_sample") 
```
As you can see, the sample is now a discrete deformation field containing only 2 deformations over the desired points.  


It seems that we are back where we started. But note that we have now choosen a completly different set of points
on which the Gaussian process is defined. This is important, as we can choose for any application that discretization of the Gaussian process, which is most useful.
  

## Marginal of a statistical mesh model

Given that a *StatisticalMeshModel* is in reality just a wrapper around a GP, it naturally allows for marginalization as well: 

```scala mdoc
val noseTipModel : StatisticalMeshModel = model.marginal(IndexedSeq(PointId(8156)))
```

Notice in this case, how the passed argument to the marginal function is an indexed sequence of point **identifiers** instead of a discrete domain.  
This is due to the fact that we are marginalizing a discrete Gaussian process. 
Since the domain of the GP is already discrete, marginalization in this case is done by selecting a subset of the discrete domain. Hence the use of point identifiers instead of 3D coordinates.


Not surprisingly, we can again sample from this nose tip model: 

```scala mdoc
val tipSample : TriangleMesh[_3D] = noseTipModel.sample
println("nb mesh points " + tipSample.numberOfPoints)
```

Given that the marginal model is a *StatisticalMeshModel*, sampling from it returns a ```TriangleMesh```. When inspecting the points of the returned sample, we see that it contains only one point, the nose tip. 

#### Nose marginal

Let's suppose that we have a full model of the face, but are only interested in the shape variations around the nose. Marginalization let's us achieve this easily.
To do so, we first need to extract the list of point identifiers of the nose points:

```scala mdoc
val middleNose = model.referenceMesh.point(PointId(8152))
val nosePtIDs = model.referenceMesh.pointSet.pointsWithId.filter( ptAndId => {
  (pt - middleNose).norm > 40
  }
)   
```
We can now use the point ids to marginalize our shape model:

```scala mdoc
val noseModel = model.marginal(nosePtIDs.toIndexedSeq)
show(noseModel, "noseModel")
```

## Probability of shapes and deformations: 

It is often interesting to assess how probable a model instance is. This can be done in Scalismo by means of the method ```pdf ``` (which stands for probability density function), 
which is defined on the class ```GaussianProcess``` and ```StatisticalMeshModel```.


```scala:silent
val defSample = noseModel.gp.sample
noseModel.gp.pdf(defSample)
```

The value of the *pdf* is often not interesting as such. But it allows us to compare the likelihood of different instances, by comparing their density value.
For numerical reasons, we usually work with the log probability:

```tut:silent
val defSample1 = noseModel.gp.sample
val defSample2 = noseModel.gp.sample

val logPDF1 = noseModel.gp.logpdf(defSample1)
val logPDF2 = noseModel.gp.logpdf(defSample2)

println("defSample1 is " + (if (logPDF1 > logPDF2) " more " else " less ") likely than defSample2 )
```


```scala mdoc:invisible
ui.close
```