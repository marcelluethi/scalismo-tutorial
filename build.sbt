organization  := "ch.unibas.cs.gravis"

name := """scalismo-tutorial"""
version       := "0.16.0"

scalaVersion  := "2.11.7"

scalacOptions := Seq("-unchecked", "-deprecation", "-encoding", "utf8")

resolvers += Resolver.bintrayRepo("unibas-gravis", "maven")

resolvers += Opts.resolver.sonatypeSnapshots

enablePlugins(TutPlugin)

libraryDependencies  ++= Seq(
            "ch.unibas.cs.gravis" % "scalismo-native-all" % "4.0.+",
            "ch.unibas.cs.gravis" %% "scalismo-ui" % "0.12.0"
)
