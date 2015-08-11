# PRIDE Spectra Clustering #

The PRIDE Spectra Clustering API is an easy-to-use Java API to cluster MS/MS spectra. It is currently used to create PRIDE Cluster - a resource clustering all identified spectra from the [PRIDE repository](http://www.ebi.ac.uk/pride).

When you use it, please cite: cite the following publication:

[J. Griss et al., Nature Methods 2013](http://dx.doi.org/10.1038/nmeth.2343). [PubMed record](http://www.ncbi.nlm.nih.gov/pubmed/23361086).

# Using the PRIDE Spectra Clustering API #

The PRIDE Spectra Clustering API was developed as a Maven project. A complete tutorial showing how the PRIDE Spectra Clustering API can be used is available in the [Wiki](https://code.google.com/p/pride-spectra-clustering/wiki/ClusteringTutorial).

## Maven ##

To use the PRIDE Spectra Clustering API in your maven project add the following dependency:
```
<dependency>
  <groupId>uk.ac.ebi.pride.tools</groupId>
  <artifactId>pride-spectra-clustering</artifactId>
  <version>1.1</version>
</dependency>
```
This dependency is available in the EBI's maven repository:
```
<repository>
  <id>ebi-repo</id>
  <name>The EBI internal repository</name>
  <url>http://www.ebi.ac.uk/~maven/m2repo</url>
  <releases>
    <enabled>true</enabled>
  </releases>
  <snapshots>
    <enabled>false</enabled>
  </snapshots>
</repository>
```

## JAR ##

Ready compiled JAR files are available in the [Downloads](https://code.google.com/p/pride-spectra-clustering/downloads/list) section.

## SVN ##

The complete source code of the PRIDE Cluster API can be checked out from the [SVN](https://code.google.com/p/pride-spectra-clustering/source/checkout).