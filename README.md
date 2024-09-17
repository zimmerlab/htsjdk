## A Java API for high-throughput sequencing data (HTS) formats. 

Custom fork that is modified to allow access to private datastructures not available in the default repo.

### Build custom Jar

Build jar:
```bash
./gradlew
```

Jar location:
```bash
build/libs/htsjdk-<version>.jar
```

### Use custom Jar in other maven project

Add to local maven repo:
```bash
mvn install:install-file \
   -Dfile=htsjdk-<version>.jar \
   -DgroupId=com.github.samtools \
   -DartifactId=htsjdk \
   -Dversion=4.1.1-SNAPSHOT \
   -Dpackaging=jar \
   -DgeneratePom=true
```