## A Java API for high-throughput sequencing data (HTS) formats. 

Custom fork that is modified to allow access to private datastructures not available in the default repo.

### How to use the Jar in a maven project

#### Option 1: Add to local maven repo:

1. Build the Jar:
```bash
./gradlew
```

Jar location:
```bash
build/libs/htsjdk-<VERSION>.jar
```

2. Install the Jar to the local maven repo:

Replace `VERSION` with the version of the previously built Jar.

```bash
mvn install:install-file \
   -Dfile=htsjdk-<version>.jar \
   -DgroupId=com.github.samtools \
   -DartifactId=htsjdk \
   -Dversion=VERSION \
   -Dpackaging=jar \
   -DgeneratePom=true
```

#### Option 2: Use remote maven repo:

1. Configure maven authentication for github:

Add this to your `~/.m2/settings.xml` file and replace `USERNAME` with your github username and `TOKEN` with a personal github access token:

```xml
<settings>
    <servers>
        <server>
            <id>github-zimmerlab-htsjdk</id>
            <username>USERNAME</username>
            <password>TOKEN</password>
        </server>
    </servers>
</settings>
```

2. Add the github package repository to the maven project:

Add this to the `pom.xml` file in the projects root directory:

```xml
<repositories>
    <repository>
        <id>github-zimmerlab-htsjdk</id>
        <url>https://maven.pkg.github.com/zimmerlab/htsjdk</url>
        <snapshots>
            <enabled>true</enabled>
        </snapshots>
    </repository>
</repositories>
```

3. Add the dependency to the `pom.xml` file in the projects root directory:

Replace `VERSION` with the version of the jar you want to use (see [releases](https://github.com/orgs/zimmerlab/packages?repo_name=htsjdk)).

```xml
<dependency>
    <groupId>com.github.zimmerlab</groupId>
    <artifactId>htsjdk</artifactId>
    <version>VERSION</version>
</dependency>
```