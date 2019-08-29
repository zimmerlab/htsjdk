package htsjdk.samtools;

import htsjdk.HtsjdkTest;
import htsjdk.samtools.cram.compression.*;
import htsjdk.samtools.cram.compression.rans.RANS;
import htsjdk.samtools.cram.encoding.ByteArrayLenEncoding;
import htsjdk.samtools.cram.encoding.external.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.*;
import htsjdk.samtools.util.Tuple;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

public class CRAMAllEncodingStrategiesTest extends HtsjdkTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/htsjdk/samtools/cram");

    // TODO: need a better test file; this has mate validation errors
    // TODO: SAM validation error: ERROR: Read name20FUKAAXX100202:2:1:20271:61529,
    // TODO: Mate Alignment start (9999748) must be <= reference sequence length (200) on reference 20
    //final File cramSourceFile = new File(TEST_DATA_DIR, "NA12878.20.21.1-100.100-SeqsPerSlice.0-unMapped.cram");
    //final File referenceFile = new File(TEST_DATA_DIR, "human_g1k_v37.20.21.1-100.fasta");
    //final File cramSourceFile = new File("/Users/cnorman/projects/gatk/src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.samtools.cram");
    //final File referenceFile = new File("/Users/cnorman/projects/gatk/src/test/resources/large/human_g1k_v37.20.21.fasta");
    final File cramSourceFile = new File("/Users/cnorman/projects/gatk/src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram");
    final File referenceFile = new File("/Users/cnorman/projects/gatk/src/test/resources/large/human_g1k_v37.20.21.fasta");
    //final File cramSourceFile = new File("/Users/cnorman/projects/testdata/samn/DDP_ATCP_265_2.cram");
    //final File referenceFile = new File("/Users/cnorman/projects/references/hg38/Homo_sapiens_assembly38.fasta");

    @Test
    public final void testBestEncodingStrategy() throws IOException {
        System.out.println(String.format("Test file size: %,d", Files.size(cramSourceFile.toPath())));
        // src/test/resources/htsjdk/samtools/cram/json/CRAMEncodingMapProfileBEST.json has the encoding map used by this strategy
        final File encodingStrategyFile = new File("src/test/resources/htsjdk/samtools/cram/json/CRAMEncodingStrategyTemplate.json");
        final CRAMEncodingStrategy testStrategy = CRAMEncodingStrategy.readFromPath(encodingStrategyFile.toPath());
        final File tempOutCRAM = File.createTempFile("bestEncodingStrategyTest", ".cram");
        System.out.println(String.format("Output file: %s", tempOutCRAM.toPath()));
        final long fileSize = testWithEncodingStrategy(testStrategy, cramSourceFile, tempOutCRAM, referenceFile);
        assertRoundTripFidelity(cramSourceFile, tempOutCRAM, referenceFile);
        tempOutCRAM.delete();
        System.out.println(String.format("Size: %,d Strategy %s", fileSize, testStrategy));
    }

    @Test
    public final void testAllEncodingStrategyCombinations() throws IOException {
        final Map<Integer, CRAMEncodingStrategy> encodingStrategyByTest = new HashMap<>();
        final Map<Integer, String> encodingParamsByTest = new HashMap<>();
        final Map<Integer, Long> fileSizeByTest = new HashMap<>();

        System.out.println(String.format("Test file size: %,d", Files.size(cramSourceFile.toPath())));
        int testCount = 1;

        // the larger reads/slice and slices/container values are only interesting if there
        // are enough records in the input file to cause these thresholds to be crossed
        //      for (final int slicesPerContainer : Arrays.asList(1, 3)) {
        for (final int gzipCompressionLevel : Arrays.asList(5, 9)) {
            for (final int readsPerSlice : Arrays.asList(10000, 20000)) {
                for (final int slicesPerContainer : Arrays.asList(1)) {
                    for (final DataSeries dataSeries : enumerateDataSeries()) {
                        for (final EncodingDescriptor encodingDescriptor : enumerateEncodingDescriptorsFor(dataSeries)) {
                            for (final ExternalCompressor compressor : enumerateExternalCompressors(gzipCompressionLevel)) {
                                final CRAMEncodingStrategy testStrategy = createEncodingStrategyForParams(
                                        gzipCompressionLevel,
                                        readsPerSlice,
                                        slicesPerContainer,
                                        dataSeries,
                                        encodingDescriptor,
                                        compressor);

                                final File tempOutCRAM = File.createTempFile("allEncodingStrategyCombinations", ".cram");
                                final long fileSize = testWithEncodingStrategy(testStrategy, cramSourceFile, tempOutCRAM, referenceFile);
                                assertRoundTripFidelity(cramSourceFile, tempOutCRAM, referenceFile);
                                tempOutCRAM.delete();
                                final File mapPath = new File(testStrategy.getCustomCompressionMapPath());
                                mapPath.delete();

                                final String testSummary = String.format(
                                        "Size: %,d Test: %,d Series: %s Encoding: %s Compressor: %s %s",
                                        fileSize,
                                        testCount,
                                        dataSeries,
                                        encodingDescriptor.getEncodingID(),
                                        compressor,
                                        testStrategy);
                                System.out.println(testSummary);
                                encodingParamsByTest.put(testCount, testSummary);
                                encodingStrategyByTest.put(testCount, testStrategy);
                                fileSizeByTest.put(testCount, fileSize);
                                testCount++;
                            }
                        }
                    }
                }
            }
            System.out.println();
        }

        // sort and display encoding strategies ordered by ascending result file size
        final List<Tuple<Long, Integer>> bestEncodingStrategies =
                fileSizeByTest.entrySet()
                .stream()
                .map(e -> new Tuple<>(e.getValue(), e.getKey()))
                .collect(Collectors.toList());
        bestEncodingStrategies.sort(Comparator.comparing(t -> t.a));
        System.out.println(String.format("%d tests, sorted by result size:", testCount));
        bestEncodingStrategies
                // take the 500 best results
                .stream().limit(50)
                .forEach((Tuple<Long, Integer> t) ->
                        System.out.println(String.format("Size: %,d Test: %d Params: %s Encoding: %s",
                                t.a,
                                t.b,
                                encodingParamsByTest.get(t.b),
                                encodingStrategyByTest.get(t.b)))
                );
    }

    public void assertRoundTripFidelity(final File cramSourceFile, final File tempOutCRAM, final File referenceFile) {
        try (final CRAMFileReader origReader = new CRAMFileReader(cramSourceFile, new ReferenceSource(referenceFile));
             final CRAMFileReader copyReader = new CRAMFileReader(tempOutCRAM, new ReferenceSource(referenceFile))) {
            final SAMRecordIterator origIterator = origReader.getIterator();
            final SAMRecordIterator copyIterator = copyReader.getIterator();
            while (origIterator.hasNext() && copyIterator.hasNext()) {
                Assert.assertEquals(copyIterator.next(), origIterator.next());
            }
            Assert.assertEquals(origIterator.hasNext(), copyIterator.hasNext());
        }
    }

    public final CRAMEncodingStrategy createEncodingStrategyForParams(
           final int  gzipCompressionLevel,
           final int readsPerSlice,
           final int slicesPerContainer,
           final DataSeries ds,
           final EncodingDescriptor encodingDescriptor,
           final ExternalCompressor compressor) throws IOException {
        final CRAMEncodingStrategy encodingStrategy = new CRAMEncodingStrategy();
        encodingStrategy.setGZIPCompressionLevel(gzipCompressionLevel);
        encodingStrategy.setRecordsPerSlice(readsPerSlice);
        encodingStrategy.setSlicesPerContainer(slicesPerContainer);
        final CompressionHeaderEncodingMap encodingMap = createEncodingMapExternalEncodingVariationFor(ds, encodingDescriptor, compressor);
        final File tempEncodingMapFile = File.createTempFile("testEncodingMap", ".json");
        tempEncodingMapFile.deleteOnExit();
        final Path encodingMapPath = tempEncodingMapFile.toPath();
        encodingMap.writeToPath(encodingMapPath);
        encodingStrategy.setEncodingMap(encodingMapPath);

        // save the map to retrieve it in case its the best one
        //final File tempStrategyFile = File.createTempFile("testEncodingMap", ".json");
        //tempStrategyFile.deleteOnExit();
        //final Path strategyPath = tempStrategyFile.toPath();
        //encodingStrategy.writeToPath(strategyPath);

        return encodingStrategy;
    }

    // write with encoding params and return the size of the generated file
    private long testWithEncodingStrategy(
            final CRAMEncodingStrategy cramEncodingStrategy,
            final File inputCRAM,
            final File tempOutputCRAM,
            final File referenceFile) throws IOException {
        final File tempOutFile = File.createTempFile("encodingStrategiesTest", ".cram");
        tempOutFile.deleteOnExit();
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(referenceFile)
                .validationStringency((ValidationStringency.SILENT))
                .open(inputCRAM);
             final FileOutputStream fos = new FileOutputStream(tempOutputCRAM)) {
            final CRAMFileWriter cramWriter = new CRAMFileWriter(
                    cramEncodingStrategy,
                    fos,
                    null,
                    true,
                    new ReferenceSource(referenceFile),
                    reader.getFileHeader(),
                    tempOutputCRAM.getName());
            final SAMRecordIterator inputIterator = reader.iterator();
            while (inputIterator.hasNext()) {
                cramWriter.addAlignment(inputIterator.next());
            }
            cramWriter.close();
        }
        return Files.size(tempOutputCRAM.toPath());
    }

    public List<ExternalCompressor> enumerateExternalCompressors(final int gzipCompressionLevel) {
        return Arrays.asList(
                new GZIPExternalCompressor(gzipCompressionLevel),
                new RANSExternalCompressor(RANS.ORDER.ZERO),
                new RANSExternalCompressor(RANS.ORDER.ONE),
                new LZMAExternalCompressor(),
                new BZIP2ExternalCompressor()
        );
    }

    private List<DataSeries> enumerateDataSeries() {
        final List<DataSeries> seriesToUse = new ArrayList<>();
        // skip the ones this implementation doesn't use
        for (final DataSeries ds : DataSeries.values()) {
            if (ds != DataSeries.TM_TestMark
                    && ds != DataSeries.TV_TestMark
                    && ds != DataSeries.QQ_scores
                    && ds != DataSeries.BB_bases) {
                seriesToUse.add(ds);
            }
        }
        return seriesToUse;
    }

    // For now, just use the descriptors that are the default for this implementation, since not all
    // encodings are interchangeable
    private Set<EncodingDescriptor> enumerateEncodingDescriptorsFor(final DataSeries ds) {
        final Set<EncodingDescriptor> descriptors = new HashSet<>();
        descriptors.add(new CompressionHeaderEncodingMap(new CRAMEncodingStrategy()).getEncodingDescriptorForDataSeries(ds));
        if (ds == DataSeries.RN_ReadName) {
            descriptors.add(
                    (new ByteArrayLenEncoding(
                        new ExternalIntegerEncoding(ds.getExternalBlockContentId()),
                        new ExternalByteArrayEncoding(ds.getExternalBlockContentId()))).toEncodingDescriptor());
        }
        return descriptors;
    }

    public CompressionHeaderEncodingMap createEncodingMapExternalEncodingVariationFor(
            final DataSeries ds,
            final EncodingDescriptor encodingDescriptor,
            final ExternalCompressor compressor) {
        // create a default encoding map and update it with the requested EXTERNAL descriptor/compressor
        final CompressionHeaderEncodingMap encodingMap = new CompressionHeaderEncodingMap(new CRAMEncodingStrategy());
        encodingMap.putExternalEncoding(ds, encodingDescriptor, compressor);
        return encodingMap;
    }

}
