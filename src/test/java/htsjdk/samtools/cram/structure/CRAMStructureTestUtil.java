package htsjdk.samtools.cram.structure;

import htsjdk.HtsjdkTest;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.build.ContainerFactory;
import htsjdk.samtools.cram.ref.ReferenceContext;
import org.testng.Assert;
import org.testng.annotations.DataProvider;

import java.util.ArrayList;
import java.util.List;

public class CRAMStructureTestUtil extends HtsjdkTest {
    public static final int READ_LENGTH_FOR_TEST_RECORDS = 123;

    private static final SAMFileHeader header = initializeSAMFileHeaderForTests();

    private static SAMFileHeader initializeSAMFileHeaderForTests() {
        // arbitrary names and length.  Just ensure we have 10 different valid refs.

        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>();
        sequenceRecords.add(new SAMSequenceRecord("0", 10));
        sequenceRecords.add(new SAMSequenceRecord("1", 10));
        sequenceRecords.add(new SAMSequenceRecord("2", 10));
        sequenceRecords.add(new SAMSequenceRecord("3", 10));
        sequenceRecords.add(new SAMSequenceRecord("4", 10));
        sequenceRecords.add(new SAMSequenceRecord("5", 10));
        sequenceRecords.add(new SAMSequenceRecord("6", 10));
        sequenceRecords.add(new SAMSequenceRecord("7", 10));
        sequenceRecords.add(new SAMSequenceRecord("8", 10));
        sequenceRecords.add(new SAMSequenceRecord("9", 10));

        final SAMFileHeader header = new SAMFileHeader(new SAMSequenceDictionary(sequenceRecords));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return header;
    }

    public static SAMFileHeader getSAMFileHeaderForTests() {
        return header;
    }

    public static CRAMRecord createMappedRecord(final int index, final int referenceIndex, final int alignmentStart) {
        return new CRAMRecord(
                1,
                index,
                SAMFlag.SECOND_OF_PAIR.intValue(),
                0,
                "A READ NAME",
                READ_LENGTH_FOR_TEST_RECORDS,
                referenceIndex,
                alignmentStart,
                0,
                30,
                "!!!".getBytes(),
                "AAA".getBytes(),
                null,
                null,
                1,
                0,
                SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                SAMRecord.NO_ALIGNMENT_START,
                -1);
    }

    public static CRAMRecord createUnmappedPlacedRecord(final int index, final int referenceIndex, final int alignmentStart) {
        return new CRAMRecord(
                1,
                index,
                SAMFlag.READ_UNMAPPED.intValue(),
                0,
                "A READ NAME",
                READ_LENGTH_FOR_TEST_RECORDS,
                referenceIndex,
                alignmentStart,
                0,
                30,
                "!!!".getBytes(),
                "AAA".getBytes(),
                null,
                null,
                1,
                0,
                SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX,
                SAMRecord.NO_ALIGNMENT_START,
                -1);
    }

    public static CRAMRecord createUnmappedUnplacedRecord(final int index) {
        return createUnmappedPlacedRecord(index, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, SAMRecord.NO_ALIGNMENT_START);
    }

    public static List<CRAMRecord> getSingleRefRecords(final int recordCount, final int singleReferenceIndex) {
        final List<CRAMRecord> records = new ArrayList<>();
        for (int i = 0; i < recordCount; i++) {
            // set half unmapped-but-placed, to show that it does not make a difference
            if (i % 2 == 0) {
                records.add(createUnmappedPlacedRecord(i, singleReferenceIndex, i + 1));
            } else {
                records.add(createMappedRecord(i, singleReferenceIndex, i + 1));
            }
        }
        return records;
    }

    public static List<CRAMRecord> getMultiRefRecords(final int recordCount) {
        final List<CRAMRecord> records = new ArrayList<>();
        for (int i = 0; i < recordCount; i++) {
            // set half unmapped-but-placed, to show that it does not make a difference
            if (i % 2 == 0) {
                records.add(createUnmappedPlacedRecord(i, i, i + 1));
            } else {
                records.add(createMappedRecord(i, i, i + 1));
            }
        }
        return records;
    }

    public static List<CRAMRecord> getUnplacedRecords(final int recordCount) {
        final List<CRAMRecord> records = new ArrayList<>();
        for (int i = 0; i < recordCount; i++) {
            final CRAMRecord record = createUnmappedUnplacedRecord(i);
            records.add(record);
        }
        return records;
    }


    // these two sets of records are "half" unplaced: they have either a valid reference index or start position,
    // but not both.  We treat these weird edge cases as unplaced.

    public static List<CRAMRecord> getHalfUnplacedNoRefRecords(final int recordCount) {
        final List<CRAMRecord> records = new ArrayList<>();
        for (int i = 0; i < recordCount; i++) {
            records.add(createUnmappedPlacedRecord(i, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, i + 1));
        }
        return records;
    }

    public static List<CRAMRecord> getHalfUnplacedNoStartRecords(final int recordCount, final int referenceIndex) {
        final List<CRAMRecord> records = new ArrayList<>();
        for (int i = 0; i < recordCount; i++) {
            records.add(createUnmappedPlacedRecord(i, referenceIndex, SAMRecord.NO_ALIGNMENT_START));
        }
        return records;
    }

    public static List<CRAMRecord> getSingleRefRecordsWithOneUnmapped(final int testRecordCount, final int mappedReferenceIndex) {
        final List<CRAMRecord> retval = getSingleRefRecords(testRecordCount - 1, mappedReferenceIndex);
        retval.add(createUnmappedUnplacedRecord(testRecordCount - 1));
        return retval;
    }

    public static List<CRAMRecord> getMultiRefRecordsWithOneUnmapped(final int testRecordCount) {
        final List<CRAMRecord> retval = getMultiRefRecords(testRecordCount - 1);
        retval.add(createUnmappedUnplacedRecord(testRecordCount - 1));
        return retval;
    }

    public static List<Container> getMultiRefContainersForStateTest(final long firstContainerByteOffset) {
        final ContainerFactory factory = new ContainerFactory(
                getSAMFileHeaderForTests(),
                new CRAMEncodingStrategy().setRecordsPerSlice(10));
        final List<Container> testContainers = new ArrayList<>(3);

        final List<CRAMRecord> records = new ArrayList<>();

        int index = 0;
        records.add(createMappedRecord(index, index, index + 1));
        final Container container0 = factory.buildContainer(records, firstContainerByteOffset);

        index++;
        records.add(createMappedRecord(index, index, index + 1));
        final Container container1 = factory.buildContainer(records, firstContainerByteOffset + 1);

        index++;
        records.add(createUnmappedUnplacedRecord(index));
        final Container container2 = factory.buildContainer(records, firstContainerByteOffset + 2);

        testContainers.add(container0);
        testContainers.add(container1);
        testContainers.add(container2);
        return testContainers;
    }

    private static Slice getIndexInitializedSlice() {
        final ReferenceContext refContext = new ReferenceContext(0);

        final Slice slice = new Slice(refContext);
        slice.byteOffsetFromCompressionHeaderStart = 1;
        slice.containerByteOffset = 1;
        slice.byteSize = 1;
        slice.index = 1;

        return slice;
    }

    private static Slice getNoContainerOffsetSlice() {
        final Slice noContainerOffset = getIndexInitializedSlice();
        noContainerOffset.containerByteOffset = Slice.UNINITIALIZED_INDEXING_PARAMETER;
        return noContainerOffset;
    }

    private static Slice getNoOffsetFromContainerSlice() {
        final Slice noOffsetFromContainer = getIndexInitializedSlice();
        noOffsetFromContainer.byteOffsetFromCompressionHeaderStart = Slice.UNINITIALIZED_INDEXING_PARAMETER;
        return noOffsetFromContainer;
    }

    private static Slice getNoSizeSlice() {
        final Slice noSize = getIndexInitializedSlice();
        noSize.byteSize = Slice.UNINITIALIZED_INDEXING_PARAMETER;
        return noSize;
    }

    private static Slice getNoIndexSlice() {
        final Slice noIndex = getIndexInitializedSlice();
        noIndex.index = Slice.UNINITIALIZED_INDEXING_PARAMETER;
        return noIndex;
    }

    @DataProvider(name = "uninitializedBAIParameterTestCases")
    static Object[][] uninitializedBAIParameterTestCases() {

        return new Object[][] {
                { getNoContainerOffsetSlice() },
                { getNoOffsetFromContainerSlice() },
                { getNoIndexSlice() }
        };
    }

    @DataProvider(name = "uninitializedCRAIParameterTestCases")
    static Object[][] uninitializedCRAIParameterTestCases() {

        return new Object[][] {
                { getNoContainerOffsetSlice() },
                { getNoOffsetFromContainerSlice() },
                { getNoSizeSlice() }
        };
    }

    // assert that slices and containers have values equal to what the caller expects

    public static void assertSliceState(final Slice slice,
                                        final ReferenceContext expectedReferenceContext,
                                        final int expectedAlignmentStart,
                                        final int expectedAlignmentSpan,
                                        final int expectedRecordCount,
                                        final int expectedBaseCount) {
        Assert.assertEquals(slice.getReferenceContext(), expectedReferenceContext);
        Assert.assertEquals(slice.alignmentStart, expectedAlignmentStart);
        Assert.assertEquals(slice.alignmentSpan, expectedAlignmentSpan);
        Assert.assertEquals(slice.nofRecords, expectedRecordCount);
        Assert.assertEquals(slice.bases, expectedBaseCount);
    }

    public static void assertSliceState(final Slice slice,
                                        final ReferenceContext expectedReferenceContext,
                                        final int expectedAlignmentStart,
                                        final int expectedAlignmentSpan,
                                        final int expectedRecordCount,
                                        final int expectedBaseCount,
                                        final int expectedGlobalRecordCounter) {
        assertSliceState(slice, expectedReferenceContext, expectedAlignmentStart, expectedAlignmentSpan, expectedRecordCount, expectedBaseCount);
        Assert.assertEquals(slice.globalRecordCounter, expectedGlobalRecordCounter);
    }

    public static void assertContainerState(final Container container,
                                            final ReferenceContext expectedReferenceContext,
                                            final int expectedAlignmentStart,
                                            final int expectedAlignmentSpan,
                                            final long expectedByteOffset) {
        Assert.assertEquals(container.getReferenceContext(), expectedReferenceContext);
        Assert.assertEquals(container.alignmentStart, expectedAlignmentStart);
        Assert.assertEquals(container.alignmentSpan, expectedAlignmentSpan);
        Assert.assertEquals(container.byteOffset, expectedByteOffset);
    }

    public static void assertContainerState(final Container container,
                                            final ReferenceContext expectedReferenceContext,
                                            final int expectedAlignmentStart,
                                            final int expectedAlignmentSpan,
                                            final int expectedRecordCount,
                                            final int expectedBaseCount,
                                            final int expectedGlobalRecordCounter,
                                            final long expectedByteOffset) {
        assertContainerState(container, expectedReferenceContext, expectedAlignmentStart, expectedAlignmentSpan, expectedByteOffset);

        Assert.assertEquals(container.nofRecords, expectedRecordCount);
        Assert.assertEquals(container.bases, expectedBaseCount);
        Assert.assertEquals(container.globalRecordCounter, expectedGlobalRecordCounter);

        Assert.assertEquals(container.getSlices().length, 1);

        // verify the underlying slice too

        assertSliceState(container.getSlices()[0], expectedReferenceContext, expectedAlignmentStart, expectedAlignmentSpan,
                expectedRecordCount, expectedBaseCount, expectedGlobalRecordCounter);
    }
}