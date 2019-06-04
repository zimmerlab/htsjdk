package htsjdk.samtools.cram.encoding;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.CramRecordTestHelper;
import htsjdk.samtools.cram.ref.ReferenceContext;
import htsjdk.samtools.cram.structure.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Map;

public class MultiRefSliceAlignmentSpanReaderTest extends CramRecordTestHelper {

    private List<CramCompressionRecord> initSpanRecords() {
        final List<CramCompressionRecord> initialRecords = createRecords();

        // note for future refactoring
        // createRecord() calls Sam2CramRecordFactory.createCramRecord()
        // which is the only way to set a record's readFeatures (except via read codec)
        // which would otherwise be null

        final String commonRead = "AAA";

        // span 1:1,3
        initialRecords.get(0).sequenceId = 1;
        initialRecords.get(0).alignmentStart = 1;
        initialRecords.get(0).readBases = commonRead.getBytes();
        initialRecords.get(0).readLength = commonRead.length();

        // span 2:2,4
        initialRecords.get(1).sequenceId = 2;
        initialRecords.get(1).alignmentStart = 2;
        initialRecords.get(1).readBases = commonRead.getBytes();
        initialRecords.get(1).readLength = commonRead.length();
        initialRecords.get(1).setSegmentUnmapped(true);     // unmapped but placed

        // span 1:3,5
        initialRecords.get(2).sequenceId = 1;
        initialRecords.get(2).alignmentStart = 3;
        initialRecords.get(2).readBases = commonRead.getBytes();
        initialRecords.get(2).readLength = commonRead.length();

        // span <unmapped>
        initialRecords.get(3).sequenceId = SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
        initialRecords.get(3).alignmentStart = 7;
        initialRecords.get(3).readBases = commonRead.getBytes();
        initialRecords.get(3).readLength = commonRead.length();

        // span totals -> 1:1,5 and 2:2,4

        return initialRecords;
    }

    @Test
    public void testSpansCoordinateSorted() {
        final List<CramCompressionRecord> initialRecords = initSpanRecords();

        // note for future refactoring
        // createHeader(records) calls CompressionHeaderBuilder.setTagIdDictionary(buildTagIdDictionary(records));
        // which is the only way to set a record's tagIdsIndex
        // which would otherwise be null

        // NOTE: multiref alignment spans are ony used for CRAI indexing, and only make sense when records are
        // coordinate sorted, so we only test with coordinateSorted = true;
        final CompressionHeader header = createHeader(initialRecords, true);
        final Slice slice = Slice.buildSlice(initialRecords, header);
        final Map<ReferenceContext, AlignmentSpan> spans = slice.getMultiRefAlignmentSpans(header, ValidationStringency.DEFAULT_STRINGENCY);

        Assert.assertEquals(spans.size(), 3);
        Assert.assertEquals(spans.get(new ReferenceContext(1)), new AlignmentSpan(1, 5, 2, 0));
        Assert.assertEquals(spans.get(new ReferenceContext(2)), new AlignmentSpan(2, 3, 0, 1));
        Assert.assertEquals(spans.get(ReferenceContext.UNMAPPED_UNPLACED_CONTEXT), AlignmentSpan.UNPLACED_SPAN);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testSpansNonCoordinateSorted() {
        final List<CramCompressionRecord> initialRecords = initSpanRecords();

        // NOTE: multiref alignment spans are ony used for CRAI indexing, and only make sense when records are
        // coordinate sorted, so test that we reject coordinateSorted = true;
        final CompressionHeader header = createHeader(initialRecords, false);
        final Slice slice = Slice.buildSlice(initialRecords, header);

        slice.getMultiRefAlignmentSpans(header, ValidationStringency.DEFAULT_STRINGENCY);
    }
}
