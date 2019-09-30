/**
 * ****************************************************************************
 * Copyright 2013 EMBL-EBI
 * <p/>
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * <p/>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p/>
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ****************************************************************************
 */
package htsjdk.samtools.cram.structure;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.BAIEntry;
import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.cram.CRAMException;
import htsjdk.samtools.cram.build.CRAMReferenceState;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.common.CramVersions;
import htsjdk.samtools.cram.common.Version;
import htsjdk.samtools.cram.io.InputStreamUtils;
import htsjdk.samtools.cram.ref.ReferenceContext;
import htsjdk.samtools.cram.structure.block.Block;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.LineReader;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.*;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Notes: Container will construct a container out of as many cramRecords as it is handed, respecting only
 * the maximum number of slices. The policy around how to break up lists of records into containers is enforced
 * by ContainerFactory.
 */
public class Container {
    private final ContainerHeader containerHeader;
    private final CompressionHeader compressionHeader;
    private final List<Slice> slices;

    // container's byte offset from the start of the containing stream, used for indexing
    private final long containerByteOffset;

    /**
     * Derive the container's {@link ReferenceContext} from its {@link Slice}s.
     *
     * A Single Reference Container contains only Single Reference Slices mapped to the same reference.
     * - set the Container's ReferenceContext to be the same as those slices
     * - set the Container's Alignment Start and Span to cover all slices
     *
     * A Multiple Reference Container contains only Multiple Reference Slices.
     * - set the Container's ReferenceContext to MULTIPLE_REFERENCE_CONTEXT
     * - unset the Container's Alignment Start and Span
     *
     * An Unmapped Container contains only Unmapped Slices.
     * - set the Container's ReferenceContext to UNMAPPED_UNPLACED_CONTEXT
     * - unset the Container's Alignment Start and Span
     *
     * Any other combination is invalid.
     *
//     * @param compressionHeader the CRAM {@link CompressionHeader} to assign to the Container
//     * @param containerSlices the constituent Slices of the Container
     * @param containerByteOffset the Container's byte offset from the start of the stream
     * @param globalRecordCounter the global record count for the first record in this container
     * @throws CRAMException for invalid Container states
     * @return the initialized Container
     */
    public Container(
            final CompressionHeader compressionHeader,
            final List<Slice> slices,
            final long containerByteOffset,
            final long globalRecordCounter) {

        this.compressionHeader = compressionHeader;
        this.slices = slices;
        this.containerByteOffset = containerByteOffset;

        final ReferenceContext commonRefContext = getDerivedReferenceContextFromSlices(slices);
        final AlignmentContext alignmentContext = getDerivedAlignmentContext(commonRefContext);

        int baseCount = 0;
        int blockCount = 0;
        int recordCount = 0;
        for (final Slice slice : slices) {
            recordCount += slice.getNumberOfRecords();
            blockCount += slice.getNumberOfBlocks();
            baseCount += slice.getBaseCount();
        }

        this.containerHeader = new ContainerHeader(
                alignmentContext,
                globalRecordCounter,
                blockCount,
                recordCount,
                baseCount);

        validateContainerReferenceContext();
        checkReferenceContexts(commonRefContext.getReferenceContextID());
    }

    private final AlignmentContext getDerivedAlignmentContext(final ReferenceContext commonRefContext) {
        int alignmentStart = SAMRecord.NO_ALIGNMENT_START;
        int alignmentSpan = AlignmentContext.NO_ALIGNMENT_SPAN;

        if (commonRefContext.isMappedSingleRef()) {
            int start = Integer.MAX_VALUE;
            // end is start + span - 1.  We can do slightly easier math instead.
            int endPlusOne = Integer.MIN_VALUE;

            for (final Slice slice : slices) {
                final AlignmentContext alignmentContext = slice.getAlignmentContext();
                start = Math.min(start, alignmentContext.getAlignmentStart());
                endPlusOne = Math.max(endPlusOne, alignmentContext.getAlignmentStart() + alignmentContext.getAlignmentSpan());
            }
            alignmentStart = start;
            alignmentSpan = endPlusOne - start;
        } else if (commonRefContext.isUnmappedUnplaced()) {
            return AlignmentContext.UNMAPPED_UNPLACED_CONTEXT;
        } else if (commonRefContext.isMultiRef()) {
            return AlignmentContext.MULTIPLE_REFERENCE_CONTEXT;
        }

        // since we're creating this container, ensure that it has a valid alignment context
        AlignmentContext.validateAlignmentContext(true, commonRefContext, alignmentStart, alignmentSpan);
        return new AlignmentContext(commonRefContext, alignmentStart, alignmentSpan);
    }

    // Note: this is the degenerate case of the CramContainerHeaderIterator, which for disq really only
    // cares about the byte offset...
    public Container(final ContainerHeader containerHeader, final long containerByteOffset) {
        this.containerHeader = containerHeader;
        this.containerByteOffset = containerByteOffset;
        compressionHeader = null;
        slices = Collections.EMPTY_LIST;
    }

    // Note: this is the case where we're reading a container from a stream. This reads in the container
    // header, and all slice blocks, but does not resolve the blocks into CRAMRecords, since we don't want
    // to do that until its necessary (and it might not be if, for example, we're indexing the container).
    public Container(final Version version, final InputStream inputStream, final long containerByteOffset) {
        containerHeader = ContainerHeader.readContainerHeader(version.major, inputStream);
        if (containerHeader.isEOF()) {
            compressionHeader = null;
            slices = Collections.EMPTY_LIST;
            this.containerByteOffset = containerByteOffset;
            return;
        }

        this.containerByteOffset = containerByteOffset;
        compressionHeader = new CompressionHeader(version.major, inputStream);

        this.slices = new ArrayList<>();
        for (int sliceCounter = 0; sliceCounter < containerHeader.getLandmarks().length; sliceCounter++) {
            final Slice slice = new Slice(
                    version.major,
                    compressionHeader,
                    inputStream,
                    containerByteOffset);
            slices.add(slice);
        }

        // sets index, byteOffset and byte size for each slice...
        distributeIndexingParametersToSlices();
        validateContainerReferenceContext();
        checkSliceReferenceContexts(getAlignmentContext().getReferenceContext().getReferenceContextID());
    }

    /**
     * Writes a complete {@link Container} with it's header to a {@link OutputStream}.
     *
     * @param version   the CRAM version to assume
     * @param outputStream  the stream to write to
     * @return the number of bytes written out
     */
    public int writeContainer(final Version version, final OutputStream outputStream) {

        // The first thing that needs to be written to the final output stream is the Container header,
        // but that contains the size of the compressed header, as well as the landmarks, which are the
        // offsets of the start of slices, and these are not known until the compressed slices blocks are
        // written out. So we first write the
        try (final ByteArrayOutputStream tempOutputStream = new ByteArrayOutputStream()) {

            // Before we can write out the container header, we need to update it with the length of the
            // compression header that follows, as well as with the landmark (slice offsets), so write the
            // compression header out to a temporary stream first to get it's length.
            getCompressionHeader().write(version, tempOutputStream);

            // ...then write out the slice blocks, computing the landmarks along the way
            final List<Integer> landmarks = new ArrayList<>();
            for (final Slice slice : getSlices()) {
                // landmark 0 = byte length of the compression header
                // landmarks after 0 = byte length of the compression header plus all slices before this one
                landmarks.add(tempOutputStream.size());
                slice.write(version.major, tempOutputStream);
            }
            getContainerHeader().setLandmarks(landmarks.stream().mapToInt(Integer::intValue).toArray());

            // compression header plus all slices, if any (EOF Containers do not; File Header Containers are handled above)
            getContainerHeader().setContainerBlocksByteSize(tempOutputStream.size());

            // Slices require the Container's landmarks and containerBlocksByteSize in case we're indexing
            distributeIndexingParametersToSlices();

            // ...then write the container header
            final int containerHeaderLength = getContainerHeader().writeContainerHeader(version.major, outputStream);

            // .. and finally, write the entire stream to the actual output stream, now that we know how big it is
            outputStream.write(tempOutputStream.toByteArray(),
                    0,
                    getContainerHeader().getContainerBlocksByteSize());

            return containerHeaderLength + getContainerHeader().getContainerBlocksByteSize();
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Reads the special container that contains the SAMFileHeader from a CRAM stream, and return just
     * the SAMFileHeader (we don't want to hand out the container since its not really a container in that
     * while it has a container header, it has compression header block, no slices, etc).
     * @param version
     * @param inputStream
     * @param id
     * @return
     */
    public static SAMFileHeader getSAMFileHeaderFromContainer(final Version version,
                                                              final InputStream inputStream,
                                                              final String id) {
        final ContainerHeader containerHeader = ContainerHeader.readContainerHeader(version.major, inputStream);
        final Block block;
        if (version.compatibleWith(CramVersions.CRAM_v3)) {
            final byte[] bytes = new byte[containerHeader.getContainerBlocksByteSize()];
            InputStreamUtils.readFully(inputStream, bytes, 0, bytes.length);
            block = Block.read(version.major, new ByteArrayInputStream(bytes));
            // ignore the rest of the container
        } else {
            /*
             * pending issue: container.containerBlocksByteSize inputStream 2 bytes shorter
             * than needed in the v21 test cram files.
             */
            block = Block.read(version.major, inputStream);
        }

        // SAMFileHeader block is prescribed by the spec to be gzip compressed
        // Use a temporary, single-use compressor cache for this block, since its
        // prescribed by the spec to be gzipped only and thus not perf sensitive.
        try (final InputStream blockStream = new ByteArrayInputStream(block.getUncompressedContent(new CompressorCache()))) {
            final ByteBuffer buffer = ByteBuffer.allocate(4);
            buffer.order(ByteOrder.LITTLE_ENDIAN);

            for (int i = 0; i < 4; i++) {
                buffer.put((byte) blockStream.read());
            }

            buffer.flip();
            final int size = buffer.asIntBuffer().get();

            final DataInputStream dataInputStream = new DataInputStream(blockStream);
            final byte[] bytes = new byte[size];
            dataInputStream.readFully(bytes);
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            try (final InputStream byteStream = new ByteArrayInputStream(bytes);
                 final LineReader lineReader = new BufferedLineReader(byteStream)) {
                return codec.decode(lineReader, id);
            }
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    public static long writeSAMFileHeaderContainer(final int major, final SAMFileHeader samFileHeader, final OutputStream os) {
        final byte[] data = CramIO.toByteArray(samFileHeader);
        // The spec recommends "reserving" 50% more space than is required by the header.
        final int length = Math.max(1024, data.length + data.length / 2);
        final byte[] blockContent = new byte[length];
        System.arraycopy(data, 0, blockContent, 0, Math.min(data.length, length));
        final Block block = Block.createRawFileHeaderBlock(blockContent);

        try (final ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream()) {
            block.write(major, byteArrayOutputStream);
            int containerBlocksByteSize = byteArrayOutputStream.size();
            // TODO: make sure this container is initialized correctly/fully (add a test for checksum ?)
            final ContainerHeader containerHeader = new ContainerHeader(
                    // we're forced to create an alignment context for this bogus container...
                    AlignmentContext.UNMAPPED_UNPLACED_CONTEXT,
                    containerBlocksByteSize,
                    0,
                    0,
                    0,
                    1,
                    //TODO: should the landmarks be specified ?
                    new int[]{},
                    0);
            final int containerHeaderByteSize = containerHeader.writeContainerHeader(major, os);
            os.write(byteArrayOutputStream.toByteArray(), 0, containerBlocksByteSize);
            return containerHeaderByteSize + containerHeader.getContainerBlocksByteSize();
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    //TODO: this unpacks all slices
    //TODO: note that this does not require a reference, which is good, since we need to be able to use it to
    // get raw CRAM records during indexing, and we don;t want that to require a reference
    //Visible for testing
    public List<CRAMRecord> getCRAMRecords(final ValidationStringency validationStringency, final CompressorCache compressorCache) {
        if (isEOF()) {
            return Collections.emptyList();
        }

        final ArrayList<CRAMRecord> records = new ArrayList<>(getContainerHeader().getNumberOfRecords());
        for (final Slice slice : getSlices()) {
            records.addAll(slice.getCRAMRecords(compressorCache, validationStringency));
        }
        return records;
    }

    // This public method should hand out slices, and let the caller get the SAMRecords out of the slice, for
    // symmetry with the constructor (which takes Slices).
    public List<SAMRecord> getSAMRecords(
            final ValidationStringency validationStringency,
            final CRAMReferenceState cramReferenceState,
            final CompressorCache compressorCache,
            final SAMFileHeader samFileHeader) {
        final List<CRAMRecord> cramRecords = getCRAMRecords(validationStringency, compressorCache);
        final List<SAMRecord> samRecords = new ArrayList<>(cramRecords.size());

        for (final Slice slice : getSlices()) {
            if (slice.getAlignmentContext().getReferenceContext().isMappedSingleRef()) {
                final byte[] referenceBases = cramReferenceState.getReferenceBases(
                        slice.getAlignmentContext().getReferenceContext().getReferenceSequenceID()
                );

                if (!slice.validateReferenceMD5(referenceBases)) {
                    final String msg = String.format(
                            "Reference sequence MD5 mismatch for slice: %s, expected MD5 %s",
                            slice.getAlignmentContext(),
                            String.format("%032x", new BigInteger(1, slice.getRefMD5())));
                    throw new CRAMException(msg);
                }
            }

            slice.normalize(cramRecords, cramReferenceState, 0, getCompressionHeader().getSubstitutionMatrix());

            for (final CRAMRecord cramRecord : cramRecords) {
                final SAMRecord samRecord = cramRecord.toSAMRecord(samFileHeader);
                samRecord.setValidationStringency(validationStringency);
                samRecords.add(samRecord);
            }
        }

        cramRecords.clear();
        return samRecords;
    }

    public ContainerHeader getContainerHeader() { return containerHeader; }
    public CompressionHeader getCompressionHeader() { return compressionHeader; }
    public AlignmentContext getAlignmentContext() { return containerHeader.getAlignmentContext(); }
    public long getContainerByteOffset() { return containerByteOffset; }
    public List<Slice> getSlices() { return slices; }
    public boolean isEOF() {
        return containerHeader.isEOF() && (getSlices() == null || getSlices().size() == 0);
    }

    /**
     * Retrieve the list of CRAI Index entries corresponding to this Container
     * @return the list of CRAI Index entries
     */
    public List<CRAIEntry> getCRAIEntries(final CompressorCache compressorCache) {
        if (isEOF()) {
            return Collections.emptyList();
        }

        //TODO: these might need to be merged, so that in the end there is only one entry per ref context ?
        // these should NOT need to be sorted ?? we need one per slice, several for multi-ref slices,
        // so no merging should be required
        return getSlices().stream()
                .map(s -> s.getCRAIEntries(compressorCache))
                .flatMap(List::stream)
                .sorted()
                .collect(Collectors.toList());
    }

    /**
     * Retrieve the list of BAIEntry Index entries corresponding to this Container
     * @return the list of BAIEntry Index entries
     */
    public List<BAIEntry> getBAIEntries(final CompressorCache compressorCache) {
        if (isEOF()) {
            return Collections.emptyList();
        }

        //TODO: these might need to be merged, so that in the end there is only one entry per ref context ?
        // these should NOT need to be sorted ??
        return getSlices().stream()
                .map(s -> s.getBAIEntries(compressorCache))
                .flatMap(List::stream)
                .sorted()
                .collect(Collectors.toList());
    }

    /**
     * Populate the indexing parameters of this Container's slices
     *
     * Requires: valid landmarks and containerBlocksByteSize
     *
     * @throws CRAMException when the Container is in an invalid state
     */
    private void distributeIndexingParametersToSlices() {
        if (slices.size() == 0) {
            return;
        }

        if (containerHeader.getLandmarks() == null) {
            throw new CRAMException("Cannot set Slice indexing parameters if this Container does not have landmarks");
        }

        if (containerHeader.getLandmarks().length != slices.size()) {
            final String format = "This Container's landmark and slice counts do not match: %d landmarks and %d slices";
            throw new CRAMException(String.format(format, containerHeader.getLandmarks().length, slices.size()));
        }

        if (containerHeader.getContainerBlocksByteSize() == 0) {
            throw new CRAMException("Cannot set Slice indexing parameters if the byte size of this Container's blocks is unknown");
        }

        final int lastSliceIndex = slices.size() - 1;
        for (int i = 0; i < lastSliceIndex; i++) {
            final Slice slice = slices.get(i);
            slice.setLandmarkIndex(i);
            slice.setByteOffsetOfSliceHeaderBlock(containerHeader.getLandmarks()[i]);
            slice.setByteSizeOfSliceBlocks(containerHeader.getLandmarks()[i + 1] - slice.getByteOffsetOfSliceHeaderBlock());
        }

        // get the last slice in the list, and
        final Slice lastSlice = slices.get(lastSliceIndex);
        lastSlice.setLandmarkIndex(lastSliceIndex);
        lastSlice.setByteOffsetOfSliceHeaderBlock(containerHeader.getLandmarks()[lastSliceIndex]);
        lastSlice.setByteSizeOfSliceBlocks(containerHeader.getContainerBlocksByteSize() - lastSlice.getByteOffsetOfSliceHeaderBlock());
    }

    private void checkReferenceContexts(final int aggregateReferenceContext) {
        final int actualReferenceContextID = getAlignmentContext().getReferenceContext().getReferenceContextID();
        if (actualReferenceContextID != aggregateReferenceContext) {
            throw new CRAMException(
                    String.format("actual container reference context (%d) doesn't match expected aggregate reference context (%d)",
                            actualReferenceContextID,
                            aggregateReferenceContext));
        }

        checkSliceReferenceContexts(actualReferenceContextID);
    }

    private void checkSliceReferenceContexts(final int actualReferenceContextID) {
        if (actualReferenceContextID == ReferenceContext.MULTIPLE_REFERENCE_ID) {
            for (final Slice slice : getSlices()) {
                if (slice.getAlignmentContext().getReferenceContext().getReferenceContextID() != ReferenceContext.MULTIPLE_REFERENCE_ID) {
                    throw new CRAMException(
                            String.format("Found slice with reference context (%d). Multi-reference container can only contain multi-ref slices.",
                                    slice.getAlignmentContext().getReferenceContext().getReferenceContextID()));
                }
            }
        }
    }

    // Compare the reference context declared by the container with the one derived from the contained slices.
    private void validateContainerReferenceContext() {
        final ReferenceContext derivedSliceReferenceContext = getDerivedReferenceContextFromSlices(getSlices());
        if (!derivedSliceReferenceContext.equals(getAlignmentContext().getReferenceContext())) {
            throw new CRAMException(String.format(
                    "Container (%s) has a reference context that doesn't match the reference context (%s) derived from it's slices.",
                    getContainerHeader(),toString(),
                    derivedSliceReferenceContext));
        }
    }

    private ReferenceContext getDerivedReferenceContextFromSlices(final List<Slice> containerSlices) {
        final Set<ReferenceContext> sliceRefContexts = containerSlices.stream()
                .map(s -> s.getAlignmentContext().getReferenceContext())
                .collect(Collectors.toSet());
        if (sliceRefContexts.isEmpty()) {
            throw new CRAMException("Cannot construct a container without any slices");
        }
        else if (sliceRefContexts.size() > 1) {
            return ReferenceContext.MULTIPLE_REFERENCE_CONTEXT;
        }

        return sliceRefContexts.iterator().next();
    }

    @Override
    public String toString() {
        return String.format("%s slices=%d",
                containerHeader.toString(),
                getSlices() == null ?
                        -1 :
                        getSlices().size());
    }

}
