package htsjdk.samtools.cram.structure;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import htsjdk.samtools.cram.compression.ExternalCompressor;
import htsjdk.samtools.cram.compression.rans.RANS;
import htsjdk.samtools.cram.encoding.external.ByteArrayStopEncoding;
import htsjdk.samtools.cram.encoding.external.ExternalByteEncoding;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.InputStreamUtils;
import htsjdk.samtools.cram.structure.block.Block;
import htsjdk.samtools.cram.structure.block.BlockCompressionMethod;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * CRAM encoding map (EncodingParams for each Data Series, and compressors for each integer ContentID)
 * for a single CRAM container compression header.
 *
 * There are two constructors; one populates the map using the encodings chosen by this (htsjdk)
 * implementation, for use writing a CRAM, and one populates the map from a serialized
 * CompressionHeader stream when reading a CRAM, resulting in encodings chosen by the implementation
 * that wrote that CRAM.
 *
 * Although the CRAM spec defines a fixed list of data series, individual CRAM implementations
 * may choose to use only a subset of these. Therefore, the actual set of encodings that are
 * instantiated can vary depending on the source.
 *
 * Notes on the CRAM write implementation: This implementation encodes ALL DataSeries to external blocks,
 * (although some of the external encodings split the data between core and external; see
 * {@link htsjdk.samtools.cram.encoding.ByteArrayLenEncoding}, and does not use the 'BB' or 'QQ'
 * DataSeries when writing CRAM at all.
 *
 * See {@link htsjdk.samtools.cram.encoding.EncodingFactory} for details on how the encodings defined
 * here are mapped to the codecs that actually distribute the data to underlying blocks.
 */
public class CompressionHeaderEncodingMap {

    // encoding params for each data series
    private Map<DataSeries, EncodingParams> encodingMap = new TreeMap<>();

    // external compressor to use for each external block, keyed by external content ID; this
    // map contains a key for each data series that is in use, plus additional ones for blocks
    // used for tags
    private final Map<Integer, ExternalCompressor> externalCompressors = new TreeMap<>();

    /**
     * Constructor used to create an encoding map for writing CRAMs
     */
    public CompressionHeaderEncodingMap() {
        // NOTE: all of these encodings use external blocks and compressors for actual CRAM
        // data. The only use of the core encodings are is encoding params for other (external)
        // encodings, i.e., ByteArrayLenEncoding uses the core blocks to store the
        addExternalRansOrderZeroEncoding(DataSeries.AP_AlignmentPositionOffset);
        addExternalRansOrderOneEncoding(DataSeries.BA_Base);
        // the BB data series is not used by this implementation when writing CRAMs
        addExternalRansOrderOneEncoding(DataSeries.BF_BitFlags);
        addExternalGzipEncoding(DataSeries.BS_BaseSubstitutionCode);
        addExternalRansOrderOneEncoding(DataSeries.CF_CompressionBitFlags);
        addExternalGzipEncoding(DataSeries.DL_DeletionLength);
        addExternalGzipEncoding(DataSeries.FC_FeatureCode);
        addExternalGzipEncoding(DataSeries.FN_NumberOfReadFeatures);
        addExternalGzipEncoding(DataSeries.FP_FeaturePosition);
        addExternalGzipEncoding(DataSeries.HC_HardClip);
        addExternalByteArrayStopTabGzipEncoding(DataSeries.IN_Insertion);
        addExternalGzipEncoding(DataSeries.MF_MateBitFlags);
        addExternalGzipEncoding(DataSeries.MQ_MappingQualityScore);
        addExternalGzipEncoding(DataSeries.NF_RecordsToNextFragment);
        addExternalGzipEncoding(DataSeries.NP_NextFragmentAlignmentStart);
        addExternalRansOrderOneEncoding(DataSeries.NS_NextFragmentReferenceSequenceID);
        addExternalGzipEncoding(DataSeries.PD_padding);
        // the QQ data series is not used by this implementation when writing CRAMs
        addExternalRansOrderOneEncoding(DataSeries.QS_QualityScore);
        addExternalRansOrderOneEncoding(DataSeries.RG_ReadGroup);
        addExternalRansOrderZeroEncoding(DataSeries.RI_RefId);
        addExternalRansOrderOneEncoding(DataSeries.RL_ReadLength);
        addExternalByteArrayStopTabGzipEncoding(DataSeries.RN_ReadName);
        addExternalGzipEncoding(DataSeries.RS_RefSkip);
        addExternalByteArrayStopTabGzipEncoding(DataSeries.SC_SoftClip);
        addExternalGzipEncoding(DataSeries.TC_TagCount);
        addExternalGzipEncoding(DataSeries.TL_TagIdList);
        addExternalGzipEncoding(DataSeries.TN_TagNameAndType);
        addExternalRansOrderOneEncoding(DataSeries.TS_InsertSize);
    }

    /**
     * Constructor used to create an encoding map when reading a CRAM.
     * @param inputStream the CRAM input stream to be consumed
     */
    public CompressionHeaderEncodingMap(final InputStream inputStream) {
        final int byteSize = ITF8.readUnsignedITF8(inputStream);
        final byte[] bytes = new byte[byteSize];
        InputStreamUtils.readFully(inputStream, bytes, 0, bytes.length);
        final ByteBuffer buffer = ByteBuffer.wrap(bytes);

        final int mapSize = ITF8.readUnsignedITF8(buffer);

        for (int i = 0; i < mapSize; i++) {
            final String dataSeriesAbbreviation = new String(new byte[]{buffer.get(), buffer.get()});
            final DataSeries dataSeries = DataSeries.byCanonicalName(dataSeriesAbbreviation);

            final EncodingID id = EncodingID.values()[buffer.get()];
            final int paramLen = ITF8.readUnsignedITF8(buffer);
            final byte[] paramBytes = new byte[paramLen];
            buffer.get(paramBytes);

            //TODO: why can't this just instantiate and store the CRAMEncoding directly, and get rid of EncodingParams
            encodingMap.put(dataSeries, new EncodingParams(id, paramBytes));
        }
    }

    /**
     * Add an external compressor for a tag block
     * @param tagId the tag as a content ID
     * @param compressor compressor to be used for this tag block
     */
    public void addTagBlockCompression(final int tagId, final ExternalCompressor compressor) {
        externalCompressors.put(tagId, compressor);
    }

    /**
     * Get the encodong params that should be used for a given DataSeries.
     * @param dataSeries
     * @return EncodingParams for the DataSeries
     */
    public EncodingParams getEncodingParamsForDataSeries(final DataSeries dataSeries) {
        return encodingMap.get(dataSeries);
    }

    /**
     * Get a list of all external IDs for this encoding map
     * @return list of all external IDs for this encoding map
     */
    public List<Integer> getExternalIDs() { return new ArrayList(externalCompressors.keySet()); }

    /**
     * Given a content ID, return a {@link Block} for that ID by obtaining the contents of the stream,
     * compressing it using the compressor for that contentID, and converting to a {@link Block}.
     * @param contentId contentID to use
     * @param outputStream stream to compress
     * @return Block containing the compressed contends of the stream
     */
    public Block createCompressedBlockForStream(final Integer contentId, final ByteArrayOutputStream outputStream) {
        final ExternalCompressor compressor = externalCompressors.get(contentId);
        final byte[] rawContent = outputStream.toByteArray();
        return Block.createExternalBlock(
                compressor.getMethod(),
                contentId,
                compressor.compress(rawContent),
                rawContent.length);
    }

    /**
     * Write the encoding map out to a CRAM Stream
     * @param outputStream
     * @throws IOException
     */
    public void write(final OutputStream outputStream) throws IOException {
        // encoding map:
        int size = 0;
        for (final DataSeries dataSeries : encodingMap.keySet()) {
            // not all DataSeries are used by this implementation
            if (encodingMap.get(dataSeries).id != EncodingID.NULL) {
                size++;
            }
        }

        final ByteBuffer mapBuffer = ByteBuffer.allocate(1024 * 100);
        ITF8.writeUnsignedITF8(size, mapBuffer);
        for (final DataSeries dataSeries : encodingMap.keySet()) {
            if (encodingMap.get(dataSeries).id == EncodingID.NULL) {
                // not all DataSeries are used by this implementation
                continue;
            }

            final String dataSeriesAbbreviation = dataSeries.getCanonicalName();
            mapBuffer.put((byte) dataSeriesAbbreviation.charAt(0));
            mapBuffer.put((byte) dataSeriesAbbreviation.charAt(1));

            final EncodingParams params = encodingMap.get(dataSeries);
            mapBuffer.put((byte) (0xFF & params.id.getId()));
            ITF8.writeUnsignedITF8(params.params.length, mapBuffer);
            mapBuffer.put(params.params);
        }
        mapBuffer.flip();
        final byte[] mapBytes = new byte[mapBuffer.limit()];
        mapBuffer.get(mapBytes);

        ITF8.writeUnsignedITF8(mapBytes.length, outputStream);
        outputStream.write(mapBytes);
    }

    private void addExternalEncoding(final DataSeries dataSeries,
                                     final EncodingParams params,
                                     final ExternalCompressor compressor) {
        externalCompressors.put(dataSeries.getExternalBlockContentId(), compressor);
        encodingMap.put(dataSeries, params);
    }

    private void addExternalByteArrayStopTabGzipEncoding(final DataSeries dataSeries) {
        addExternalEncoding(dataSeries,
                new ByteArrayStopEncoding((byte) '\t', dataSeries.getExternalBlockContentId()).toParam(),
                ExternalCompressor.createGZIP());
    }

    // Visible for testing, because without this we have no way to unit test round tripping an
    // encoding map that contains the handful of data series that htsjdk generally doesn't use
    // when writing, since there is no code add those data series to the map.
    void addExternalEncoding(final DataSeries dataSeries, final ExternalCompressor compressor) {
        // we need a concrete type; the choice of Byte is arbitrary.
        // params are equal for all External Encoding value types
        final EncodingParams params = new ExternalByteEncoding(dataSeries.getExternalBlockContentId()).toParam();
        addExternalEncoding(dataSeries, params, compressor);
    }

    private void addExternalGzipEncoding(final DataSeries dataSeries) {
        addExternalEncoding(dataSeries, ExternalCompressor.createGZIP());
    }

    private void addExternalRansOrderOneEncoding(final DataSeries dataSeries) {
        addExternalEncoding(dataSeries, ExternalCompressor.createRANS(RANS.ORDER.ONE));
    }

    private void addExternalRansOrderZeroEncoding(final DataSeries dataSeries) {
        addExternalEncoding(dataSeries, ExternalCompressor.createRANS(RANS.ORDER.ZERO));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CompressionHeaderEncodingMap that = (CompressionHeaderEncodingMap) o;

        if (!this.encodingMap.equals(that.encodingMap)) return false;
        return this.externalCompressors.equals(that.externalCompressors);
    }

    @Override
    public int hashCode() {
        int result = encodingMap.hashCode();
        result = 31 * result + externalCompressors.hashCode();
        return result;
    }

    private static class SerializedEncodingMap {
        final List<EncodingMapEntry> entries = new ArrayList<>();

        public SerializedEncodingMap(final CompressionHeaderEncodingMap targetMap) {
            targetMap.encodingMap.entrySet().forEach(
                    e -> entries.add(
                                new EncodingMapEntry(
                                        e.getKey().getExternalBlockContentId(),
                                        e.getKey(),
                                        e.getValue(),
                                        targetMap.externalCompressors.get(e.getKey().getExternalBlockContentId())))
                    );
        }
    }

    private static class EncodingMapEntry {
        final Integer contentID;
        final DataSeries dataSeries;
        final EncodingParams encodingParams;
        final BlockCompressionMethod compressionMethod;
        final RANS.ORDER order;

        public EncodingMapEntry(final Integer contentID,
                                final DataSeries dataSeries,
                                final EncodingParams encodingParams,
                                final ExternalCompressor compressor){
            this.contentID = contentID;
            this.dataSeries = dataSeries;
            this.encodingParams = encodingParams;
            this.compressionMethod = compressor.getMethod();
            //TODO: fix this
            if (compressionMethod == BlockCompressionMethod.RANS) {
                this.order = compressor.getClass().equals(ExternalCompressor.RANSOrder0ExternalCompressor.class) ?
                        RANS.ORDER.ZERO :
                        RANS.ORDER.ONE;
            } else {
                this.order = null;
            }
        }
    }

    public CompressionHeaderEncodingMap(final SerializedEncodingMap serializedMap) {
        serializedMap.entries.forEach(e -> {
            this.encodingMap.put(e.dataSeries, e.encodingParams);
            this.externalCompressors.put(e.contentID, ExternalCompressor.getCompressorForMethod(e.compressionMethod, e.order));
        });
    }

    /**
     * Constructor used to create an encoding map from a serialized JSON file when writing a CRAM.
     * @param encodingMapPath the CRAM input stream to be consumed
     */
    public static CompressionHeaderEncodingMap readFromPath(final Path encodingMapPath) {
        try (final BufferedReader fr = Files.newBufferedReader(encodingMapPath)) {
            final Gson gson = new Gson();
            final SerializedEncodingMap serializedMap = gson.fromJson(fr, SerializedEncodingMap.class);
            final CompressionHeaderEncodingMap encodingMap = new CompressionHeaderEncodingMap(serializedMap);
            return encodingMap;
        } catch (final IOException e) {
            throw new RuntimeIOException("Failed reading encoding map json file", e);
        }
    }

    /**
     * Constructor used to write a serialized (JSON) encoding map to record the encoding map
     * used to create a particular CRAM stream..
     * @param encodingMapPath the path to an encoding map JSON file
     */
    public void writeToPath(final Path encodingMapPath) {
        //TODO: add a version #
        try (final BufferedWriter fileWriter = Files.newBufferedWriter(encodingMapPath)) {
            final GsonBuilder gsonBuilder = new GsonBuilder().enableComplexMapKeySerialization().setPrettyPrinting();
            final SerializedEncodingMap serializedMap = new SerializedEncodingMap(this);
            final Gson gson = gsonBuilder.create();
            final String jsonEncodingString = gson.toJson(serializedMap, SerializedEncodingMap.class);
            fileWriter.write(jsonEncodingString);
        } catch (final IOException e) {
            throw new RuntimeIOException("Failed writing json file for encoding map", e);
        }
    }
}
