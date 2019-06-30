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
package htsjdk.samtools.cram.encoding.core;

import htsjdk.samtools.cram.encoding.CRAMCodec;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.EncodingID;
import htsjdk.samtools.cram.structure.SliceBlocksReadStreams;
import htsjdk.samtools.cram.structure.SliceBlocksWriteStreams;

import java.nio.ByteBuffer;
import java.util.Arrays;

public class CanonicalHuffmanByteEncoding extends CRAMEncoding<Byte> {
    private final byte[] values;
    private final int[] bitLengths;
    private final ByteBuffer buf;

    public CanonicalHuffmanByteEncoding(final byte[] values, final int[] bitLengths) {
        super(EncodingID.HUFFMAN);
        this.values = values;
        this.bitLengths = bitLengths;
        this.buf = ByteBuffer.allocate(ITF8.MAX_BYTES * (values.length + bitLengths.length));
    }

    /**
     * Create a new instance of this encoding using the (ITF8 encoded) serializedParams.
     * @param serializedParams
     * @return CanonicalHuffmanByteEncoding with parameters populated from serializedParams
     */
    public static CanonicalHuffmanByteEncoding fromSerializedEncodingParams(final byte[] serializedParams) {
        final ByteBuffer buf = ByteBuffer.wrap(serializedParams);

        final int valueSize = ITF8.readUnsignedITF8(buf);
        final byte[] values = new byte[valueSize];
        buf.get(values);

        final int lengthSize = ITF8.readUnsignedITF8(buf);
        final int[] bitLengths = new int[lengthSize];
        for (int i = 0; i < lengthSize; i++) {
            bitLengths[i] = ITF8.readUnsignedITF8(buf);
        }

        return new CanonicalHuffmanByteEncoding(values, bitLengths);
    }

    @Override
    public byte[] toSerializedEncodingParams() {
        buf.clear();
        ITF8.writeUnsignedITF8(values.length, buf);
        for (final byte value : values) {
            buf.put(value);
        }

        ITF8.writeUnsignedITF8(bitLengths.length, buf);
        for (final int value : bitLengths) {
            ITF8.writeUnsignedITF8(value, buf);
        }

        buf.flip();
        final byte[] array = new byte[buf.limit()];
        buf.get(array);

        return array;
    }

    @Override
    public CRAMCodec<Byte> buildCodec(final SliceBlocksReadStreams sliceBlocksReadStreams, final SliceBlocksWriteStreams sliceBlocksWriteStreams) {
        return new CanonicalHuffmanByteCodec(
                sliceBlocksReadStreams == null ? null : sliceBlocksReadStreams.getCoreBlockInputStream(),
                sliceBlocksWriteStreams == null ? null : sliceBlocksWriteStreams.getCoreOutputStream(),
                values,
                bitLengths);
    }

    @Override
    public String toString() {
        return String.format("Values: %s BitLengths %s",
                Arrays.toString(values),
                Arrays.toString(bitLengths));
    }
}
