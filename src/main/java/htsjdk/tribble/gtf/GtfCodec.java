package htsjdk.tribble.gtf;

import htsjdk.samtools.util.*;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.*;
import htsjdk.tribble.util.ParsingUtils;

import java.io.*;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Predicate;
import java.util.zip.GZIPInputStream;

public class GtfCodec extends AbstractFeatureCodec<GtfFeature, LineIterator> {

    private final static Log logger = Log.getInstance(GtfCodec.class);

    private static final int NUM_FIELDS = 9;

    private static final int CONTIG_INDEX = 0;
    private static final int SOURCE_INDEX = 1;
    private static final int TYPE_INDEX = 2;
    private static final int START_INDEX = 3;
    private static final int END_INDEX = 4;
    private static final int SCORE_INDEX = 5;
    private static final int STRAND_INDEX = 6;
    private static final int FRAME_INDEX = 7;
    private static final int ATTRIBUTES_INDEX = 8;

    private final Queue<GtfFeatureImpl> activeFeatures = new ArrayDeque<>();
    private final Queue<GtfFeatureImpl> featuresToFlush = new ArrayDeque<>();

    private GtfFeatureImpl geneFeature;
    private boolean geneFeatureArtificial = false;
    private GtfFeatureImpl transcriptFeature;
    private boolean transcriptFeatureArtificial = false;
    private GtfFeatureImpl lastFeature;

    private final DecodeDepth decodeDepth;

    private int currentLine = 0;

    private final Predicate<String> filterOutAttribute;

    public GtfCodec() {
        this(DecodeDepth.DEEP);
    }

    public GtfCodec(final DecodeDepth decodeDepth) {
        this(decodeDepth, KEY -> false);
    }

    public GtfCodec(final DecodeDepth decodeDepth, final Predicate<String> filterOutAttribute) {
        super(GtfFeature.class);
        this.decodeDepth = decodeDepth;
        this.filterOutAttribute = filterOutAttribute;

        for (final String key : new String[] {GtfConstants.GENE_ID_ATTRIBUTE_KEY, GtfConstants.TRANSCRIPT_ID_ATTRIBUTE_KEY}) {
            if (filterOutAttribute.test(key)) {
                throw new IllegalArgumentException(String.format("Filtering out attribute <%s> is not allowed", key));
            }
        }
    }

    public enum DecodeDepth {
        DEEP,
        SHALLOW
    }

    @Override
    public GtfFeature decode(final LineIterator lineIterator) throws IOException {
        return decode(lineIterator, decodeDepth);
    }

    private GtfFeature decode(final LineIterator lineIterator, final DecodeDepth decodeDepth) throws IOException {
        currentLine++;

        if (!lineIterator.hasNext()) {
            // no more lines, flush whatever is active
            prepareToFlushFeatures();
            return featuresToFlush.poll();
        }

        final String line = lineIterator.next();

        if (line.startsWith(GtfConstants.COMMENT_START)) {
            return featuresToFlush.poll();
        }

        final GtfFeatureImpl thisFeature = new GtfFeatureImpl(parseLine(line, currentLine, this.filterOutAttribute));

        if (decodeDepth == DecodeDepth.DEEP) {
            // link to parents / children / co-features

            if (thisFeature.getType().equalsIgnoreCase("gene")) {
                // encountered a new gene, flush the previous gene
                prepareToFlushFeatures();

                activeFeatures.add(thisFeature);
                geneFeature = thisFeature;
                geneFeatureArtificial = false;
                transcriptFeature = null;
                lastFeature = null;

            } else if (thisFeature.getType().equalsIgnoreCase("transcript")) {
                // encountered a new transcript

                if (geneFeature == null) {
                    // encountered a transcript without a gene, generate a gene feature

                    prepareToFlushFeatures();

                    HashMap<String, List<String>> geneAttributes = new HashMap<>();
                    // copy gene attributes from the transcript to the gene attributes
                    for (Map.Entry<String, List<String>> entry : thisFeature.getAttributes().entrySet()) {
                        if (entry.getKey().toLowerCase().startsWith("gene")) {
                            geneAttributes.put(entry.getKey(), entry.getValue());
                        }
                    }

                    geneFeature = new GtfFeatureImpl(
                            new GtfBaseData(
                                    thisFeature.getContig(),
                                    thisFeature.getSource(),
                                    "gene",
                                    thisFeature.getStart(),
                                    thisFeature.getEnd(),
                                    thisFeature.getScore(),
                                    thisFeature.getStrand(),
                                    thisFeature.getFrame(),
                                    geneAttributes
                            )
                    );

                    activeFeatures.add(geneFeature);
                    geneFeatureArtificial = true;
                }

                if (geneFeatureArtificial) {
                    // update the genomic coordinates of the gene feature (useful when the gene line was not encountered)
                    if (thisFeature.getStart() < geneFeature.getStart()) {
                        geneFeature.getBaseData().setStart(thisFeature.getStart());
                    }
                    if (thisFeature.getEnd() > geneFeature.getEnd()) {
                        geneFeature.getBaseData().setEnd(thisFeature.getEnd());
                    }
                }

                activeFeatures.add(thisFeature);

                geneFeature.addChild(thisFeature);
                thisFeature.addParent(geneFeature);

                if (transcriptFeature != null) {
                    transcriptFeature.addCoFeature(thisFeature);
                }

                transcriptFeature = thisFeature;
                transcriptFeatureArtificial = false;
                lastFeature = null;

            } else {
                // any other feature requires that both a gene and a transcript have been encountered

                final String geneId = extractSingleAttribute(thisFeature.getAttribute(GtfConstants.GENE_ID_ATTRIBUTE_KEY));
                final String transcriptId = extractSingleAttribute(thisFeature.getAttribute(GtfConstants.TRANSCRIPT_ID_ATTRIBUTE_KEY));

                if (geneId == null) {
                    throw new TribbleException(String.format("Missing gene_id attribute in GTF line <%s>", line));
                }
                if (transcriptId == null) {
                    throw new TribbleException(String.format("Missing transcript_id attribute in GTF line <%s>", line));
                }

                if (geneFeature == null || !extractSingleAttribute(geneFeature.getAttribute(GtfConstants.GENE_ID_ATTRIBUTE_KEY))
                        .equalsIgnoreCase(geneId)) {

                    prepareToFlushFeatures();

                    transcriptFeature = null;
                    lastFeature = null;

                    HashMap<String, List<String>> geneAttributes = new HashMap<>();
                    // copy gene attributes from the transcript to the gene attributes
                    for (Map.Entry<String, List<String>> entry : thisFeature.getAttributes().entrySet()) {
                        if (entry.getKey().toLowerCase().startsWith("gene")) {
                            geneAttributes.put(entry.getKey(), entry.getValue());
                        }
                    }

                    geneFeature = new GtfFeatureImpl(
                            new GtfBaseData(
                                    thisFeature.getContig(),
                                    thisFeature.getSource(),
                                    "gene",
                                    thisFeature.getStart(),
                                    thisFeature.getEnd(),
                                    thisFeature.getScore(),
                                    thisFeature.getStrand(),
                                    thisFeature.getFrame(),
                                    geneAttributes
                            )
                    );

                    activeFeatures.add(geneFeature);

                    geneFeatureArtificial = true;
                }

                if (transcriptFeature == null || !extractSingleAttribute(transcriptFeature.getAttribute(GtfConstants.TRANSCRIPT_ID_ATTRIBUTE_KEY))
                        .equalsIgnoreCase(transcriptId)) {

                    HashMap<String, List<String>> transcriptAttributes = new HashMap<>();
                    // copy gene attributes from the transcript to the gene attributes
                    for (Map.Entry<String, List<String>> entry : thisFeature.getAttributes().entrySet()) {
                        if (entry.getKey().toLowerCase().startsWith("transcript")) {
                            transcriptAttributes.put(entry.getKey(), entry.getValue());
                        }
                    }

                    GtfFeatureImpl newTranscriptFeature = new GtfFeatureImpl(
                            new GtfBaseData(
                                    thisFeature.getContig(),
                                    thisFeature.getSource(),
                                    "transcript",
                                    thisFeature.getStart(),
                                    thisFeature.getEnd(),
                                    thisFeature.getScore(),
                                    thisFeature.getStrand(),
                                    thisFeature.getFrame(),
                                    transcriptAttributes
                            )
                    );

                    lastFeature = null;

                    activeFeatures.add(newTranscriptFeature);

                    geneFeature.addChild(newTranscriptFeature);
                    newTranscriptFeature.addParent(geneFeature);

                    if (transcriptFeature != null) {
                        transcriptFeature.addCoFeature(newTranscriptFeature);
                    }

                    transcriptFeature = newTranscriptFeature;

                    transcriptFeatureArtificial = true;
                }

                if (geneFeatureArtificial) {
                    // update the genomic coordinates of the gene feature (useful when the gene line was not encountered)
                    if (thisFeature.getStart() < geneFeature.getStart()) {
                        geneFeature.getBaseData().setStart(thisFeature.getStart());
                    }
                    if (thisFeature.getEnd() > geneFeature.getEnd()) {
                        geneFeature.getBaseData().setEnd(thisFeature.getEnd());
                    }
                }
                if (transcriptFeatureArtificial) {
                    // update the genomic coordinates of the transcript feature (useful when the gene line was not encountered)
                    if (thisFeature.getStart() < transcriptFeature.getStart()) {
                        transcriptFeature.getBaseData().setStart(thisFeature.getStart());
                    }
                    if (thisFeature.getEnd() > transcriptFeature.getEnd()) {
                        transcriptFeature.getBaseData().setEnd(thisFeature.getEnd());
                    }
                }

                activeFeatures.add(thisFeature);

                transcriptFeature.addChild(thisFeature);
                thisFeature.addParent(transcriptFeature);

                if (lastFeature != null) {
                    lastFeature.addCoFeature(thisFeature);
                }

                lastFeature = thisFeature;
            }

        } else if (decodeDepth == DecodeDepth.SHALLOW) {
            // submit every feature as soon as it is encountered without linking to parents / children / co-features
            prepareToFlushFeatures();
        }

        return featuresToFlush.poll();
    }

    private static Map<String, List<String>> parseAttributes(final String attributesString) throws UnsupportedEncodingException {

        if (attributesString.equals(GtfConstants.UNDEFINED_FIELD_VALUE)) {
            return Collections.emptyMap();
        }

        final Map<String, List<String>> attributes = new LinkedHashMap<>();

        final List<String> splitLine = ParsingUtils.split(attributesString, GtfConstants.ATTRIBUTE_DELIMITER);

        for (String attribute : splitLine) {

            attribute = attribute.trim();

            if (attribute.isEmpty()) {
                continue;
            }

            final List<String> keyValue = ParsingUtils.split(attribute.trim(), GtfConstants.KEY_VALUE_SEPARATOR);

            if (keyValue.size() < 2) {
                throw new TribbleException(String.format("Attribute <%s> is invalid", attribute));
            }

            final String key = URLDecoder.decode(keyValue.get(0).trim(), StandardCharsets.UTF_8);
            final String value = keyValue.subList(1, keyValue.size()).stream()
                    .reduce((a, b) -> a + b)
                    .orElse("")
                    .replaceAll("\"", "");

            attributes.put(key, decodeAttributeValue(value));
        }

        return attributes;
    }

    private static List<String> decodeAttributeValue(final String attributeValue) {
        final List<String> splitValues = ParsingUtils.split(attributeValue, GtfConstants.VALUE_DELIMITER);

        final List<String> decodedValues = new ArrayList<>();

        for (final String encodedValue : splitValues) {
            decodedValues.add(URLDecoder.decode(encodedValue.trim(), StandardCharsets.UTF_8));
        }

        return decodedValues;
    }

    static String extractSingleAttribute(final List<String> values) {
        if (values == null || values.isEmpty()) {
            return null;
        }

        if (values.size() != 1) {
            throw new TribbleException(String.format("Expected a single value, but found multiple values: <%s>", values));
        }

        return values.get(0);
    }

    private static GtfBaseData parseLine(final String line, final int currentLine, final Predicate<String> filterOutAttribute) {

        final List<String> splitLine = ParsingUtils.split(line, GtfConstants.FIELD_DELIMITER);

        if (splitLine.size() != NUM_FIELDS) {
            throw new TribbleException(String.format("GTF line <%s> does not have the expected number of fields", line));
        }

        try {
            final String contig = URLDecoder.decode(splitLine.get(CONTIG_INDEX), StandardCharsets.UTF_8);
            final String source = URLDecoder.decode(splitLine.get(SOURCE_INDEX), StandardCharsets.UTF_8);
            final String type = URLDecoder.decode(splitLine.get(TYPE_INDEX), StandardCharsets.UTF_8);
            final int start = Integer.parseInt(splitLine.get(START_INDEX));
            final int end = Integer.parseInt(splitLine.get(END_INDEX));
            final double score = splitLine.get(SCORE_INDEX).equals(GtfConstants.UNDEFINED_FIELD_VALUE) ? -1 : Double.parseDouble(splitLine.get(SCORE_INDEX));
            final int frame = splitLine.get(FRAME_INDEX).equals(GtfConstants.UNDEFINED_FIELD_VALUE) ? -1 : Integer.parseInt(splitLine.get(FRAME_INDEX));
            final Strand strand = Strand.decode(splitLine.get(STRAND_INDEX));
            final Map<String, List<String>> attributes = parseAttributes(splitLine.get(ATTRIBUTES_INDEX));

            // remove attributes matching the filter
            attributes.keySet().removeIf(filterOutAttribute);

            return new GtfBaseData(contig, source, type, start, end, score, strand, frame, attributes);

        } catch (NumberFormatException e) {
            throw new TribbleException(String.format("Error parsing GTF line <%s> at line %d", line, currentLine), e);
        } catch (IOException e) {
            throw new TribbleException(String.format("Error decoding GTF line <%s> at line %d", line, currentLine), e);
        }
    }

    @Override
    public Feature decodeLoc(LineIterator lineIterator) throws IOException {
        return decode(lineIterator, DecodeDepth.SHALLOW);
    }

    @Override
    public boolean canDecode(final String inputFilePath) {

        boolean canDecode;

        try {
            Path p = IOUtil.getPath(inputFilePath);
            canDecode = FileExtensions.GTF.stream().anyMatch(ext -> p.toString().endsWith(ext));

            if (canDecode) {

                final InputStream inputStream = IOUtil.hasGzipFileExtension(p) ?
                        new GZIPInputStream(Files.newInputStream(p)) :
                        Files.newInputStream(p);

                try (BufferedReader br = new BufferedReader(new InputStreamReader(inputStream))) {

                    String line = br.readLine();

                    while (line.startsWith(GtfConstants.COMMENT_START)) {
                        line = br.readLine();
                        if (line == null) {
                            return false;
                        }
                    }

                    final List<String> fields = ParsingUtils.split(line, GtfConstants.FIELD_DELIMITER);

                    canDecode &= fields.size() == NUM_FIELDS;

                    if (canDecode) {
                        try {
                            final int start = Integer.parseInt(fields.get(START_INDEX));
                            final int end = Integer.parseInt(fields.get(END_INDEX));
                        } catch (NumberFormatException e) {
                            return false;
                        }

                        if (!fields.get(ATTRIBUTES_INDEX).toLowerCase().contains("gene_id")) {
                            logger.error("GTF file does not contain gene_id attribute");
                            return false;
                        }
                        if (!fields.get(ATTRIBUTES_INDEX).toLowerCase().contains("transcript_id")) {
                            logger.error("GTF file does not contain transcript_id attribute");
                            return false;
                        }

                        final String strand = fields.get(STRAND_INDEX);

                        canDecode &= strand.equals(Strand.POSITIVE.toString()) ||
                                strand.equals(Strand.NEGATIVE.toString()) ||
                                strand.equals(Strand.NONE.toString()) ||
                                strand.equals("?");
                    }
                }
            }
        } catch (FileNotFoundException e) {
            logger.error(String.format("File not found: %s", inputFilePath));
            return false;
        } catch (IOException e) {
            logger.error(String.format("Error reading file: %s", inputFilePath));
            return false;
        }

        return canDecode;
    }

    @Override
    public FeatureCodecHeader readHeader(LineIterator lineIterator) {
        List<String> header = new ArrayList<>();

        while (lineIterator.hasNext()) {
            String line = lineIterator.peek();
            if (line.startsWith(GtfConstants.COMMENT_START)) {
                header.add(line);
                lineIterator.next();
            } else {
                break;
            }
        }

        return new FeatureCodecHeader(header, FeatureCodecHeader.NO_HEADER_END);
    }

    private void prepareToFlushFeatures() {
        featuresToFlush.addAll(activeFeatures);
        activeFeatures.clear();
    }

    @Override
    public LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(new SynchronousLineReader(bufferedInputStream));
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return new AsciiLineReaderIterator(AsciiLineReader.from(bufferedInputStream));
    }

    @Override
    public boolean isDone(final LineIterator lineIterator) {
        return !lineIterator.hasNext() && activeFeatures.isEmpty() && featuresToFlush.isEmpty();
    }

    @Override
    public void close(final LineIterator lineIterator) {
        featuresToFlush.clear();
        activeFeatures.clear();
        CloserUtil.close(lineIterator);
    }
}
