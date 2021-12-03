package htsjdk.variant.vcf;

import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.IntGenotypeFieldAccessors;

import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;

/**
 * Functions specific to encoding VCF records.
 */
public class VCFEncoder {

    public static final Charset VCF_CHARSET = StandardCharsets.UTF_8;
    private static final String QUAL_FORMAT_STRING = "%.2f";
    private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

    private final IntGenotypeFieldAccessors GENOTYPE_FIELD_ACCESSORS = new IntGenotypeFieldAccessors();

    private VCFHeader header;

    private boolean allowMissingFieldsInHeader = false;

    private boolean outputTrailingFormatFields = false;

    private final VCFTextTransformer vcfTextTransformer;

    /**
     * Prepare a VCFEncoder that will encode records appropriate to the given VCF header, optionally
     * allowing missing fields in the header.
     */
    public VCFEncoder(final VCFHeader header, final boolean allowMissingFieldsInHeader, final boolean outputTrailingFormatFields) {
        if (header == null) {
            throw new NullPointerException("The VCF header must not be null.");
        }
        this.header = header;
        this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;
        this.outputTrailingFormatFields = outputTrailingFormatFields;
        this.vcfTextTransformer = header.getVCFHeaderVersion().isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_3)
            ? new VCFPercentEncodedTextTransformer()
            : new VCFPassThruTextTransformer();
    }

    /**
     * @deprecated since 10/24/13 use the constructor
     */
    @Deprecated
    public void setVCFHeader(final VCFHeader header) {
        this.header = header;
    }

    /**
     * @deprecated since 10/24/13 use the constructor
     */
    @Deprecated
    public void setAllowMissingFieldsInHeader(final boolean allow) {
        this.allowMissingFieldsInHeader = allow;
    }

    /**
     * encodes a {@link VariantContext} as a VCF line
     *
     * Depending on the use case it may be more efficient to {@link #write(Appendable, VariantContext)} directly
     * instead of creating an intermediate string.
     *
     * @return the VCF line
     */
    public String encode(final VariantContext context) {
        try {
            final StringBuilder stringBuilder = new StringBuilder(1000);
            write(stringBuilder, context);
            return stringBuilder.toString();
        } catch (final IOException error) {
            throw new RuntimeIOException("Cannot encode variant", error);
        }
    }


    /**
     * encodes a {@link VariantContext} context as VCF, and writes it directly to an {@link Appendable}
     *
     * This may be more efficient than calling {@link #encode(VariantContext)} and then writing the result since it
     * avoids creating an intermediate string.
     *
     * @param vcfOutput the {@link Appendable} to write to
     * @param context the variant
     * @return the java.lang.Appendable 'vcfOutput'
     * @throws IOException
     */
    public void write(final Appendable vcfOutput, VariantContext context) throws IOException {
        if (this.header == null) {
            throw new NullPointerException("The header field must be set on the VCFEncoder before encoding records.");
        }
        // If this context came from a version of VCF different from that of the header which we wrote
        // we need to decode the VC then re-encode its in-memory representation
        if (context.getVersion() != header.getVCFHeaderVersion()) {
            context = context.fullyDecode(header, true);
        }

        // CHROM
        vcfOutput.append(context.getContig()).append(VCFConstants.FIELD_SEPARATOR)
                // POS
                .append(String.valueOf(context.getStart())).append(VCFConstants.FIELD_SEPARATOR)
                // ID
                .append(context.getID()).append(VCFConstants.FIELD_SEPARATOR)
                // REF
                .append(context.getReference().getDisplayString()).append(VCFConstants.FIELD_SEPARATOR);

        // ALT
        if ( context.isVariant() ) {
            Allele altAllele = context.getAlternateAllele(0);
            String alt = altAllele.getDisplayString();
            vcfOutput.append(alt);

            for (int i = 1; i < context.getAlternateAlleles().size(); i++) {
                altAllele = context.getAlternateAllele(i);
                alt = altAllele.getDisplayString();
                vcfOutput.append(',');
                vcfOutput.append(alt);
            }
        } else {
            vcfOutput.append(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
        }

        vcfOutput.append(VCFConstants.FIELD_SEPARATOR);

        // QUAL
        if ( !context.hasLog10PError()) {
            vcfOutput.append(VCFConstants.MISSING_VALUE_v4);
        } else {
            vcfOutput.append(formatQualValue(context.getPhredScaledQual()));
        }
        vcfOutput.append(VCFConstants.FIELD_SEPARATOR)
                // FILTER
                .append(getFilterString(context)).append(VCFConstants.FIELD_SEPARATOR);

        // INFO
        final Map<String, String> infoFields = new TreeMap<>();
        for (final Map.Entry<String, Object> field : context.getAttributes().entrySet()) {
            if (!this.header.hasInfoLine(field.getKey())) {
                fieldIsMissingFromHeaderError(context, field.getKey(), "INFO");
            }

            final String outputValue = formatVCFField(field.getValue(), context.isFullyDecoded());
            if (outputValue != null) {
                infoFields.put(field.getKey(), outputValue);
            }
        }
        writeInfoString(infoFields, vcfOutput);

        // FORMAT
        final GenotypesContext gc = context.getGenotypes();
        if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
            vcfOutput.append(VCFConstants.FIELD_SEPARATOR);
            vcfOutput.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
        } else {
            final List<String> genotypeAttributeKeys = context.calcVCFGenotypeKeys(this.header);
            if ( !genotypeAttributeKeys.isEmpty()) {
                for (final String format : genotypeAttributeKeys) {
                    if (!this.header.hasFormatLine(format)) {
                        fieldIsMissingFromHeaderError(context, format, "FORMAT");
                    }
                }
                final String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);

                vcfOutput.append(VCFConstants.FIELD_SEPARATOR);
                vcfOutput.append(genotypeFormatString);

                final Map<Allele, String> alleleStrings = buildAlleleStrings(context);
                appendGenotypeData(context, alleleStrings, genotypeAttributeKeys, vcfOutput);
            }
        }
    }

    VCFHeader getVCFHeader() {
        return this.header;
    }

    boolean getAllowMissingFieldsInHeader() {
        return this.allowMissingFieldsInHeader;
    }

    private String getFilterString(final VariantContext vc) {
        if (vc.isFiltered()) {
            for (final String filter : vc.getFilters()) {
                if (!this.header.hasFilterLine(filter)) {
                    fieldIsMissingFromHeaderError(vc, filter, "FILTER");
                }
            }

            return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
        } else {
            return vc.filtersWereApplied() ? VCFConstants.PASSES_FILTERS_v4 : VCFConstants.UNFILTERED;
        }
    }

    private static String formatQualValue(final double qual) {
        String s = String.format(Locale.US, QUAL_FORMAT_STRING, qual);
        if (s.endsWith(QUAL_FORMAT_EXTENSION_TO_TRIM)) {
            s = s.substring(0, s.length() - QUAL_FORMAT_EXTENSION_TO_TRIM.length());
        }
        return s;
    }

    private void fieldIsMissingFromHeaderError(final VariantContext vc, final String id, final String field) {
        if (!allowMissingFieldsInHeader) {
            throw new IllegalStateException("Key " + id + " found in VariantContext field " + field
                    + " at " + vc.getContig() + ":" + vc.getStart()
                    + " but this key isn't defined in the VCFHeader.  We require all VCFs to have"
                    + " complete VCF headers by default.");
        }
    }

    String formatVCFField(final Object val, final boolean fullyDecoded) {
        if (val == null) {
            return VCFConstants.MISSING_VALUE_v4;
        } else if (val instanceof Double) {
            return formatVCFDouble((Double) val);
        } else if (val instanceof Boolean) {
            return (Boolean) val ? "" : null; // empty string for true, null for false
        } else if (val instanceof List) {
            return formatList((List<?>) val, fullyDecoded);
        } else if (val.getClass().isArray()) {
            return val.getClass().getComponentType().isPrimitive()
                ? formatPrimitiveArray(val)
                : formatList(Arrays.asList((Object[]) val), fullyDecoded);
        } else if (val instanceof String) {
            final String s = val.toString();
            // If the VariantContext from which this string was obtained was already fully decoded,
            // its in-memory representation may contain special characters which must be re-encoded,
            // while strings which have not been decoded yet represent the field as read directly
            // from the source VCF, so they are written back out without encoding
            return fullyDecoded ? vcfTextTransformer.encodeText(s) : s;
        } else {
            return val.toString();
        }
    }

    private String formatPrimitiveArray(final Object v) {
        final int len = Array.getLength(v);
        if (len == 0) return VCFConstants.MISSING_VALUE_v4;
        int i = 0;
        final StringBuilder s = new StringBuilder();
        if (v instanceof int[]) {
            final int[] a = (int[]) v;
            for (;;) {
                s.append(a[i++]);
                if (i == len) break;
                s.append(',');
            }
        } else if (v instanceof double[]) {
            final double[] a = (double[]) v;
            for (;;) {
                s.append(formatVCFDouble(a[i++]));
                if (i == len) break;
                s.append(',');
            }
        } else if (v instanceof long[]) {
            final long[] a = (long[]) v;
            for (;;) {
                s.append(a[i++]);
                if (i == len) break;
                s.append(',');
            }
        } else {
            for (;;) {
                s.append(formatVCFField(Array.get(v, i++), false));
                if (i == len) break;
                s.append(',');
            }
        }
        return s.toString();
    }

    private String formatList(final List<?> list, final boolean fullyDecoded) {
        if (list.isEmpty()) return VCFConstants.MISSING_VALUE_v4;
        final StringBuilder s = new StringBuilder();
        final Iterator<?> it = list.iterator();
        for (;;) {
            s.append(formatVCFField(it.next(), fullyDecoded));
            if (!it.hasNext()) break;
            s.append(',');
        }
        return s.toString();
    }

    private String formatGPKey(final Object val) {
        if (val instanceof List) return formatList((List<?>) val, true);
        if (!(val instanceof String)) {
            throw new TribbleException("Value for GP key was of unexpected type: " + val.getClass());
        }
        final String[] splits = ((String) val).split(",");
        // We need to special-case GP because there is a discrepancy in the scale used to record
        // its values between pre-4.3 and 4.3+ VCF. Pre-4.3 GP is phred scale encoded while
        // 4.3+ GP is a linear probability, bringing it in line with other standard keys that
        // use the P suffix (c.f. VCF 4.3 spec section 7.2).

        // Some tools in the wild apparently already use linear scaled GP, so we have to
        // be careful about converting inputs. We check whether GP values are already linear
        // scaled by seeing if the values' sum is approximately equal to 1, like we
        // would expect if the values were linear scale probabilities.
        // c.f. https://sourceforge.net/p/vcftools/mailman/vcftools-spec/thread/CEBCD558.FA29%25browning%40u.washington.edu/
        double sum = 0;

        final List<Double> rawGPValues = new ArrayList<>(splits.length);
        for (final String s : splits) {
            final double GP = VCFUtils.parseVcfDouble(s);
            rawGPValues.add(GP);
            sum += GP;
        }

        final boolean wasLinearScale = GeneralUtils.compareDoubles(sum, 1, VCFConstants.VCF_ENCODING_EPSILON) == 0;
        if (!wasLinearScale && header.getVCFHeaderVersion().isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_3)) {
            rawGPValues.replaceAll(GP -> QualityUtil.getErrorProbabilityFromPhredScore((int) Math.round(GP)));
        }
        return formatList(rawGPValues, true);
    }

    /**
     * Takes a double value and pretty prints it to a String for display
     * <p>
     * Large doubles =&gt; gets %.2f style formatting
     * Doubles &lt; 1 / 10 but &gt; 1/100 =&gt; get %.3f style formatting
     * Double &lt; 1/100 =&gt; %.3e formatting
     *
     * @param d
     * @return
     */
    public static String formatVCFDouble(final double d) {
        final String format;
        if (d < 1) {
            if (d < 0.01) {
                if (Math.abs(d) >= 1e-20) {
                    format = "%.3e";
                } else {
                    // return a zero format
                    return "0.00";
                }
            } else {
                format = "%.3f";
            }
        } else {
            format = "%.2f";
        }

        return String.format(Locale.US, format, d);
    }

    static int countOccurrences(final char c, final String s) {
        int count = 0;
        for (int i = 0; i < s.length(); i++) {
            count += s.charAt(i) == c ? 1 : 0;
        }
        return count;
    }

    static boolean isMissingValue(final String s) {
        // we need to deal with the case that it's a list of missing values
        return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
    }

    /*
     * Add the genotype data
     */
    public void addGenotypeData(final VariantContext vc, final Map<Allele, String> alleleMap, final List<String> genotypeFormatKeys, final StringBuilder builder) {
        try {
            appendGenotypeData(vc,alleleMap,genotypeFormatKeys,builder);
        } catch (final IOException err) {
            throw new RuntimeIOException("addGenotypeData failed",err);
        }
    }

    /**
     * Add the genotype Data to a java.lang.Appendable
     * @param vc the variant
     * @param alleleMap
     * @param genotypeFormatKeys
     * @param vcfoutput VCF output
     * @throws IOException
     */
    private void appendGenotypeData(final VariantContext vc, final Map<Allele, String> alleleMap, final List<String> genotypeFormatKeys, final Appendable vcfoutput) throws IOException {
        final int ploidy = vc.getMaxPloidy(2);

        for (final String sample : this.header.getGenotypeSamples()) {
            vcfoutput.append(VCFConstants.FIELD_SEPARATOR);

            Genotype g = vc.getGenotype(sample);
            if (g == null) {
                g = GenotypeBuilder.createMissing(sample, ploidy);
            }

            final List<String> attrs = new ArrayList<>(genotypeFormatKeys.size());
            for (final String field : genotypeFormatKeys) {
                if (field.equals(VCFConstants.GENOTYPE_KEY)) {
                    if (!g.isAvailable()) {
                        throw new IllegalStateException("GTs cannot be missing for some samples if they are available for others in the record");
                    }

                    writeGtField(alleleMap, vcfoutput, g);
                    continue;

                } else {
                    final String outputValue;
                    if (field.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
                        outputValue = g.isFiltered() ? g.getFilters() : VCFConstants.PASSES_FILTERS_v4;
                    } else if (field.equals(VCFConstants.GENOTYPE_POSTERIORS_KEY)) {
                        outputValue = g.hasExtendedAttribute(field)
                            ? formatGPKey(g.getExtendedAttribute(field))
                            : VCFConstants.MISSING_VALUE_v4;
                    } else {
                        final IntGenotypeFieldAccessors.Accessor accessor = GENOTYPE_FIELD_ACCESSORS.getAccessor(field);
                        if (accessor != null) {
                            final int[] intValues = accessor.getValues(g);
                            if (intValues == null) {
                                outputValue = VCFConstants.MISSING_VALUE_v4;
                            } else if (intValues.length == 1) { // fast path
                                outputValue = Integer.toString(intValues[0]);
                            } else {
                                final StringBuilder sb = new StringBuilder();
                                sb.append(intValues[0]);
                                for (int i = 1; i < intValues.length; i++) {
                                    sb.append(',');
                                    sb.append(intValues[i]);
                                }
                                outputValue = sb.toString();
                            }
                        } else {
                            Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : VCFConstants.MISSING_VALUE_v4;
                            outputValue = formatVCFField(val, vc.isFullyDecoded());
                        }
                    }

                    if (outputValue != null) {
                        attrs.add(outputValue);
                    }
                }
            }

            // strip off trailing missing values
            if (!outputTrailingFormatFields) {
                for (int i = attrs.size() - 1; i >= 0; i--) {
                    if (isMissingValue(attrs.get(i))) {
                        attrs.remove(i);
                    } else {
                        break;
                    }
                }
            }

            for (int i = 0; i < attrs.size(); i++) {
                if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY)) {
                    vcfoutput.append(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                }
                vcfoutput.append(attrs.get(i));
            }
        }
    }

    /**
     * write the encoded GT field for a Genotype
     * @param alleleMap a mapping of Allele -> GT allele value (from {@link this#buildAlleleStrings(VariantContext)}
     * @param vcfoutput the appendable to write to, to avoid inefficiency due to string copying
     * @param g the genotoype to encode
     * @throws IOException if appending fails with an IOException
     */
    public static void writeGtField(final Map<Allele, String> alleleMap, final Appendable vcfoutput, final Genotype g) throws IOException {
        writeAllele(g.getAllele(0), alleleMap, vcfoutput);
        for (int i = 1; i < g.getPloidy(); i++) {
            vcfoutput.append(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
            writeAllele(g.getAllele(i), alleleMap, vcfoutput);
        }
    }

    /*
     * Create the info string; assumes that no values are null
     */
    private void writeInfoString(final Map<String, String> infoFields, final Appendable vcfoutput) throws IOException {
        if ( infoFields.isEmpty() ) {
            vcfoutput.append(VCFConstants.EMPTY_INFO_FIELD);
            return;
        }

        boolean isFirst = true;
        for (final Map.Entry<String, String> entry : infoFields.entrySet()) {
            if (isFirst) {
                isFirst = false;
            } else {
                vcfoutput.append(VCFConstants.INFO_FIELD_SEPARATOR);
            }

            vcfoutput.append(entry.getKey());

            if ( ! entry.getValue().isEmpty()) {
                final VCFInfoHeaderLine metaData = this.header.getInfoHeaderLine(entry.getKey());
                if ( metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0 ) {
                    vcfoutput.append('=');
                    vcfoutput.append(entry.getValue());
                }
            }
        }
    }

    /**
     * Easy way to generate the GT field for a Genotype.  This will be less efficient than using
     * {@link this#writeGtField(Map, Appendable, Genotype)} because of redundant Map initializations
     * @param vc a VariantContext which must contain g or the results are likely to be incorrect
     * @param g a Genotype in vc
     * @return a String containing the encoding of the GT field of g
     */
    public static String encodeGtField(VariantContext vc, Genotype g) {
      final StringBuilder builder = new StringBuilder();
        try {
            writeGtField(VCFEncoder.buildAlleleStrings(vc), builder, g);
        } catch (final IOException e) {
            throw new RuntimeException("Somehow we failed to append to a StringBuilder, this shouldn't happen.", e);
        }
        return builder.toString();
    }

    /**
     * return a Map containing Allele -> String(allele position) for all Alleles in VC
     * (as well as NO_CALL)
     * ex: A,T,TC -> { A:0, T:1, TC:2, NO_CALL:EMPTY_ALLELE}
     * This may be efficient when looking up values for many genotypes per VC
     */
    public static Map<Allele, String> buildAlleleStrings(final VariantContext vc) {
        final Map<Allele, String> alleleMap = new HashMap<>(vc.getAlleles().size() + 1);
        alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for (int i = 0; i < alleles.size(); i++) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

    private static void writeAllele(final Allele allele, final Map<Allele, String> alleleMap, final Appendable vcfOutput) throws IOException {
        final String encoding = alleleMap.get(allele);
        if (encoding == null) {
            throw new RuntimeException("Allele " + allele + " is not an allele in the variant context");
        }
        vcfOutput.append(encoding);
    }
}
