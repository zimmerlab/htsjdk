package htsjdk.tribble.gtf;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;

import java.util.*;

public class GtfBaseData {

    private final int hashCode;

    private final String contig;
    private final String source;
    private final String type;
    private final int start;
    private final int end;
    private final double score;
    private final Strand strand;
    private final int frame;
    private final Map<String, List<String>> attributes;

    public GtfBaseData(final String contig, final String source, final String type, final int start, final int end,
                       final double score, final Strand strand, final int frame,
                       final Map<String, List<String>> attributes) {
        this.contig = contig;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = Gff3BaseData.copyAttributesSafely(attributes);

        this.hashCode = computeHashCode();
    }

    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!other.getClass().equals(GtfBaseData.class)) {
            return false;
        }
        final GtfBaseData otherBaseData = (GtfBaseData) other;

        return otherBaseData.getContig().equals(getContig()) &&
                otherBaseData.getSource().equals(getSource()) &&
                otherBaseData.getType().equals(getType()) &&
                otherBaseData.getStart() == getStart() &&
                otherBaseData.getEnd() == getEnd() &&
                ((Double)otherBaseData.getScore()).equals(score) &&
                otherBaseData.getFrame() == getFrame() &&
                otherBaseData.getStrand().equals(getStrand()) &&
                otherBaseData.getAttributes().equals(getAttributes());
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    private int computeHashCode() {
        int hash = getContig().hashCode();
        hash = 31 * hash + getSource().hashCode();
        hash = 31 * hash + getType().hashCode();
        hash = 31 * hash + getStart();
        hash = 31 * hash + getEnd();
        hash = 31 * hash + Double.hashCode(getScore());
        hash = 31 * hash + getFrame();
        hash = 31 * hash + getStrand().hashCode();
        hash = 31 * hash + getAttributes().hashCode();

        return hash;
    }

    public String getContig() {
        return contig;
    }

    public String getSource() {
        return source;
    }

    public String getType() {
        return type;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public double getScore() {
        return score;
    }

    public Strand getStrand() {
        return strand;
    }

    public int getFrame() {
        return frame;
    }

    public Map<String, List<String>> getAttributes() {
        return attributes;
    }

    public List<String> getAttribute(final String key) {
        return attributes.getOrDefault(key, Collections.emptyList());
    }
}
