package htsjdk.tribble.gtf;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

import java.util.List;
import java.util.Set;

public interface GtfFeature extends Feature {

    GtfBaseData getBaseData();

    @Override
    default String getContig() {
        return getBaseData().getContig();
    }

    default String getSource() {
        return getBaseData().getSource();
    }

    default String getType() {
        return getBaseData().getType();
    }

    @Override
    default int getStart() {
        return getBaseData().getStart();
    }

    @Override
    default int getEnd() {
        return getBaseData().getEnd();
    }

    default double getScore() {
        return getBaseData().getScore();
    }

    default Strand getStrand() {
        return getBaseData().getStrand();
    }

    default int getFrame() {
        return getBaseData().getFrame();
    }

    default List<String> getAttribute(final String key) {
        return getBaseData().getAttribute(key);
    }

    /**
     * Get the set of top level features from which this feature is descended.
     * Top level features are features with no linked parents
     * @return set of top level feature from which this feature is descended
     */
    Set<? extends GtfFeature> getTopLevelFeatures();

    boolean isTopLevelFeature();

    GtfFeature getParent();

    Set<? extends GtfFeature> getChildren();

    Set<? extends GtfFeature> getAncestors();

    Set<? extends GtfFeature> getDescendents();

    Set<? extends GtfFeature> getCoFeatures();

    Set<? extends GtfFeature> flatten();

    boolean hasParent();

    boolean hasChildren();

    boolean hasCoFeatures();
}
