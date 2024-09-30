package htsjdk.tribble.gtf;

import htsjdk.samtools.util.Tuple;
import htsjdk.tribble.annotation.Strand;

import java.util.*;
import java.util.stream.Collectors;

public class GtfFeatureImpl implements GtfFeature {

    private final GtfBaseData baseData;

    private GtfFeatureImpl parent;
    private final LinkedHashSet<GtfFeatureImpl> children = new LinkedHashSet<>();
    private final LinkedHashSet<GtfFeatureImpl> coFeatures = new LinkedHashSet<>();

    private final Set<GtfFeatureImpl> topLevelFeatures = new HashSet<>();

    public GtfFeatureImpl(final String contig, final String source, final String type,
                          final int start, final int end, final Double score, final Strand strand, final int frame,
                          final Map<String, List<String>> attributes) {
        this.baseData = new GtfBaseData(contig, source, type, start, end, score, strand, frame, attributes);
    }

    public GtfFeatureImpl(final GtfBaseData baseData) {
        this.baseData = baseData;
    }

    @Override
    public Set<GtfFeatureImpl> getTopLevelFeatures() {
        if (isTopLevelFeature()) {
            return Collections.singleton(this);
        }
        return topLevelFeatures;
    }

    @Override
    public boolean isTopLevelFeature() {
        return this.topLevelFeatures.isEmpty();
    }

    @Override
    public GtfFeatureImpl getParent() {
        return this.parent;
    }

    @Override
    public Set<GtfFeatureImpl> getChildren() {
        return this.children;
    }

    @Override
    public GtfBaseData getBaseData() {
        return this.baseData;
    }

    @Override
    public Set<GtfFeatureImpl> getAncestors() {
        final List<GtfFeatureImpl> ancestors = new ArrayList<>(Collections.singleton(this.parent));
        ancestors.addAll(parent.getAncestors());
        return new LinkedHashSet<>(ancestors);
    }

    @Override
    public Set<GtfFeatureImpl> getDescendents() {
        final List<GtfFeatureImpl> descendants = new ArrayList<>(children);
        for (final GtfFeatureImpl child : children) {
            descendants.addAll(child.getDescendents());
        }
        return new LinkedHashSet<>(descendants);
    }

    @Override
    public Set<GtfFeatureImpl> getCoFeatures() {
        return this.coFeatures;
    }

    @Override
    public boolean hasParent() {
        return this.parent != null;
    }

    @Override
    public boolean hasChildren() {
        return !this.children.isEmpty();
    }

    @Override
    public boolean hasCoFeatures() {
        return !this.coFeatures.isEmpty();
    }

    void addParent(final GtfFeatureImpl parent) {
        this.parent = parent;
        parent.children.add(this);

        addTopLevelFeatures(new HashSet<>(parent.getTopLevelFeatures()));
    }

    void addChild(final GtfFeatureImpl child) {
        children.add(child);
    }

    private void addTopLevelFeatures(final Collection<GtfFeatureImpl> topLevelFeaturesToAdd) {
        this.topLevelFeatures.addAll(topLevelFeaturesToAdd);

        for (final GtfFeatureImpl child : children) {
            child.addTopLevelFeatures(topLevelFeaturesToAdd);
            child.removeTopLevelFeatures(this);
        }
    }

    private void removeTopLevelFeatures(final GtfFeatureImpl topLevelFeatureToRemove) {
        this.topLevelFeatures.remove(topLevelFeatureToRemove);

        for (final GtfFeatureImpl child : children) {
            child.removeTopLevelFeatures(topLevelFeatureToRemove);
        }
    }

    public void addCoFeature(final GtfFeatureImpl coFeature) {
        if (!parent.equals(coFeature.parent)) {
            throw new IllegalArgumentException("Co-features must have the same parent");
        }
        for (final GtfFeatureImpl feature : coFeatures) {
            feature.addCoFeatureShallow(coFeature);
            coFeature.addCoFeatureShallow(feature);
        }
        addCoFeatureShallow(coFeature);
        coFeature.addCoFeatureShallow(this);
    }

    public void addCoFeatureShallow(final GtfFeatureImpl coFeature) {
        coFeatures.add(coFeature);
    }

    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof GtfFeature)) {
            return false;
        }
        /*
        To test for equality, the doubly linked list representation used to represent feature relationships is
        replaced with a graph representation. Equality for between two features is tested by testing equality
        between their base data fields, and equality between the graphs they are part of.
         */
        return baseData.equals(((GtfFeature) other).getBaseData()) &&
                new GtfGraph(this).equals(new GtfGraph((GtfFeature) other));
    }

    @Override
    public int hashCode() {
        return baseData.hashCode();
    }

    @Override
    public Set<GtfFeatureImpl> flatten() {
        final LinkedHashSet<GtfFeatureImpl> features = new LinkedHashSet<>(Collections.singleton(this));

        features.addAll(this.getDescendents());
        return features;
    }

    private static class GtfGraph {

        private final Set<GtfBaseData> nodes = new HashSet<>();
        private final Set<Tuple<GtfBaseData, GtfBaseData>> parentEdges = new HashSet<>();
        private final Set<Tuple<GtfBaseData, GtfBaseData>> childEdges = new HashSet<>();
        private final Set<Set<GtfBaseData>> coFeatureSets = new HashSet<>();

        GtfGraph(final GtfFeature feature) {
            feature.getTopLevelFeatures().stream().flatMap(f -> f.flatten().stream()).forEach(this::addFeature);
        }

        private void addFeature(final GtfFeature feature) {
            addNode(feature);
            addParentEdges(feature);
            addChildEdges(feature);
            addCoFeatureSet(feature);
        }

        private void addNode(final GtfFeature feature) {
            this.nodes.add(feature.getBaseData());
        }

        private void addParentEdges(final GtfFeature feature) {
            this.parentEdges.add(new Tuple<>(feature.getBaseData(), feature.getParent().getBaseData()));
        }

        private void addChildEdges(final GtfFeature feature) {
            for (final GtfFeature child : feature.getChildren()) {
                this.childEdges.add(new Tuple<>(feature.getBaseData(), child.getBaseData()));
            }
        }

        private void addCoFeatureSet(final GtfFeature feature) {
            if (feature.hasCoFeatures()) {
                final Set<GtfBaseData> coFeaturesBaseData = feature.getCoFeatures().stream().map(
                        GtfFeature::getBaseData).collect(Collectors.toSet());
                coFeaturesBaseData.add(feature.getBaseData());
                coFeatureSets.add(coFeaturesBaseData);
            }
        }

        @Override
        public boolean equals(Object other) {
            if (other == this) {
                return true;
            }
            if (!other.getClass().equals(GtfGraph.class)) {
                return false;
            }

            return nodes.equals(((GtfGraph) other).nodes) &&
                    parentEdges.equals(((GtfGraph) other).parentEdges) &&
                    childEdges.equals(((GtfGraph) other).childEdges) &&
                    coFeatureSets.equals(((GtfGraph) other).coFeatureSets);
        }

        @Override
        public int hashCode() {
            int hash = nodes.hashCode();
            hash = 31 * hash + parentEdges.hashCode();
            hash = 31 * hash + childEdges.hashCode();
            hash = 31 * hash + coFeatureSets.hashCode();

            return hash;
        }
    }
}
