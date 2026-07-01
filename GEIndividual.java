import java.util.Arrays;

public class GEIndividual {

    int[] codons; // the chromosome integers
    int codonLength; // number of codons

    int codonIndex; // next codon to consume
    int wraps; // how many times we have wrapped around
    String phenotype; // expression string
    boolean valid; // false
    double fitness;

    public GEIndividual(int[] codons) {
        this.codons = Arrays.copyOf(codons, codons.length);
        this.codonLength = codons.length;
        this.phenotype = null;
        this.valid = false;
        this.fitness = 1e9;
        this.codonIndex = 0;
        this.wraps = 0;
    }

    // Deep copy
    public GEIndividual copy() {
        GEIndividual clone = new GEIndividual(Arrays.copyOf(codons, codonLength));
        clone.phenotype = this.phenotype;
        clone.valid = this.valid;
        clone.fitness = this.fitness;
        return clone;
    }

    // Consume next codon, with wrapping
    public int nextCodon(int maxWraps) {
        if (codonIndex >= codonLength) {
            wraps++;
            codonIndex = 0;
            if (wraps > maxWraps)
                return -1; // exceeded wrap limit
        }
        return codons[codonIndex++];
    }

    // Reset mapping state
    public void resetMapping() {
        codonIndex = 0;
        wraps = 0;
        phenotype = null;
        valid = false;
    }

    // toString – for tracing / debug
    @Override
    public String toString() {
        return "GEIndividual{"
                + "fitness=" + String.format("%.6f", fitness)
                + ", valid=" + valid
                + ", phenotype=" + phenotype
                + ", codons=" + Arrays.toString(codons)
                + "}";
    }
}
