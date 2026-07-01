public class BNFGrammar {

    // grammar symbols
    static final int NT_EXPR = 0;
    static final int NT_OP = 1;
    static final int NT_PREOP = 2;
    static final int NT_VAR = 3;
    static final int NT_CONST = 4;

    // grammar rules per non-terminal
    static final int EXPR_RULES = 4;
    static final int OP_RULES = 4;
    static final int PREOP_RULES = 2;

    static final int CONST_RULES = 4;

    static int numVars; // number of load input variables
    static int maxDepth; // maximum recursion depth during mapping
    static int maxWraps; // maximum chromosome wraps during mapping

    // constant table
    static final double[] CONST_VALUES = { 0.5, 1.0, 1.5, 2.0 };

    // Initialise grammar with run-time parameters
    public static void init(int numloadVars, int mappingMaxDepth, int mappingMaxWraps) {
        numVars = numloadVars;
        maxDepth = mappingMaxDepth;
        maxWraps = mappingMaxWraps;

    }

    // Maps the codon chromosome of an individual to a phenotype string.

    public static String mapGenotype(GEIndividual ind) {
        ind.resetMapping();

        String result = expandNT(NT_EXPR, 0, ind);

        if (result == null) {
            ind.phenotype = null;
            ind.valid = false;

        } else {
            ind.phenotype = result;
            ind.valid = true;

        }
        return result;
    }

    // Recursively expands non-terminal 'nt' at the current recursion depth.
    private static String expandNT(int nt, int depth, GEIndividual ind) {

        // EXPR
        if (nt == NT_EXPR) {

            int rule;
            if (depth >= maxDepth) {
                // Force a terminal rule to stop runaway recursion.
                int c = ind.nextCodon(maxWraps);
                if (c < 0)
                    return null;
                rule = 2 + (c % 2); // 2 or 3

            } else {
                int c = ind.nextCodon(maxWraps);
                if (c < 0)
                    return null;
                if (depth == 0)
                    rule = 0;
                else
                    rule = c % EXPR_RULES;

            }

            switch (rule) {
                case 0: {
                    // ( <expr> <op> <expr> )
                    String left = expandNT(NT_EXPR, depth + 1, ind);
                    if (left == null)
                        return null;
                    String op = expandNT(NT_OP, depth + 1, ind);
                    if (op == null)
                        return null;
                    String right = expandNT(NT_EXPR, depth + 1, ind);
                    if (right == null)
                        return null;
                    return "(" + left + " " + op + " " + right + ")";
                }
                case 1: {
                    // <pre-op> ( <expr> )
                    String preop = expandNT(NT_PREOP, depth + 1, ind);
                    if (preop == null)
                        return null;
                    String inner = expandNT(NT_EXPR, depth + 1, ind);
                    if (inner == null)
                        return null;
                    return preop + "(" + inner + ")";
                }
                case 2:
                    // <var>
                    return expandNT(NT_VAR, depth + 1, ind);
                case 3:
                    // <const>
                    return expandNT(NT_CONST, depth + 1, ind);
                default:
                    return null;
            }
        }

        // OP
        if (nt == NT_OP) {
            int c = ind.nextCodon(maxWraps);
            if (c < 0)
                return null;
            int rule = c % OP_RULES;
            String[] ops = { "+", "-", "*", "/" };
            return ops[rule];
        }

        // PRE-OP
        if (nt == NT_PREOP) {
            int c = ind.nextCodon(maxWraps);
            if (c < 0)
                return null;
            int rule = c % PREOP_RULES;
            String[] preops = { "abs", "neg" };
            return preops[rule];
        }

        // VAR
        if (nt == NT_VAR) {
            int c = ind.nextCodon(maxWraps);
            if (c < 0)
                return null;
            int rule = c % numVars;
            String varName = "x" + rule;
            return varName;
        }

        // CONST
        if (nt == NT_CONST) {
            int c = ind.nextCodon(maxWraps);
            if (c < 0)
                return null;
            int rule = c % CONST_RULES;
            double val = CONST_VALUES[rule];
            return String.valueOf(val);
        }

        return null; // unknown non-terminal
    }
}
