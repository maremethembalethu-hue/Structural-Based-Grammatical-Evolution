
public class ExpressionEvaluator {

    // shared mutable state during one parse call
    private String expr;
    private int pos;
    private double[] inputs;

    // Evaluates the expression string with the given input variable values.
    public double evaluate(String expression, double[] inputValues) {
        if (expression == null || expression.isEmpty())
            return 0.0;
        this.expr = expression.trim();
        this.pos = 0;
        this.inputs = inputValues;
        try {
            double result = parseExpr();
            if (Double.isNaN(result) || Double.isInfinite(result)
                    || Math.abs(result) > 1e6) {
                return 0.0;
            }
            return result;
        } catch (Exception e) {
            return 0.0; // any parse error to treat as invalid
        }
    }

    // Recursive-descent parser
    private double parseExpr() {
        double result = parseTerm();
        while (pos < expr.length()) {
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == '+') {
                pos++;
                result += parseTerm();
            } else if (pos < expr.length() && expr.charAt(pos) == '-') {
                pos++;
                result -= parseTerm();
            } else {
                break;
            }
        }
        return result;
    }

    private double parseTerm() {
        double result = parseFactor();
        while (pos < expr.length()) {
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == '*') {
                pos++;
                result *= parseFactor();
            } else if (pos < expr.length() && expr.charAt(pos) == '/') {
                pos++;
                double denom = parseFactor();
                // Protected division
                result = (Math.abs(denom) > 1e-9) ? result / denom : 1.0;
            } else {
                break;
            }
        }
        return result;
    }

    private double parseFactor() {
        skipWhitespace();

        if (pos < expr.length() && expr.charAt(pos) == '(') {
            pos++;
            double val = parseExpr();
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == ')')
                pos++;
            return val;
        }

        if (expr.startsWith("abs", pos)) {
            pos += 3;
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == '(')
                pos++;
            double val = parseExpr();
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == ')')
                pos++;
            return Math.abs(val);
        }

        if (expr.startsWith("neg", pos)) {
            pos += 3;
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == '(')
                pos++;
            double val = parseExpr();
            skipWhitespace();
            if (pos < expr.length() && expr.charAt(pos) == ')')
                pos++;
            return -val;
        }

        if (pos < expr.length() && expr.charAt(pos) == 'x') {
            int start = pos;
            pos++; // skip 'x'
            while (pos < expr.length() && Character.isDigit(expr.charAt(pos)))
                pos++;
            int idx = Integer.parseInt(expr.substring(start + 1, pos));
            if (idx < inputs.length)
                return inputs[idx];
            return 0.0;
        }

        // number
        return parseNumber();
    }

    // Parse a numeric
    private double parseNumber() {
        skipWhitespace();
        int start = pos;
        if (pos < expr.length() && expr.charAt(pos) == '-')
            pos++; // leading minus
        while (pos < expr.length()
                && (Character.isDigit(expr.charAt(pos)) || expr.charAt(pos) == '.')) {
            pos++;
        }
        if (start == pos)
            return 0.0;
        try {
            return Double.parseDouble(expr.substring(start, pos));
        } catch (NumberFormatException e) {
            return 0.0;
        }
    }

    // Skip any whitespace at the current position.
    private void skipWhitespace() {
        while (pos < expr.length() && Character.isWhitespace(expr.charAt(pos)))
            pos++;
    }
}
