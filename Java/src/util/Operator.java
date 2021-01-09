package util;

/**
 * A comparison operator between two doubles.
 *
 * It will always return a boolean.
 */
public interface Operator {
    boolean operation(double a, double b);

    /**
     * Less than operator.
     * <p>
     * (a, b): true if a < b, else false.
     */
    Operator lt = (a, b) -> a < b;

    /**
     * Greater than operator.
     * <p>
     * (a, b): true if a > b, else false.
     */
    Operator gt = (a, b) -> a > b;

    /**
     * Equals operator.
     * <p>
     * (a, b): true if a = b, else false.
     */
    Operator eq = (a, b) -> a == b;
}
