package util;

public interface FitnessFunction {
    double calculate(int[][] bitmatrix,
                     String[] seq);

    double calculateMax(String[] seq);

    double calculateMin();
}
