package util;

import java.util.ArrayList;

public interface FitnessFunction {
    double calculate(int[][] bitmatrix,
                     String[] seq);
}
