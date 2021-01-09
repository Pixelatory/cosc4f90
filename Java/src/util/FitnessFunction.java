package util;

import java.util.ArrayList;

public interface FitnessFunction {
    double calculate(ArrayList<ArrayList<Integer>> bitmatrix,
                     ArrayList<String> seq);
}
