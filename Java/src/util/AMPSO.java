package util;

import java.util.ArrayList;

public class AMPSO extends PSO {

    public AMPSO(ArrayList<String> seq, int n, double w, double c1, double c2, double vmax, int vmaxiterlimit, double term, int maxIter, ArrayList<Operator> ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
    }
}
