package util;

import java.util.ArrayList;
import java.util.Comparator;

/*
    Nicholas Aksamit
    2020

    Basic Multi-Guided PSO Constructor
 */
public class MGPSO extends PSO {
    protected double c3;

    public MGPSO(ArrayList<String> seq,
                 int n,
                 double w,
                 double c1,
                 double c2,
                 double c3,
                 double vmax,
                 int vmaxiterlimit,
                 double term,
                 int maxIter,
                 ArrayList<Operator> ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
        this.c3 = c3;
    }
}
