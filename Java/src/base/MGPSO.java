package base;

import base.PSO;
import util.FitnessFunction;
import util.Operator;
import util.Pair;

import java.util.ArrayList;

/*
    Nicholas Aksamit
    2020

    Basic Multi-Guided PSO Constructor
 */
public abstract class MGPSO extends PSO {
    protected final double c3;
    protected final Pair<FitnessFunction, Operator>[] f;

    public MGPSO(String[] seq,
                 int n,
                 double w,
                 double c1,
                 double c2,
                 double c3,
                 double vmax,
                 int vmaxiterlimit,
                 double[] term,
                 int maxIter,
                 Pair<FitnessFunction, Operator>[] f,
                 Operator[] ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
        this.c3 = c3;
        this.f = f;
    }
}
