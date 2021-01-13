package base;

import base.PSO;
import util.FitnessFunction;
import util.Operator;
import util.Pair;
import util.Triplet;

import java.util.ArrayList;

/*
    Nicholas Aksamit
    2020

    Basic Multi-Guided PSO Constructor.

    FitnessFunction's are expected to be
    in minimization format. If maximization,
    then simply input -f(x).
 */
public abstract class MGPSO extends PSO {
    protected final double c3;
    protected final ArrayList<Triplet<FitnessFunction, Double, Double>> f;

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
                 FitnessFunction[] f,
                 Operator[] ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
        this.c3 = c3;
        this.f = new ArrayList<>();

        for(FitnessFunction ff : f) {
            this.f.add(new Triplet<>(ff, ff.calculateMax(seq), ff.calculateMin()));
        }
    }
}
