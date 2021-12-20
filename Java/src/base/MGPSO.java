package base;

import util.FitnessFunction;
import util.Operator;

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
    protected final FitnessFunction[] f;
    protected final int[] n;

    public MGPSO(String[] seq,
                 int[] n,
                 double w,
                 double c1,
                 double c2,
                 double c3,
                 double vmax,
                 int vmaxiterlimit,
                 double[] term,
                 int maxIter,
                 FitnessFunction[] f,
                 Operator[] ops) throws Exception {
        super(seq, 0, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
        this.n = n;
        this.c3 = c3;
        this.f = f;

        if(n.length != f.length)
            throw new Exception("Amount of swarms must match amount of functions");
    }
}
