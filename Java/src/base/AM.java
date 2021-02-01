package base;

import util.Operator;

/*
    Nicholas Aksamit
    2020

    Angular Modulated PSO for the MSA problem.
 */
public abstract class AM extends PSO {

    public AM(String[] seq,
              int n,
              double w,
              double c1,
              double c2,
              double vmax,
              int vmaxiterlimit,
              double[] term,
              int maxIter,
              Operator[] ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
    }
}
