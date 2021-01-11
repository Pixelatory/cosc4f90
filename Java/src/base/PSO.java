package base;

import util.Operator;

import java.util.ArrayList;
import java.util.Comparator;

/*
    Nicholas Aksamit
    2020

    Basic PSO Constructor
 */
public abstract class PSO extends Thread {
    protected final String[] seq;
    protected final int n;
    protected final double w;
    protected final double c1;
    protected final double c2;
    protected final double vmax;
    protected final int vmaxiterlimit;
    protected final double[] term;
    protected final int maxIter;
    protected final Operator[] ops;

    public PSO(String[] seq,
               int n,
               double w,
               double c1,
               double c2,
               double vmax,
               int vmaxiterlimit,
               double[] term,
               int maxIter,
               Operator[] ops) {
        this.seq = seq;
        this.n = n;
        this.w = w;
        this.c1 = c1;
        this.c2 = c2;
        this.vmax = vmax;
        this.vmaxiterlimit = vmaxiterlimit;
        this.term = term;
        this.maxIter = maxIter;
        this.ops = ops;
    }
}
