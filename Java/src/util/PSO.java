package util;

import java.util.ArrayList;
import java.util.Comparator;

/*
    Nicholas Aksamit
    2020

    Basic PSO Constructor
 */
public class PSO {
    protected ArrayList<String> seq;
    protected int n;
    protected double w;
    protected double c1;
    protected double c2;
    protected double vmax;
    protected int vmaxiterlimit;
    protected double term;
    protected int maxIter;
    protected ArrayList<Comparator<Double>> ops;

    public PSO(ArrayList<String> seq,
               int n,
               double w,
               double c1,
               double c2,
               double vmax,
               int vmaxiterlimit,
               double term,
               int maxIter,
               ArrayList<Comparator<Double>> ops) {
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