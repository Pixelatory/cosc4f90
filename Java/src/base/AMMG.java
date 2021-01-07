package base;

import MGPSO;
import util.Operator;

import java.util.ArrayList;

/*
    Nicholas Aksamit
    2020

    Angular Modulated Multi-Guided PSO for the MSA problem.
 */
class AMMG extends MGPSO {
    private int n;
    private int w;
    private double c1;
    private double c2;
    private double c3;
    private double vmax;
    private int vmaxiterlimit;
    private double term;
    private int maxIter;

    public AMMG(ArrayList<String> seq,
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
        super(seq, n, w, c1, c2, c3, vmax, vmaxiterlimit, term, maxIter, ops);
    }


}
