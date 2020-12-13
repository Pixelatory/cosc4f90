import util.PSO;

import java.util.ArrayList;
import java.util.Comparator;

/*
    Nicholas Aksamit
    2020

    Angular Modulated PSO for the MSA problem.
 */
class AM extends PSO {

    public AM(ArrayList<String> seq,
              int n,
              double w,
              double c1,
              double c2,
              double vmax,
              int vmaxiterlimit,
              double term,
              int maxIter,
              ArrayList<Comparator<Double>> ops) {
        super(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, ops);
    }
}
