import util.Helper;
import util.Pair;

import java.util.concurrent.ThreadLocalRandom;

public class Testing {
    public static boolean test(double x, double y, double z) {
        return x + y < (24 * (1 - Math.pow(z, 2))) / (7 - 5 * z);
    }

    public static boolean test2(double c1, double c2, double c3, double l, double w) {
        double val = c1 + l * c2 + (1 - l) * c3;
        return 0 < val &&
                val < (4 * (1 - Math.pow(w, 2))) / (1 - w + (((Math.pow(c1, 2) + Math.pow(l, 2) * Math.pow(c2, 2) + Math.pow(1 - l, 2) * Math.pow(c3, 2)) * (1 + w)) / (3 * Math.pow(val, 2))));
    }

    public static void main(String[] args) {
        //System.out.println(test(1.496180, 1.496180, 0.729844));

        /*for (int i = 0; i < 10000; i++) {
            double c1 = ThreadLocalRandom.current().nextDouble(0, 2);
            double c2 = ThreadLocalRandom.current().nextDouble(0, 2);
            double w = ThreadLocalRandom.current().nextDouble(0, 1);
            if (test(c1, c2, w) && w >= 0.1)
                System.out.println(c1 + " " + c2 + " " + w);
        }*/

        /*for (int i = 0; i < 10000; i++) {
            double c1 = ThreadLocalRandom.current().nextDouble(0, 2);
            double c2 = ThreadLocalRandom.current().nextDouble(0, 2);
            double c3 = ThreadLocalRandom.current().nextDouble(0, 2);
            double w = ThreadLocalRandom.current().nextDouble(0, 1);
            boolean tmp = true;
            for (double l = 0; l <= 1; l+= 0.01) {
                if (!test2(c1, c2, c3, l, w)) {
                    tmp = false;
                    break;
                }
            }
            if (tmp)
                System.out.println(c1 + " " + c2 + " " + c3 + " " + w);
        }*/
        Pair<Integer, Integer>[] p = new Pair[]{new Pair<>(-25, 5), new Pair<>(-15, 10)};
        Pair<Integer, Integer>[] p2 = new Pair[]{new Pair<>(-20, 5), new Pair<>(-10, 10)};
        System.out.println(Helper.hypervolumeIndicator(p, new Pair<>(-50, 20)));
        System.out.println(Helper.hypervolumeIndicator(p2, new Pair<>(-50, 20)));
    }
}
