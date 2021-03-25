import java.util.concurrent.ThreadLocalRandom;

public class Testing {
    public static boolean test(double x, double y, double z) {
        return x + y < (24 * (1 - Math.pow(z, 2))) / (7 - 5 * z);
    }

    public static void main(String[] args) {
        System.out.println(test(1.496180, 1.496180, 0.729844));

        for (int i = 0; i < 10000; i++) {
            double c1 = ThreadLocalRandom.current().nextDouble(0, 2);
            double c2 = ThreadLocalRandom.current().nextDouble(0, 2);
            double w = ThreadLocalRandom.current().nextDouble(0, 1);
            if (test(c1, c2, w) && w >= 0.1)
                System.out.println(c1 + " " + c2 + " " + w);
        }
    }
}
