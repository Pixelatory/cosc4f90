public class Testing2 {
    public static void main(String[] args) {
        int lower_bound = 0;
        int upper_bound = 30;
        int num = 0;
        num -= 1;
        System.out.println((num - lower_bound) % (upper_bound-lower_bound) + lower_bound);

        int x = 30;
        if (x < lower_bound)
            x = upper_bound - (lower_bound - x) % (upper_bound - lower_bound);
        else
            x = lower_bound + (x - lower_bound) % (upper_bound - lower_bound);
        System.out.println(x);
    }
}
