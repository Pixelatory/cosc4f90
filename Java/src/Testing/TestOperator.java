import com.rits.cloning.Cloner;
import org.junit.Test;
import util.Operator;

import java.util.ArrayList;

public class TestOperator {
    @Test
    public void testInfeasible() {
        Operator lt = (a, b) -> a < b;
        Operator gt = (a, b) -> a > b;

        // Testing less than variants
        assert(!lt.operation(5, 5));
        assert(lt.operation(5, 6));
        assert(!lt.operation(6, 5));

        assert(!gt.operation(5, 5));
        assert(gt.operation(6, 5));
        assert(!gt.operation(5, 6));
    }

    @Test
    public void testCloning() {
        ArrayList<ArrayList<String>> tmp = new ArrayList<>();
        tmp.add(new ArrayList<>());
        tmp.get(0).add("wow");

        Cloner cloner = new Cloner();
        ArrayList<ArrayList<String>> tmp2 = cloner.deepClone(tmp);
        System.out.println(Boolean.toString(tmp2 == tmp) + " " + Boolean.toString(tmp2.get(0) == tmp.get(0)) + " " + Boolean.toString(tmp2.get(0).get(0) == tmp.get(0).get(0)));

        ArrayList<String> tmp3 = new ArrayList<>();
        tmp3.add("wowwwww");

        ArrayList<String> tmp4 = cloner.deepClone(tmp3);
        System.out.println(Boolean.toString(tmp3 == tmp4) + " " + Boolean.toString(tmp3.get(0) == tmp4.get(0)));

        String tmp5 = "wowwwww";
        String tmp6 = new String(tmp5);
        System.out.println(tmp5 == tmp6);
    }
}
