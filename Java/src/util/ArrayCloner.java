package util;

import org.kamranzafar.commons.cloner.ObjectCloner;

public class ArrayCloner {
    public static <E> E[] deepcopy(E[] arr) {
        ObjectCloner<E[]> cloner = new ObjectCloner<>();
        return cloner.deepClone(arr);
    }

    public static <E> E[][] deepcopy(E[][] arr) {
        E[][] newArr = (E[][]) new Object[arr.length][];

        for (int i = 0; i < newArr.length; i++) {
            newArr[i] = ArrayCloner.deepcopy(arr[i]);
        }

        return newArr;
    }

    public static <E> E[][][] deepcopy(E[][][] arr) {
        E[][][] newArr = (E[][][]) new Object[arr.length][][];

        for (int i = 0; i < newArr.length; i++) {
            newArr[i] = ArrayCloner.deepcopy(arr[i]);
        }

        return newArr;
    }

}
