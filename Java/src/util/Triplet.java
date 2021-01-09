package util;

public class Triplet<K, V, T> extends Pair<K, V> implements FixedCollection{
    private T t;

    public Triplet(K k, V v, T t) {
        super(k, v);
        this.t = t;
    }

    public T getThird() {
        return t;
    }

    public void setThird(T t) {
        this.t = t;
    }

    @Override
    public Pair<V, K> swap() {
        throw new UnsupportedOperationException("Triplet cannot swap values.");
    }
}
