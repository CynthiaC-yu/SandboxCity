import java.util.Random;

/** Generating and save a set of transit level and its thresholds */
public class TransitLevels {
    public static class Level {
        public final String name;
        public double lower;   // inclusive
        public double upper;   // exclusive (Double.POSITIVE_INFINITY 表示无上界)

        Level(String name, double lower, double upper) {
            this.name  = name;
            this.lower = lower;
            this.upper = upper;
        }
        /** Determine which cate does the original transit value belongs to */
        public boolean contains(double v) {
            return v >= lower && v < upper;
        }
        @Override public String toString() {
            String lo = String.format("%.1f", lower);
            String up = (upper==Double.POSITIVE_INFINITY)? "∞" : String.format("%.1f", upper);
            return name + ": ["+lo+", "+up+")";
        }
    }

    public final Level[] levels;

    /** Creating n-1 thresholds，decreasing order；auto add -∞ / +∞ to create n categories */
    public TransitLevels(int n, Random rnd) {
//        System.out.println("DEBUG: Run Transit levels");
        if (n<2) throw new IllegalArgumentException("need ≥2 levels");
        double[] th = new double[n-1];
        for(int i=0;i<th.length;i++) th[i]= rnd.nextDouble();
        java.util.Arrays.sort(th);

        levels = new Level[n];
        double prev = Double.POSITIVE_INFINITY;
        for (int idx=0; idx<n; idx++) {
            double cur = (idx==n-1)? Double.NEGATIVE_INFINITY : th[th.length-1-idx];
            String name = "T"+(idx+1);
            levels[idx] = new Level(name, cur, prev);
            prev = cur;
        }
    }
    public TransitLevels(Level[] predefinedLevels) {
        this.levels = predefinedLevels;
    }

    /** According to the original value to find corresponding index */
    public int indexOf(double value){
        for(int i=0;i<levels.length;i++)
            if(levels[i].contains(value)) return i;
        return -1;
    }

    /* ===== 新增：更可控的等级构造 ===== */

    /** 等分分位：T1 为最高 1/n，T2 为 (1-2/n, 1-1/n] ... */
    public static TransitLevels equalQuantiles(int n) {
        Level[] lv = new Level[n];
        double prev = Double.POSITIVE_INFINITY;
        for (int idx = 0; idx < n; idx++) {
            double lo = (idx == n - 1) ? Double.NEGATIVE_INFINITY : (1.0 - (idx + 1) / (double) n);
            String name = "T" + (idx + 1);
            lv[idx] = new Level(name, lo, prev);
            prev = lo;
        }
        return new TransitLevels(lv);
    }

    /** 随机阈值 + 最小间隔（单位是 [0,1] 的绝对宽度），T1 最大。 */
    public static TransitLevels randomWithMinGap(int n, double minGap, Random rnd) {
        if (n < 2) throw new IllegalArgumentException("need ≥2 levels");
        if (minGap < 0 || minGap >= 1) throw new IllegalArgumentException("minGap in [0,1)");
        int k = n - 1;
        java.util.List<Double> th = new java.util.ArrayList<>();
        int guard = 0, MAX_TRY = 100000;

        while (th.size() < k && guard++ < MAX_TRY) {
            double v = rnd.nextDouble();
            boolean ok = true;
            for (double u : th) if (Math.abs(u - v) < minGap) { ok = false; break; }
            if (ok) th.add(v);
        }
        if (th.size() < k) {  // 兜底：退回等分分位
            return equalQuantiles(n);
        }
        th.sort(java.util.Collections.reverseOrder());
        Level[] lv = new Level[n];
        double prev = Double.POSITIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            double lo = (i == n - 1) ? Double.NEGATIVE_INFINITY : th.get(i);
            lv[i] = new Level("T" + (i + 1), lo, prev);
            prev = lo;
        }
        return new TransitLevels(lv);
    }

}




