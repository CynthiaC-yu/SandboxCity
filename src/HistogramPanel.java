import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.Arrays;

public class HistogramPanel extends JPanel {
    private double[] data;
    private int bins = 20;
    private boolean auto = true;
    private String title = "";
    private double min, max;
    private int[] counts;

    public HistogramPanel() {
        setPreferredSize(new Dimension(340, 220));
    }

    public void setData(String title, double[] data, boolean auto, int fixedBins) {
        this.title = title;
        this.data  = (data == null)? new double[0] : data.clone();
        this.auto  = auto;
        this.bins  = Math.max(1, fixedBins);
        recompute();
        repaint();
    }

    private static double zfix(double v){ return Math.abs(v) < 1e-12 ? 0.0 : v; }

    private static final DecimalFormat SCI = new DecimalFormat("0.###E0");
    private static String fmt(double v){
        v = zfix(v);
        double a = Math.abs(v);
        if ((a > 0 && a < 1e-3) || a >= 1e3) return SCI.format(v);   // 科学计数
        return String.format("%.4f", v);                              // 否则保留 4 位小数
    }

    private void recompute() {
        if (data.length == 0) { counts = new int[1]; min = 0; max = 1; return; }
        double[] clean = Arrays.stream(data).filter(d -> !Double.isNaN(d) && !Double.isInfinite(d)).toArray();
        if (clean.length == 0) { counts = new int[1]; min = 0; max = 1; return; }

        min = Arrays.stream(clean).min().orElse(0);
        max = Arrays.stream(clean).max().orElse(1);

        if (auto) {
            Arrays.sort(clean);
            int n = clean.length;
            int q1i = n/4, q3i = (3*n)/4;
            double iqr = Math.max(1e-12, clean[q3i] - clean[q1i]);
            double width = 2 * iqr / Math.cbrt(n);                // Freedman–Diaconis
            bins = (width <= 0) ? Math.max(5, (int)Math.round(Math.sqrt(n)))
                    : Math.max(1, (int)Math.ceil((max - min) / width));
            bins = Math.min(80, bins);
        }
        if (zfix(max - min) == 0.0) { counts = new int[]{ clean.length }; return; } // 常数分布

        counts = new int[bins];
        for (double v : clean) {
            int idx = (int) Math.floor((v - min) / (max - min) * bins);
            if (idx >= bins) idx = bins - 1;
            if (idx < 0) idx = 0;
            counts[idx]++;
        }
    }

    @Override protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        var g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        int W = getWidth(), H = getHeight();
        int top = 28, left = 60, right = 12, bottom = 40;

        // Title
        g2.setFont(getFont().deriveFont(Font.BOLD, 12f));
        g2.drawString(title, 10, 18);

        // Axis box
        int x0 = left, y0 = H - bottom, w = W - left - right, h = H - top - bottom;
        g2.drawRect(x0, y0 - h, w, h);

        if (counts == null) return;

        if (counts.length == 1) {
            // 常数分布：画一根居中的条，避免“整屏铺满”
            int bh = (int) (h * 0.85);
            int bw = Math.max(6, Math.min(w/8, 40));
            int x = x0 + w/2 - bw/2;
            int y = y0 - bh;
            g2.setColor(new Color(100, 140, 220));
            g2.fillRect(x, y, bw, bh);
        } else {
            int maxC = Arrays.stream(counts).max().orElse(1);
            double barW = (double) w / counts.length;
            for (int i = 0; i < counts.length; i++) {
                int bh = (int) Math.round(h * (counts[i] / (double) maxC));
                int x = (int) Math.round(x0 + i * barW) + 1;
                int bw = (int) Math.round(barW) - 2;
                int y = y0 - bh;
                g2.setColor(new Color(100, 140, 220));
                g2.fillRect(x, y, Math.max(1, bw), bh);
            }
        }

        // X ticks: min / mid / max —— 科学计数法或 4 位小数
        g2.setColor(Color.DARK_GRAY);
        g2.setFont(getFont().deriveFont(11f));
        String sMin = fmt(min);
        String sMid = fmt((min+max)/2.0);
        String sMax = fmt(max);
        g2.drawString(sMin, x0 - 6, y0 + 18);
        g2.drawString(sMid, x0 + w/2 - 18, y0 + 18);
        g2.drawString(sMax, x0 + w - 56, y0 + 18);

        // mean ± sd
        double[] clean = Arrays.stream(data).filter(d -> !Double.isNaN(d) && !Double.isInfinite(d)).toArray();
        double mean = Arrays.stream(clean).average().orElse(Double.NaN);
        double sd = Double.isNaN(mean) ? Double.NaN :
                Math.sqrt(Arrays.stream(clean).map(v->(v-mean)*(v-mean)).average().orElse(0));
        String stats = String.format("mean=%s, sd=%s, n=%d", fmt(mean), fmt(sd), clean.length);
        g2.drawString(stats, 10, H - 10);
    }
}
