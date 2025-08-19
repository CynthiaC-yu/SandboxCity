import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;



public class SandboxCityGUI extends JFrame {
    private final int[] SIZES = {4, 5};
    private SandboxCity city;
    private JPanel gridPanel, legendPanel, infoPanel;
    private JLabel lblHx, lblHxy, lblRatio;
    private JComboBox<Integer> cmbSize, cmbRaceCount, cmbTransitCount;
    private SandboxCity.Neighborhood selectedNb = null;
    private JPanel selectedCell   = null;
    private final Border THICK_BORDER = new LineBorder(Color.RED, 3);
    private final Border THIN_BORDER  = new LineBorder(Color.BLACK);
    private JPanel metricPanel;
    private JLabel lblD, lblG, lblC, lblR, lblP;
    // --- Variance panel widgets ---
    private JLabel lblVarH, lblVarD, lblVarG, lblVarC, lblVarR, lblVarP;
    private JPanel variancePanel;
    private JSpinner spSweepRes;
    private JPanel rightCol;
    // Variance chart
    private VarianceMiniChart varChart;
    private JCheckBox[] cbVars = new JCheckBox[6];
    private boolean[] varSelected = new boolean[]{true,true,true,true,true,true};

    private JButton btnSweep;
//    private JProgressBar pbSweep;
    private volatile boolean sweepRunning = false;

    private JSpinner betaSpinner;
    private JSpinner kappaSpinner;
    private JComboBox<String> modeCombo;
    private JButton applyParamsBtn;
    private JButton printSnapshotBtn;

    private static final double DEFAULT_BETA  = 0.3;
    private static final double DEFAULT_KAPPA = 0.6;

    private static final String MODE_HIER = "Hierarchical (π~Dir(β), p~Dir(κπ))";
    private static final String MODE_DIR1 = "Dirichlet(1)";
    private static final String MODE_UNIF = "Uniform-Normalize (legacy)";

    private JButton btnInitSeg;  // ← 新增：右侧条带要用





    public SandboxCityGUI() {
        super("Sandbox City Segregation Simulator");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLayout(new BorderLayout(8,8));

        // --- Top control panel ---
        JPanel top = new JPanel();

        // Grid size
        cmbSize = new JComboBox<>();
        for(int s: SIZES) cmbSize.addItem(s);
        cmbSize.setSelectedItem(5);
        top.add(new JLabel("Grid size:"));
        top.add(cmbSize);

        // Race count selection (min 2, max 8)
        cmbRaceCount = new JComboBox<>();
        for(int rc=2; rc<=8; rc++) cmbRaceCount.addItem(rc);
        cmbRaceCount.setSelectedItem(4);
        top.add(new JLabel("Race count:"));
        top.add(cmbRaceCount);

        // Transit levels selection (min 2, max 8)
        cmbTransitCount = new JComboBox<>();
        for(int tc=2; tc<=8; tc++) cmbTransitCount.addItem(tc);
        cmbTransitCount.setSelectedItem(5);
        top.add(new JLabel("Transit levels:"));
        top.add(cmbTransitCount);

        // Init & Compute buttons
        JButton btnInit = new JButton("Init Random Distribution");
        btnInit.addActionListener(e -> initCity());
        top.add(btnInit);

        btnInitSeg = new JButton("Init Max Segregation");
        btnInitSeg.addActionListener(e -> initMaxSegregated());
        // 不要 add 到 top，右侧条带会统一摆放


//        JButton btnInitSeg = new JButton("Init Max Segregation");
//        btnInitSeg.addActionListener(e -> initMaxSegregated());
//        top.add(btnInitSeg);

//        top.add(buildParamControls());

        // 在你已有 top 控件区初始化处（比如添加完 “Init …” 按钮之后）
//        JButton btnEnsemble = new JButton("Ensemble Analysis…");
//        btnEnsemble.addActionListener(e -> {
//            EnsembleDialog dlg = new EnsembleDialog(this);
//            dlg.open();
//        });
//        top.add(btnEnsemble);




        add(top, BorderLayout.NORTH);

        // --- Center grid panel ---
        gridPanel = new JPanel();
        add(gridPanel, BorderLayout.CENTER);

        // --- East side: legend + info vertical ---
        JPanel east = new JPanel();
        east.setLayout(new BorderLayout(0,8));   // up legend，down info
        add(east, BorderLayout.EAST);

        /* 2-a Legend */
        legendPanel = new JPanel();
        legendPanel.setLayout(new BoxLayout(legendPanel, BoxLayout.Y_AXIS));
        legendPanel.setBorder(BorderFactory.createTitledBorder("Legend"));
        east.add(legendPanel, BorderLayout.NORTH);

        /* 2-b InfoPanel */
        infoPanel = new JPanel();
        infoPanel.setLayout(new BoxLayout(infoPanel, BoxLayout.Y_AXIS));
        infoPanel.setBorder(BorderFactory.createTitledBorder("Selected Block"));
        east.add(infoPanel, BorderLayout.CENTER);

        /* ---- Metric block ---- */
        metricPanel = new JPanel();
        metricPanel.setLayout(new BoxLayout(metricPanel, BoxLayout.Y_AXIS));
        metricPanel.setBorder(BorderFactory.createTitledBorder("Metric"));

        JPanel eastWrap = new JPanel(new BorderLayout(8,0));
        eastWrap.add(east, BorderLayout.WEST);          // 左列：Legend + Info

        // right col：Metric + Variance
        rightCol = new JPanel();
        rightCol.setLayout(new BoxLayout(rightCol, BoxLayout.Y_AXIS));
        eastWrap.add(rightCol, BorderLayout.EAST);

        add(eastWrap, BorderLayout.EAST);

//        rightCol.add(metricPanel);
//        rightCol.add(Box.createVerticalStrut(8));

        buildParamControls();        // ensure spinners/buttons/combos exist

        JPanel metricHeader = new JPanel(new BorderLayout(8, 0));
        metricHeader.add(metricPanel, BorderLayout.CENTER);
        metricHeader.add(buildRightControlsStrip(), BorderLayout.EAST);


        rightCol.add(metricHeader);
        rightCol.add(Box.createVerticalStrut(8));  // 与下方 Variance 分隔


        lblHx    = new JLabel("H(X): –");
        lblHxy   = new JLabel("H(X|Y): –");
        lblRatio = new JLabel("Ratio: –");
        lblD = new JLabel("Dissimilarity D: –");
        lblG = new JLabel("Gini G: –");
        lblC = new JLabel("Variation C: –");
        lblR = new JLabel("Rel. Diversity R: –");
        lblP = new JLabel("Exposure P: –");



        initCity();

        setSize(900, 700);
        setLocationRelativeTo(null);
        setVisible(true);
    }

    private void initCity() {
        int size         = (Integer)cmbSize.getSelectedItem();
        int raceCount    = (Integer)cmbRaceCount.getSelectedItem();
        int transitCount = (Integer)cmbTransitCount.getSelectedItem();

        double beta  = (betaSpinner  != null) ? ((Number) betaSpinner.getValue()).doubleValue()  : DEFAULT_BETA;
        double kappa = (kappaSpinner != null) ? ((Number) kappaSpinner.getValue()).doubleValue() : DEFAULT_KAPPA;


        city = new SandboxCity(size, raceCount, transitCount, beta, kappa, System.nanoTime());
        city.initializeRandom();



        selectedNb = null;
        buildGridUI(size);
        updateMetrics();
        buildLegend();
        buildInfo(null);
        buildMetric();
        updateMetrics();
        buildVariance();

        syncParamSpinnersFromCity();



    }

    public SandboxCity getCity(){ return city; }


    private void syncParamSpinnersFromCity() {
        if (city == null) return;
        if (betaSpinner != null)  betaSpinner.setValue(city.getBeta());
        if (kappaSpinner != null) kappaSpinner.setValue(city.getKappa());
        if (modeCombo != null) {
            SandboxCity.GenerationMode m = city.getMode();
            modeCombo.setSelectedItem(
                    m == SandboxCity.GenerationMode.DIRICHLET1 ? MODE_DIR1 :
                            m == SandboxCity.GenerationMode.UNIFORM_NORMALIZE ? MODE_UNIF :
                                    MODE_HIER
            );
        }
    }

    private JPanel coloredSwatch(int paletteIndex, int total, float s, float b) {
        JPanel sw = new JPanel();
        sw.setPreferredSize(new Dimension(20, 20));

        float shift = 120f / 360f;
        float hue   = (paletteIndex/(float) total + shift) % 1f;

        sw.setBackground(Color.getHSBColor(hue, s, b));
        sw.setBorder(LineBorder.createBlackLineBorder());
        return sw;
    }


    private double[] computeAverageRaceProps() {
        int N = city.grid.length;
        int totalCells = N * N;
        int R = city.races.length;
        double[] sum = new double[R];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                double[] rp = city.grid[i][j].raceProps;
                for (int k = 0; k < R; k++) sum[k] += rp[k];
            }
        for (int k = 0; k < R; k++) sum[k] /= totalCells;
        return sum;
    }

    private static void leftAlign(JComponent... comps) {
        for (JComponent c : comps) {
            c.setAlignmentX(Component.LEFT_ALIGNMENT);
            Dimension pref = c.getPreferredSize();
            c.setMaximumSize(new Dimension(Integer.MAX_VALUE, pref.height));
        }
    }

    private static void commitSpinner(JSpinner sp) {
        try { sp.commitEdit(); } catch (java.text.ParseException ignore) {}
    }

    private void applyParams() {
        if (city == null) {
            JOptionPane.showMessageDialog(this, "Please initialize the city（Init）。",
                    "Hint", JOptionPane.INFORMATION_MESSAGE);
            return;
        }

        // 确保文本框里的新输入已提交
        commitSpinner(betaSpinner);
        commitSpinner(kappaSpinner);

        double beta  = ((Number) betaSpinner.getValue()).doubleValue();
        double kappa = ((Number) kappaSpinner.getValue()).doubleValue();
        String sel   = (String) modeCombo.getSelectedItem();

        SandboxCity.GenerationMode m = SandboxCity.GenerationMode.HIERARCHICAL;
        if (MODE_DIR1.equals(sel)) m = SandboxCity.GenerationMode.DIRICHLET1;
        if (MODE_UNIF.equals(sel)) m = SandboxCity.GenerationMode.UNIFORM_NORMALIZE;

        city.setParams(beta, kappa);
        city.setMode(m);
        city.initializeRandom();

        rebuildLevelsAndUI();
        // 让 spinner 反映可能的夹紧/修正
        syncParamSpinnersFromCity();
    }



    private JPanel buildParamControls() {
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT, 8, 0));

        // 只在第一次创建
        if (betaSpinner == null) {
            double betaInit  = (city != null) ? city.getBeta()  : DEFAULT_BETA;
            double kappaInit = (city != null) ? city.getKappa() : DEFAULT_KAPPA;

            betaSpinner = new JSpinner(new SpinnerNumberModel(betaInit, 1e-6, 1e6, 0.05));
            ((JSpinner.NumberEditor) betaSpinner.getEditor()).getFormat().setMinimumFractionDigits(2);

            kappaSpinner = new JSpinner(new SpinnerNumberModel(kappaInit, 1e-6, 1e6, 0.05));
            ((JSpinner.NumberEditor) kappaSpinner.getEditor()).getFormat().setMinimumFractionDigits(2);

            modeCombo = new JComboBox<>(new String[]{ MODE_HIER, MODE_DIR1, MODE_UNIF });

            applyParamsBtn = new JButton("Apply β, κ, mode");
            applyParamsBtn.setToolTipText("Apply Parameters and randomize the Sandbox City");
            applyParamsBtn.addActionListener(e -> applyParams());  // 抽成单独方法，见下

            printSnapshotBtn = new JButton("Print snapshot");
            printSnapshotBtn.setToolTipText("Print current π/avg/unit maximum ratio range");
            printSnapshotBtn.addActionListener(e -> {
                if (city == null) {
                    JOptionPane.showMessageDialog(this, "City is not initialized.", "Hint", JOptionPane.INFORMATION_MESSAGE);
                    return;
                }
                String s = city.snapshotString();
                System.out.println(s);
                JOptionPane.showMessageDialog(
                        this,
                        "<html><pre style='font-family:monospace'>" + s + "</pre></html>",
                        "Snapshot",
                        JOptionPane.INFORMATION_MESSAGE
                );
            });
        }



        // 用“同一套控件”搭 panel（可多次调用，但不会 new 新控件）
        panel.add(new JLabel("β")); panel.add(betaSpinner);
        panel.add(Box.createHorizontalStrut(8));
        panel.add(new JLabel("κ")); panel.add(kappaSpinner);
        panel.add(Box.createHorizontalStrut(8));
        panel.add(new JLabel("Mode")); panel.add(modeCombo);
        panel.add(Box.createHorizontalStrut(8));
        panel.add(applyParamsBtn);
        panel.add(printSnapshotBtn);
        return panel;
    }



    /** Build the Variance panel (sweep over thresholds) and attach under metricPanel */
    private static final String VAR_PANEL_NAME = "VARIANCE_PANEL";

    private JPanel makeVarRow(int idx, JLabel label){
        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT,5,0));
        JCheckBox cb = new JCheckBox();
        cb.setSelected(varSelected[idx]);       // read the choice
        cb.addItemListener(e -> {
            varSelected[idx] = cb.isSelected(); // update memory
            if (varChart != null) varChart.setMask(varSelected);
            label.setEnabled(cb.isSelected());  // making the label to be gray when its not selected
        });
        cbVars[idx] = cb;

        row.add(cb);
        row.add(label);
        leftAlign(row);
        return row;
    }

    private void setVarianceControlsEnabled(boolean en){
        if (spSweepRes != null) spSweepRes.setEnabled(en);
        if (btnSweep   != null) btnSweep.setEnabled(en);
        if (cbVars     != null) for (JCheckBox cb : cbVars) if (cb != null) cb.setEnabled(en);
    }

    private void runSweepAsync(){
        if (city == null || sweepRunning) return;

        double res = ((Number) spSweepRes.getValue()).doubleValue();
        if (res <= 0.0 || res > 1.0) {
            JOptionPane.showMessageDialog(this, "Resolution must be in (0, 1].",
                    "Invalid resolution", JOptionPane.WARNING_MESSAGE);
            return;
        }

        sweepRunning = true;
        setVarianceControlsEnabled(false);

        SwingWorker<ThresholdSweepAnalyzer.SweepResult, Void> worker =
                new SwingWorker<ThresholdSweepAnalyzer.SweepResult, Void>() {
                    @Override protected ThresholdSweepAnalyzer.SweepResult doInBackground() {
                        //Another thread to run sweep, so it won't block the UI thread
                        return ThresholdSweepAnalyzer.sweep(city, res);
                    }
                    @Override protected void done() {
                        try {
                            ThresholdSweepAnalyzer.SweepResult r = get();
                            updateVariancePanel(r);                  // 刷 label + 图表
                        } catch (Exception ex) {
                            JOptionPane.showMessageDialog(
                                    SandboxCityGUI.this,
                                    "Sweep failed: " + ex.getMessage(),
                                    "Error",
                                    JOptionPane.ERROR_MESSAGE
                            );
                        } finally {
                            sweepRunning = false;
                            setVarianceControlsEnabled(true);
//                            if (pbSweep != null) pbSweep.setVisible(false);
                        }
                    }
                };
        worker.execute();
    }


    private void buildVariance(){
        if (variancePanel != null) rightCol.remove(variancePanel);

        variancePanel = new JPanel();
        variancePanel.setLayout(new BoxLayout(variancePanel, BoxLayout.Y_AXIS));
        variancePanel.setBorder(BorderFactory.createTitledBorder("Variance"));

        JLabel title = new JLabel("<html><b>sweep over thresholds</b></html>");
        lblVarH = new JLabel("Var(H): –");
        lblVarD = new JLabel("Var(D): –");
        lblVarG = new JLabel("Var(G): –");
        lblVarC = new JLabel("Var(C): –");
        lblVarR = new JLabel("Var(R): –");
        lblVarP = new JLabel("Var(P): –");

        variancePanel.add(title);
        variancePanel.add(makeVarRow(0, lblVarH));
        variancePanel.add(makeVarRow(1, lblVarD));
        variancePanel.add(makeVarRow(2, lblVarG));
        variancePanel.add(makeVarRow(3, lblVarC));
        variancePanel.add(makeVarRow(4, lblVarR));
        variancePanel.add(makeVarRow(5, lblVarP));

// Resolution
        JPanel resRow = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        resRow.add(new JLabel("Resolution:"));
        spSweepRes = new JSpinner(new SpinnerNumberModel(0.05, 0.01, 1.00, 0.01));
        ((JSpinner.DefaultEditor) spSweepRes.getEditor()).getTextField().setColumns(4);
        resRow.add(spSweepRes);

        btnSweep = new JButton("Sweep");
        btnSweep.addActionListener(e -> runSweepAsync());
        resRow.add(btnSweep);

        variancePanel.add(Box.createVerticalStrut(4));
        variancePanel.add(resRow);

        // Chart
        variancePanel.add(Box.createVerticalStrut(6));
        varChart = new VarianceMiniChart();
        variancePanel.add(varChart);
        // Initial Mask
        varChart.setMask(varSelected);
        varChart.setValues(null);

        leftAlign(title, variancePanel, resRow, varChart);

        rightCol.add(variancePanel);
        rightCol.revalidate();
        rightCol.repaint();

    }


    private void buildMetric(){
        metricPanel.removeAll();

        metricPanel.add(new JLabel("<html><b>Segregation Score:</b></html>"));
        metricPanel.add(Box.createVerticalStrut(8));

        metricPanel.add(lblRatio);
        metricPanel.add(lblD);
        metricPanel.add(lblG);
        metricPanel.add(lblC);
        metricPanel.add(lblR);
        metricPanel.add(lblP);

        leftAlign( lblRatio, lblD, lblG, lblC, lblR, lblP);

        metricPanel.revalidate();
        metricPanel.repaint();



    }


    private void buildInfo(SandboxCity.Neighborhood nb) {
        infoPanel.removeAll();
        if (nb==null) { infoPanel.add(new JLabel("Nothing selected")); return; }

        int raceLength = city.races.length;


        /* Transit value spinner */
        JSpinner spVal = new JSpinner(
                new SpinnerNumberModel(nb.transitValue, 0, 1.0, 0.01));


        spVal.addChangeListener(e -> {
            nb.transitValue = ((Number)spVal.getValue()).doubleValue();
            nb.level = city.transit.indexOf(nb.transitValue);
            rebuildLevelsAndUI();
        });

        JPanel raceBox = new JPanel(new GridLayout(raceLength,3,4,2));
        double[] rp = nb.raceProps;
        JSpinner[] spRace = new JSpinner[rp.length];

        for(int k=0;k<rp.length;k++){

            spRace[k] = new JSpinner(
                    new SpinnerNumberModel(rp[k], 0.0, 1.0, 0.01));

            raceBox.add(coloredSwatch(k, raceLength, 1f, 0.7f));
            raceBox.add(new JLabel(city.races[k]));
            raceBox.add(spRace[k]);
        }

        JButton btnApply = new JButton("Apply");
        btnApply.addActionListener(e -> {
            double sum = 0;
            for(int k=0;k<rp.length;k++){
                rp[k] = ((Number)spRace[k].getValue()).doubleValue();
                sum  += rp[k];
            }
            // Normalized to ensure percentages sum to 1
            for(int k=0;k<rp.length;k++) rp[k] /= sum;
            rebuildLevelsAndUI();
        });


        JPanel itemTransit = new JPanel(new GridLayout(0,1));
            String levelName = city.transitCategories[nb.level];
            JLabel transitInfo = new JLabel(
                    String.format("Level: " + levelName));
            transitInfo.setFont(transitInfo.getFont().deriveFont(Font.BOLD));
            itemTransit.add(transitInfo);

            itemTransit.add(new JLabel("<html><b>Transit Value: </b></html>"));

        JPanel titlePanelRace = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
            titlePanelRace.add(new JLabel("<html><b>Race Proportion:</b></html>"));



        /* Setting up infoPanel */
        infoPanel.add(itemTransit);
        infoPanel.add(spVal);
        infoPanel.add(Box.createVerticalStrut(6));
        infoPanel.add(titlePanelRace);
        infoPanel.add(raceBox);
        infoPanel.add(btnApply);
        infoPanel.revalidate(); infoPanel.repaint();

    }


    private void buildLegend() {

        legendPanel.removeAll();

        JPanel titlePanelTransit = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        titlePanelTransit.add(new JLabel("<html><b>Transit Level:</b></html>\""));
        legendPanel.add(titlePanelTransit);

        int transitLength = city.transit.levels.length;

        for (int i = 0; i < transitLength - 1; i++) {   // T1 has no upper bd, Tn no lower bd
            TransitLevels.Level L = city.transit.levels[i];
            double init = L.lower;
            JSpinner sp = new JSpinner(new SpinnerNumberModel(init, 0, 1.001, 0.01));
            int idx = i;                                         // for lambda 捕获
            sp.addChangeListener(e -> {
                double bound = ((Number)sp.getValue()).doubleValue();
                city.transit.levels[idx].lower = bound;
                city.transit.levels[idx+1].upper = bound;
                rebuildLevelsAndUI();
            });

            JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT,5,2));
            row.add(coloredSwatch(i, transitLength, 0.7f, 1f));
            row.add(new JLabel("T"+(i+1)+" ≥ "));
            row.add(sp);
            legendPanel.add(row);
        }

        TransitLevels.Level last = city.transit.levels[transitLength-1];
        JPanel rowN = new JPanel(new FlowLayout(FlowLayout.LEFT,5,2));
        rowN.add(coloredSwatch(transitLength-1, transitLength, 0.7f, 1));
        String txt_Tn= String.format("T%d < %.2f", transitLength, last.upper);
        rowN.add(new JLabel(txt_Tn));
        legendPanel.add(rowN);

        legendPanel.add(Box.createVerticalStrut(10));

        JPanel titlePanelRace = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 0));
        titlePanelRace.add(new JLabel("<html><b>Average Race Proportion:</b></html>\""));
        legendPanel.add(titlePanelRace);


        double[] display = computeAverageRaceProps();

        for (int k = 0; k < city.races.length; k++) {
            JPanel item = new JPanel(new FlowLayout(FlowLayout.LEFT,5,2));
            JPanel sw = new JPanel(); sw.setPreferredSize(new Dimension(20,20));

            int n = city.races.length;
            float shift = 120f / 360f;
            float hue   = (k/(float)n + shift) % 1f;
            sw.setBackground(Color.getHSBColor(hue,1f,0.7f));

            sw.setBorder(LineBorder.createBlackLineBorder());
            String txt = String.format("%s: %.1f%%", city.races[k], display[k]*100);
            item.add(sw);
            item.add(new JLabel(txt));
            legendPanel.add(item);
        }

        legendPanel.revalidate();
        legendPanel.repaint();
    }

    private void buildGridUI(int size) {
        gridPanel.removeAll();
        gridPanel.setLayout(new GridLayout(size, size, 4, 4));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                SandboxCity.Neighborhood nb = city.grid[i][j];
                JPanel cell = new JPanel(new BorderLayout());
                cell.setBorder(THIN_BORDER);

                if (nb == selectedNb) {
                    cell.setBorder(THICK_BORDER);
                    selectedCell = cell;
                }
//                cell.setBorder(new LineBorder(Color.BLACK));

                int n =  city.transitCategories.length;
                float shift = 120f / 360f;
                float hue = (nb.level / (float) n + shift) % 1;

                cell.setBackground(Color.getHSBColor(hue, 0.7f, 1.0f));
                cell.setOpaque(true);
                cell.add(new RaceBar(nb.raceProps), BorderLayout.SOUTH);



                cell.addMouseListener(new MouseAdapter() {
                    @Override public void mouseClicked(MouseEvent e) {

                        if (selectedCell != null) {
                            selectedCell.setBorder(THIN_BORDER);
                        }
                        selectedCell = cell;
                        selectedCell.setBorder(THICK_BORDER);

                        selectedNb = nb;
                        buildInfo(selectedNb);

                    }
                });


                gridPanel.add(cell);
            }
        }
        gridPanel.revalidate();
        gridPanel.repaint();
    }

    private void initMaxSegregated() {
        int size         = (Integer) cmbSize.getSelectedItem();
        int raceCount    = (Integer) cmbRaceCount.getSelectedItem();
        int transitCount = (Integer) cmbTransitCount.getSelectedItem();

        // 新建城市并生成网格对象（用随机的阈值，但我们会覆盖每个格子的内容）
        city = new SandboxCity(size, raceCount, transitCount);
        city.initializeRandom(); // 仅用于把 grid[][] 填好 Neighborhood 对象

        // 若 race 数与 transit 档数不同，提示一下；映射仍然会循环分配
        if (raceCount != transitCount) {
            JOptionPane.showMessageDialog(
                    this,
                    "Race count ("+raceCount+") and Transit levels ("+transitCount+") is not equal; recycle bins.",
                    "Info",
                    JOptionPane.INFORMATION_MESSAGE
            );
        }

        // ① For every transit level randomly assign a race
        int T = city.transit.levels.length;
        int R = city.races.length;
        List<Integer> order = new ArrayList<>();
        for (int r = 0; r < R; r++) order.add(r);
        Collections.shuffle(order, new Random());
        int[] levelToRace = new int[T];
        for (int t = 0; t < T; t++) levelToRace[t] = order.get(t % R);

        // ② Making every level has about the number of bins
        int total = size * size;
        int base = total / T, extra = total % T;
        int[] quota = new int[T];
        for (int t = 0; t < T; t++) quota[t] = base + (t < extra ? 1 : 0);

        Random rnd = new Random();
        int curLevel = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                while (curLevel < T && quota[curLevel] == 0) curLevel++;
                if (curLevel >= T) curLevel = T - 1; // 兜底

                SandboxCity.Neighborhood nb = city.grid[i][j];

                // 1-hot 的 race 向量
                Arrays.fill(nb.raceProps, 0.0);
                int raceIdx = levelToRace[curLevel];
                nb.raceProps[raceIdx] = 1.0;

                // transit 值：落在该 level 的 [lower, upper) 区间内
                TransitLevels.Level L = city.transit.levels[curLevel];
                double lo = Double.isInfinite(L.lower) ? 0.0 : L.lower;
                double up = Double.isInfinite(L.upper) ? 1.0 : L.upper;
                if (up <= lo) { // 阈值异常时的兜底
                    nb.transitValue = Math.max(0.0, Math.min(1.0, lo));
                } else {
                    nb.transitValue = lo + (up - lo) * rnd.nextDouble();
                }
                nb.level = curLevel;

                quota[curLevel]--;
            }
        }

        selectedNb = null;
        buildGridUI(size);
        buildLegend();
        buildInfo(null);
        buildMetric();
        updateMetrics();
    }

    private JPanel buildRightControlsStrip() {
//        if (betaSpinner == null || kappaSpinner == null || modeCombo == null
//                || applyParamsBtn == null || printSnapshotBtn == null) {
//            buildParamControls();
//        }

        JPanel wrap = new JPanel();
        wrap.setLayout(new BoxLayout(wrap, BoxLayout.Y_AXIS));
        wrap.setBorder(BorderFactory.createTitledBorder("Controls"));

        JPanel r0 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        r0.add(btnInitSeg);

        JPanel r1 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        r1.add(new JLabel("β"));
        r1.add(betaSpinner);
        r1.add(Box.createHorizontalStrut(6));
        r1.add(new JLabel("κ"));
        r1.add(kappaSpinner);

        JPanel r2 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        r2.add(new JLabel("Mode"));
        r2.add(modeCombo);

        JPanel r3 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        r3.add(applyParamsBtn);

        JPanel r4 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        r4.add(printSnapshotBtn);

        JPanel r5 = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 0));
        JButton btnEnsemble = new JButton("Ensemble Analysis…");
        btnEnsemble.addActionListener(e -> new EnsembleDialog(this).open());
        r5.add(btnEnsemble);

        leftAlign(r0, r1, r2, r3, r4, r5);
        wrap.add(r0); wrap.add(Box.createVerticalStrut(4));
        wrap.add(r1); wrap.add(r2); wrap.add(r3); wrap.add(r4); wrap.add(r5);
        return wrap;
    }



//    private JPanel buildRightControlsStrip() {
//        JPanel strip = new JPanel(new FlowLayout(FlowLayout.RIGHT, 8, 0));
//
//////        JButton btnInitSeg = this.btnInitSeg;         // 你原来的“Init Max Segregation”按钮
////        strip.add(btnInitSeg);
////
////        strip.add(new JLabel("β"));
////        strip.add(this.betaSpinner);
////        strip.add(new JLabel("κ"));
////        strip.add(this.kappaSpinner);
////
////        strip.add(new JLabel("Mode"));
////        strip.add(this.modeCombo);
////
////        strip.add(this.applyParamsBtn);
////        strip.add(this.printSnapshotBtn);
//
//        JButton btnEnsemble = new JButton("Ensemble Analysis…");
//        btnEnsemble.addActionListener(e -> new EnsembleDialog(this).open());
//        strip.add(btnEnsemble);
//
//        return strip;
//    }





    // Updating UI：
    private void rebuildLevelsAndUI() {
        /* Updating the blocks*/
        for (int i = 0; i < city.getSize(); i++) {
            for (int j = 0; j < city.getSize(); j++) {
                SandboxCity.Neighborhood nb = city.grid[i][j];
                nb.level = city.transit.indexOf(nb.transitValue);   // <- 更新 level
            }
        }

        buildGridUI(city.getSize());
        updateMetrics();
        buildLegend();

        if (selectedNb != null) buildInfo(selectedNb);

        gridPanel.repaint();
    }



    private void updateMetrics() {
//        lblHx   .setText(String.format("H(X): %.4f", city.computeTotalEntropy()));
//        lblHxy  .setText(String.format("H(X|Y): %.4f", city.computeConditionalEntropy()));
        lblRatio.setText(String.format("Theil H: %.4f", city.computeEntropyRatio()));

        lblD.setText(String.format("Dissimilarity D: %.4f",  city.computeDissimilarity()));
        lblG.setText(String.format("Gini G: %.4f",           city.computeGini()));
        lblC.setText(String.format("Variation C: %.4f",      city.computeVariation()));
        lblR.setText(String.format("Rel. Diversity R: %.4f", city.computeRelativeDiversity()));
        lblP.setText(String.format("Exposure P: %.4f",       city.computeExposure()));
    }

    private void updateVariancePanel(ThresholdSweepAnalyzer.SweepResult r){
        lblVarH.setText(String.format("Var(H): %.6e", r.varH));
        lblVarD.setText(String.format("Var(D): %.6e", r.varD));
        lblVarG.setText(String.format("Var(G): %.6e", r.varG));
        lblVarC.setText(String.format("Var(C): %.6e", r.varC));
        lblVarR.setText(String.format("Var(R): %.6e", r.varR));
        lblVarP.setText(String.format("Var(P): %.6e", r.varP));
        if (varChart != null) {
            varChart.setValues(new double[]{ r.varH, r.varD, r.varG, r.varC, r.varR, r.varP });
            varChart.setMask(varSelected);
        }
        variancePanel.revalidate();
        variancePanel.repaint();
    }




    class RaceBar extends JComponent {
        double[] props;
        RaceBar(double[] props) {
            this.props = props;
            setPreferredSize(new Dimension(0,20));
        }
        protected void paintComponent(Graphics g) {
            int w = getWidth(), h = getHeight(), x0 = 0;
            for (int k = 0; k < props.length; k++) {

                int n = props.length;
                float shift = 120f / 360f;
                float hue   = (k/(float)n + shift) % 1f;
                int wseg = (int)Math.round(props[k]*w);
                g.setColor(Color.getHSBColor(hue, 1f, 0.7f));
                g.fillRect(x0, 0, wseg, h);
                x0 += wseg;
            }
        }
    }

    private static class VarianceMiniChart extends JComponent {
        private final String[] labels = {"H","D","G","C","R","P"};
        private final double[] vals   = new double[6];
        private final boolean[] mask  = new boolean[]{true,true,true,true,true,true};

        VarianceMiniChart() { setPreferredSize(new Dimension(180, 90)); }

        void setValues(double[] v){
            if (v != null) {
                int n = Math.min(vals.length, v.length);
                for (int i=0;i<n;i++) vals[i] = v[i];
            } else {
                for (int i=0;i<vals.length;i++) vals[i] = Double.NaN;
            }
            repaint();
        }

        void setMask(boolean[] m){
            if (m != null) {
                for (int i=0;i<mask.length;i++) mask[i] = (i < m.length) ? m[i] : true;
            }
            repaint();
        }

        @Override protected void paintComponent(Graphics g) {
            Graphics2D g2 = (Graphics2D) g.create();
            try {
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                int w = getWidth(), h = getHeight();
                int pad = 6, bottom = h - pad - 12;
                int top = pad, availH = Math.max(1, bottom - top);

                // counts the selected metrics
                int sel = 0; for (boolean b : mask) if (b) sel++;
                if (sel == 0) {
                    g2.setColor(new Color(0,0,0,120));
                    g2.setFont(getFont().deriveFont(12f));
                    String msg = "No bars selected";
                    int sw = g2.getFontMetrics().stringWidth(msg);
                    g2.drawString(msg, (w-sw)/2, (h+g2.getFontMetrics().getAscent())/2 - 4);
                    return;
                }

                int gap = 6;
                int bw = Math.max(4, (w - pad*2 - gap*(sel-1)) / sel);

                double max = 0;
                for (int i=0;i<vals.length;i++)
                    if (mask[i] && !Double.isNaN(vals[i])) max = Math.max(max, vals[i]);
                if (max <= 0) max = 1;

                g2.setColor(new Color(0,0,0,60));
                g2.drawLine(pad, bottom, w - pad, bottom);

                int j = 0; // 选中序号
                for (int i=0;i<labels.length;i++){
                    if (!mask[i]) continue;
                    double v = Double.isNaN(vals[i]) ? 0 : vals[i];
                    int bh = (int)Math.round(v / max * (availH - 2));
                    int x = pad + j * (bw + gap);
                    int y = bottom - bh;

                    float hue = (i/(float)labels.length + 0.60f) % 1f;
                    g2.setColor(Color.getHSBColor(hue, 0.50f, 0.95f));
                    g2.fillRoundRect(x, y, bw, bh, 6, 6);
                    g2.setColor(new Color(0,0,0,80));
                    g2.drawRoundRect(x, y, bw, bh, 6, 6);

                    g2.setFont(getFont().deriveFont(10f));
                    FontMetrics fm = g2.getFontMetrics();
                    int lx = x + (bw - fm.stringWidth(labels[i]))/2;
                    g2.setColor(new Color(0,0,0,120));
                    g2.drawString(labels[i], lx, bottom + fm.getAscent());
                    j++;
                }
            } finally { g2.dispose(); }
        }
    }


}
