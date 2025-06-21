/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package FiniteDifference1;

/**
 *
 * @author Roddier
 */
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.util.List;

public class ODEFiniteDifference {
    
    // Interface pour les fonctions f(x)
    interface Function {
        double apply(double x);
        String getName();
    }
    
    // Implémentations des fonctions
    static class SinFunction implements Function {
        public double apply(double x) { return Math.sin(Math.PI * x); }
        public String getName() { return "sin(πx)"; }
    }
    
    static class CubicFunction implements Function {
        public double apply(double x) { return x * x * x; }
        public String getName() { return "x³"; }
    }
    
    // Types d'équations différentielles
    enum EquationType {
        TYPE1("-u'' + u = f"),
        TYPE2("-u'' + u' = f"),
        TYPE3("-u'' = f");
        
        private final String description;
        EquationType(String description) { this.description = description; }
        public String getDescription() { return description; }
    }
    
    // Classe pour stocker les résultats
    static class Solution {
        double[] x;
        double[] numerical;
        double[] analytical;
        double error;
        int n;
        EquationType type;
        Function function;
        
        Solution(double[] x, double[] numerical, double[] analytical, double error, int n, EquationType type, Function function) {
            this.x = x;
            this.numerical = numerical;
            this.analytical = analytical;
            this.error = error;
            this.n = n;
            this.type = type;
            this.function = function;
        }
    }
    
    // Résolution de système tridiagonal
    private static double[] solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
        int n = b.length;
        double[] cp = new double[n];
        double[] dp = new double[n];
        double[] x = new double[n];
        
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        
        for (int i = 1; i < n; i++) {
            double denom = b[i] - a[i] * cp[i-1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i-1]) / denom;
        }
        
        x[n-1] = dp[n-1];
        for (int i = n-2; i >= 0; i--) {
            x[i] = dp[i] - cp[i] * x[i+1];
        }
        
        return x;
    }
    
    // Solutions analytiques
    private static double getAnalyticalSolution(double x, EquationType type, Function function) {
        if (function instanceof SinFunction) {
            switch (type) {
                case TYPE1: // -u'' + u = sin(πx)
                    return Math.sin(Math.PI * x) / (1 + Math.PI * Math.PI);
                case TYPE2: // -u'' + u' = sin(πx)
                    return (Math.sin(Math.PI * x) - Math.PI * Math.cos(Math.PI * x)) / (1 + Math.PI * Math.PI);
                case TYPE3: // -u'' = sin(πx)
                    return -Math.sin(Math.PI * x) / (Math.PI * Math.PI);
            }
        } else if (function instanceof CubicFunction) {
            switch (type) {
                case TYPE1: // -u'' + u = x³
                    return x * x * x - 6 * x;
                case TYPE2: // -u'' + u' = x³
                    return x * x * x - 6 * x + 6;
                case TYPE3: // -u'' = x³
                    return -x * x * x * x * x / 20 + x * x * x / 6;
            }
        }
        return 0;
    }
    
    // Résolution numérique
    public static Solution solve(int n, EquationType type, Function function) {
        double h = 1.0 / n;
        double[] x = new double[n+1];
        for (int i = 0; i <= n; i++) {
            x[i] = i * h;
        }
        
        // Matrices pour le système Ax = b
        double[] a = new double[n-1]; // sous-diagonale
        double[] b = new double[n-1]; // diagonale
        double[] c = new double[n-1]; // sur-diagonale
        double[] d = new double[n-1]; // second membre
        
        // Construction du système selon le type d'équation
        switch (type) {
            case TYPE1: // -u'' + u = f
                for (int i = 0; i < n-1; i++) {
                    a[i] = (i > 0) ? -1.0/(h*h) : 0;
                    b[i] = 2.0/(h*h) + 1.0;
                    c[i] = (i < n-2) ? -1.0/(h*h) : 0;
                    d[i] = function.apply(x[i+1]);
                }
                break;
                
            case TYPE2: // -u'' + u' = f
                for (int i = 0; i < n-1; i++) {
                    a[i] = (i > 0) ? -1.0/(h*h) - 1.0/(2*h) : 0;
                    b[i] = 2.0/(h*h);
                    c[i] = (i < n-2) ? -1.0/(h*h) + 1.0/(2*h) : 0;
                    d[i] = function.apply(x[i+1]);
                }
                break;
                
            case TYPE3: // -u'' = f
                for (int i = 0; i < n-1; i++) {
                    a[i] = (i > 0) ? -1.0/(h*h) : 0;
                    b[i] = 2.0/(h*h);
                    c[i] = (i < n-2) ? -1.0/(h*h) : 0;
                    d[i] = function.apply(x[i+1]);
                }
                break;
        }
        
        // Résolution
        double[] u_inner = solveTridiagonal(a, b, c, d);
        
        // Construction de la solution complète (avec conditions aux limites u(0) = u(1) = 0)
        double[] numerical = new double[n+1];
        double[] analytical = new double[n+1];
        
        numerical[0] = 0;
        numerical[n] = 0;
        for (int i = 1; i < n; i++) {
            numerical[i] = u_inner[i-1];
        }
        
        // Solution analytique
        for (int i = 0; i <= n; i++) {
            analytical[i] = getAnalyticalSolution(x[i], type, function);
        }
        
        // Calcul de l'erreur L2
        double error = 0;
        for (int i = 0; i <= n; i++) {
            error += Math.pow(numerical[i] - analytical[i], 2);
        }
        error = Math.sqrt(error * h);
        
        return new Solution(x, numerical, analytical, error, n, type, function);
    }
    
    // Calcul de l'ordre de convergence
    public static double calculateConvergenceOrder(double error1, double error2, int n1, int n2) {
        return Math.log(error1 / error2) / Math.log((double)n2 / n1);
    }
    
    // Classe pour la visualisation graphique
    static class GraphPanel extends JPanel {
        private List<Solution> solutions;
        private boolean showError;
        
        public GraphPanel(List<Solution> solutions, boolean showError) {
            this.solutions = solutions;
            this.showError = showError;
            setPreferredSize(new Dimension(800, 600));
            setBackground(Color.WHITE);
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int margin = 50;
            
            // Axes
            g2.setColor(Color.BLACK);
            g2.drawLine(margin, height + margin, width + margin, height + margin); // x-axis
            g2.drawLine(margin, margin, margin, height + margin); // y-axis
            
            if (showError) {
                drawErrorGraph(g2, margin, width, height);
            } else {
                drawSolutionGraph(g2, margin, width, height);
            }
        }
        
        private void drawSolutionGraph(Graphics2D g2, int margin, int width, int height) {
            if (solutions.isEmpty()) return;
            
            Solution sol = solutions.get(0);
            
            // Trouver les limites
            double minY = Double.MAX_VALUE, maxY = Double.MIN_VALUE;
            for (int i = 0; i < sol.x.length; i++) {
                minY = Math.min(minY, Math.min(sol.numerical[i], sol.analytical[i]));
                maxY = Math.max(maxY, Math.max(sol.numerical[i], sol.analytical[i]));
            }
            
            double range = maxY - minY;
            if (range == 0) range = 1;
            
            // Dessiner la solution analytique
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * width);
                int y1 = margin + height - (int)((sol.analytical[i] - minY) / range * height);
                int x2 = margin + (int)(sol.x[i+1] * width);
                int y2 = margin + height - (int)((sol.analytical[i+1] - minY) / range * height);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            // Dessiner la solution numérique
            g2.setColor(Color.RED);
            g2.setStroke(new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{5}, 0));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * width);
                int y1 = margin + height - (int)((sol.numerical[i] - minY) / range * height);
                int x2 = margin + (int)(sol.x[i+1] * width);
                int y2 = margin + height - (int)((sol.numerical[i+1] - minY) / range * height);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            // Légende
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            g2.drawString(sol.type.getDescription() + " avec f = " + sol.function.getName(), margin, 25);
            g2.drawString("n = " + sol.n, margin, 45);
            
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.setColor(Color.BLUE);
            g2.drawString("— Solution analytique", width - 150, margin + 20);
            g2.setColor(Color.RED);
            g2.drawString("--- Solution numérique", width - 150, margin + 40);
            g2.setColor(Color.BLACK);
            g2.drawString(String.format("Erreur L2: %.2e", sol.error), width - 150, margin + 60);
        }
        
        private void drawErrorGraph(Graphics2D g2, int margin, int width, int height) {
            if (solutions.size() < 2) return;
            
            double minError = solutions.stream().mapToDouble(s -> s.error).min().orElse(1e-10);
            double maxError = solutions.stream().mapToDouble(s -> s.error).max().orElse(1);
            
            // Échelle logarithmique
            double logMinError = Math.log10(minError);
            double logMaxError = Math.log10(maxError);
            double logRange = logMaxError - logMinError;
            if (logRange == 0) logRange = 1;
            
            int minN = solutions.stream().mapToInt(s -> s.n).min().orElse(10);
            int maxN = solutions.stream().mapToInt(s -> s.n).max().orElse(1000);
            double logMinN = Math.log10(minN);
            double logMaxN = Math.log10(maxN);
            double logNRange = logMaxN - logMinN;
            
            // Dessiner les points et les relier
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            
            for (int i = 0; i < solutions.size() - 1; i++) {
                Solution s1 = solutions.get(i);
                Solution s2 = solutions.get(i + 1);
                
                int x1 = margin + (int)((Math.log10(s1.n) - logMinN) / logNRange * width);
                int y1 = margin + height - (int)((Math.log10(s1.error) - logMinError) / logRange * height);
                int x2 = margin + (int)((Math.log10(s2.n) - logMinN) / logNRange * width);
                int y2 = margin + height - (int)((Math.log10(s2.error) - logMinError) / logRange * height);
                
                g2.drawLine(x1, y1, x2, y2);
                g2.fillOval(x1 - 3, y1 - 3, 6, 6);
            }
            
            // Dernier point
            Solution lastSol = solutions.get(solutions.size() - 1);
            int lastX = margin + (int)((Math.log10(lastSol.n) - logMinN) / logNRange * width);
            int lastY = margin + height - (int)((Math.log10(lastSol.error) - logMinError) / logRange * height);
            g2.fillOval(lastX - 3, lastY - 3, 6, 6);
            
            // Titre
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            Solution firstSol = solutions.get(0);
            g2.drawString("Évolution de l'erreur - " + firstSol.type.getDescription(), margin, 25);
            g2.drawString("f = " + firstSol.function.getName(), margin, 45);
            
            // Labels des axes
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.drawString("Nombre de mailles (log)", width/2, height + margin + 40);
            
            // Rotation pour le label Y
            AffineTransform orig = g2.getTransform();
            g2.rotate(-Math.PI/2);
            g2.drawString("Erreur L2 (log)", -height/2 - 50, 20);
            g2.setTransform(orig);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            // Choix de la fonction et du type d'équation
            Function[] functions = {new SinFunction(), new CubicFunction()};
            EquationType[] types = {EquationType.TYPE1, EquationType.TYPE2, EquationType.TYPE3};
            
            JFrame mainFrame = new JFrame("Résolveur d'équations différentielles");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout());
            
            JPanel controlPanel = new JPanel();
            JComboBox<Function> functionCombo = new JComboBox<>(functions);
            JComboBox<EquationType> typeCombo = new JComboBox<>(types);
            JButton solveButton = new JButton("Résoudre");
            JButton errorButton = new JButton("Courbe d'erreur");
            
            functionCombo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof Function) {
                        setText(((Function) value).getName());
                    }
                    return this;
                }
            });
            
            controlPanel.add(new JLabel("Fonction:"));
            controlPanel.add(functionCombo);
            controlPanel.add(new JLabel("Équation:"));
            controlPanel.add(typeCombo);
            controlPanel.add(solveButton);
            controlPanel.add(errorButton);
            
            JTabbedPane tabbedPane = new JTabbedPane();
            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(tabbedPane, BorderLayout.CENTER);
            
            solveButton.addActionListener(e -> {
                Function selectedFunction = (Function) functionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                tabbedPane.removeAll();
                
                // Résolution pour différents nombres de mailles
                int[] meshSizes = {50, 500, 5000};
                List<Solution> solutions = new ArrayList<>();
                
                System.out.println("\n=== Résultats pour " + selectedType.getDescription() + " avec f = " + selectedFunction.getName() + " ===");
                
                for (int n : meshSizes) {
                    Solution sol = solve(n, selectedType, selectedFunction);
                    solutions.add(sol);
                    
                    System.out.printf("n = %d: Erreur L2 = %.6e\n", n, sol.error);
                    
                    // Calcul de l'ordre de convergence
                    if (solutions.size() > 1) {
                        Solution prevSol = solutions.get(solutions.size() - 2);
                        double order = calculateConvergenceOrder(prevSol.error, sol.error, prevSol.n, sol.n);
                        System.out.printf("Ordre de convergence: %.2f\n", order);
                    }
                    
                    // Ajouter un onglet pour chaque solution
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false);
                    tabbedPane.addTab("n = " + n, panel);
                }
                
                System.out.println();
            });
            
            errorButton.addActionListener(e -> {
                Function selectedFunction = (Function) functionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                // Calcul pour plusieurs tailles de mailles pour la courbe d'erreur
                int[] meshSizes = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000};
                List<Solution> errorSolutions = new ArrayList<>();
                
                for (int n : meshSizes) {
                    Solution sol = solve(n, selectedType, selectedFunction);
                    errorSolutions.add(sol);
                }
                
                GraphPanel errorPanel = new GraphPanel(errorSolutions, true);
                
                JFrame errorFrame = new JFrame("Évolution de l'erreur");
                errorFrame.add(errorPanel);
                errorFrame.setSize(900, 700);
                errorFrame.setLocationRelativeTo(mainFrame);
                errorFrame.setVisible(true);
            });
            
            mainFrame.setSize(1000, 800);
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            // Exemple initial
            functionCombo.setSelectedIndex(0);
            typeCombo.setSelectedIndex(0);
            solveButton.doClick();
        });
    }
}