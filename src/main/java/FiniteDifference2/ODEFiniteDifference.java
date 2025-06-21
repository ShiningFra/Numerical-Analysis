package FiniteDifference2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;

public class ODEFiniteDifference extends JFrame {
    private static final int PANEL_WIDTH = 800;
    private static final int PANEL_HEIGHT = 600;
    private static final double PI = Math.PI;
    
    // Domaine de calcul [0,1] x [0,1]
    private static final double X_MIN = 0.0, X_MAX = 1.0;
    private static final double Y_MIN = 0.0, Y_MAX = 1.0;
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            new ODEFiniteDifference().runSimulations();
        });
    }
    
    public void runSimulations() {
        System.out.println("=== Résolution d'EDP par différences finies 2D ===\n");
        
        // Test avec différents nombres de mailles
        int[] meshSizes = {50, 100, 200, 500, 1000, 2000, 5000};
        
        // Résoudre les trois équations
        EDPSolver solver1 = new EDPSolver1(); // -u'' + u = f
        EDPSolver solver2 = new EDPSolver2(); // -u'' + u' = f  
        EDPSolver solver3 = new EDPSolver3(); // -u'' = f
        
        System.out.println("1. Équation: -u'' + u = f");
        runConvergenceAnalysis(solver1, meshSizes);
        
        System.out.println("\n2. Équation: -u'' + u' = f");
        runConvergenceAnalysis(solver2, meshSizes);
        
        System.out.println("\n3. Équation: -u'' = f");
        runConvergenceAnalysis(solver3, meshSizes);
        
        // Visualisation pour N=50, 500, 5000
        int[] visualMeshSizes = {50, 500, 5000};
        
        for (int N : visualMeshSizes) {
            System.out.println("\n=== Visualisation pour N = " + N + " ===");
            
            // Résoudre avec chaque équation
            Result result1 = solver1.solve(N);
            Result result2 = solver2.solve(N);
            Result result3 = solver3.solve(N);
            
            // Afficher les graphiques
            showVisualization("Equation 1: -u'' + u = f (N=" + N + ")", result1, N);
            showVisualization("Equation 2: -u'' + u' = f (N=" + N + ")", result2, N);
            showVisualization("Equation 3: -u'' = f (N=" + N + ")", result3, N);
        }
        
        // Courbes de convergence
        showConvergencePlots(solver1, solver2, solver3, meshSizes);
    }
    
    private void runConvergenceAnalysis(EDPSolver solver, int[] meshSizes) {
        List<Double> errors = new ArrayList<>();
        
        for (int N : meshSizes) {
            Result result = solver.solve(N);
            errors.add(result.error);
            System.out.printf("N=%d, Erreur=%.6e\n", N, result.error);
        }
        
        // Calcul de l'ordre de convergence
        System.out.println("Ordres de convergence:");
        for (int i = 1; i < errors.size(); i++) {
            double h1 = 1.0 / meshSizes[i-1];
            double h2 = 1.0 / meshSizes[i];
            double order = Math.log(errors.get(i-1) / errors.get(i)) / Math.log(h1 / h2);
            System.out.printf("Entre N=%d et N=%d: ordre = %.3f\n", 
                meshSizes[i-1], meshSizes[i], order);
        }
    }
    
    private void showVisualization(String title, Result result, int N) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(PANEL_WIDTH + 100, PANEL_HEIGHT + 100);
        
        SolutionPanel panel = new SolutionPanel(result, N);
        frame.add(panel);
        frame.setVisible(true);
    }
    
    private void showConvergencePlots(EDPSolver solver1, EDPSolver solver2, EDPSolver solver3, int[] meshSizes) {
        JFrame frame = new JFrame("Courbes de convergence");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(PANEL_WIDTH, PANEL_HEIGHT);
        
        ConvergencePlotPanel panel = new ConvergencePlotPanel(solver1, solver2, solver3, meshSizes);
        frame.add(panel);
        frame.setVisible(true);
    }
    
    // Fonction f(x,y) = sin(pi*x) + y^3
    public static double f(double x, double y) {
        return Math.sin(PI * x) + y * y * y;
    }
    
    // Solution exacte pour -u'' + u = f (approximation)
    public static double exactSolution1(double x, double y) {
        return Math.sin(PI * x) / (PI * PI + 1) + y * y * y;
    }
    
    // Solution exacte pour -u'' + u' = f (approximation)
    public static double exactSolution2(double x, double y) {
        return Math.sin(PI * x) / (PI * PI) + y * y * y / 2;
    }
    
    // Solution exacte pour -u'' = f (approximation)
    public static double exactSolution3(double x, double y) {
        return -Math.sin(PI * x) / (PI * PI) + y * y * y * y / 12;
    }
    
    // Classe abstraite pour les solveurs
    abstract static class EDPSolver {
        public abstract Result solve(int N);
        protected abstract double exactSolution(double x, double y);
        
        protected double[][] buildMatrix(int N) {
            int size = (N-1) * (N-1);
            double[][] A = new double[size][size];
            double h = 1.0 / N;
            
            for (int i = 0; i < N-1; i++) {
                for (int j = 0; j < N-1; j++) {
                    int idx = i * (N-1) + j;
                    setupMatrixRow(A, idx, i, j, N, h);
                }
            }
            return A;
        }
        
        protected abstract void setupMatrixRow(double[][] A, int idx, int i, int j, int N, double h);
        
        protected double[] buildRHS(int N) {
            double[] b = new double[(N-1) * (N-1)];
            double h = 1.0 / N;
            
            for (int i = 0; i < N-1; i++) {
                for (int j = 0; j < N-1; j++) {
                    double x = (i + 1) * h;
                    double y = (j + 1) * h;
                    b[i * (N-1) + j] = f(x, y);
                }
            }
            return b;
        }
        
        protected double[] solveLinearSystem(double[][] A, double[] b) {
            return gaussianElimination(A, b);
        }
        
        private double[] gaussianElimination(double[][] A, double[] b) {
            int n = A.length;
            double[][] augmented = new double[n][n + 1];
            
            // Créer la matrice augmentée
            for (int i = 0; i < n; i++) {
                System.arraycopy(A[i], 0, augmented[i], 0, n);
                augmented[i][n] = b[i];
            }
            
            // Élimination de Gauss
            for (int i = 0; i < n; i++) {
                // Pivot partiel
                int maxRow = i;
                for (int k = i + 1; k < n; k++) {
                    if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
                        maxRow = k;
                    }
                }
                
                // Échanger les lignes
                double[] temp = augmented[i];
                augmented[i] = augmented[maxRow];
                augmented[maxRow] = temp;
                
                // Élimination
                for (int k = i + 1; k < n; k++) {
                    double factor = augmented[k][i] / augmented[i][i];
                    for (int j = i; j <= n; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
            
            // Substitution arrière
            double[] x = new double[n];
            for (int i = n - 1; i >= 0; i--) {
                x[i] = augmented[i][n];
                for (int j = i + 1; j < n; j++) {
                    x[i] -= augmented[i][j] * x[j];
                }
                x[i] /= augmented[i][i];
            }
            
            return x;
        }
        
        protected double calculateError(double[] solution, int N) {
            double error = 0.0;
            double h = 1.0 / N;
            
            for (int i = 0; i < N-1; i++) {
                for (int j = 0; j < N-1; j++) {
                    double x = (i + 1) * h;
                    double y = (j + 1) * h;
                    double numerical = solution[i * (N-1) + j];
                    double exact = exactSolution(x, y);
                    error += (numerical - exact) * (numerical - exact);
                }
            }
            return Math.sqrt(error * h * h);
        }
    }
    
    // Solveur pour -u'' + u = f
    static class EDPSolver1 extends EDPSolver {
        @Override
        protected double exactSolution(double x, double y) {
            return exactSolution1(x, y);
        }
        
        @Override
        protected void setupMatrixRow(double[][] A, int idx, int i, int j, int N, double h) {
            double coeff = 1.0 / (h * h);
            
            // Point central
            A[idx][idx] = -4 * coeff + 1.0;
            
            // Voisins
            if (i > 0) A[idx][idx - (N-1)] = coeff;
            if (i < N-2) A[idx][idx + (N-1)] = coeff;
            if (j > 0) A[idx][idx - 1] = coeff;
            if (j < N-2) A[idx][idx + 1] = coeff;
        }
        
        @Override
        public Result solve(int N) {
            double[][] A = buildMatrix(N);
            double[] b = buildRHS(N);
            double[] solution = solveLinearSystem(A, b);
            double error = calculateError(solution, N);
            
            return new Result(solution, error, N);
        }
    }
    
    // Solveur pour -u'' + u' = f
    static class EDPSolver2 extends EDPSolver {
        @Override
        protected double exactSolution(double x, double y) {
            return exactSolution2(x, y);
        }
        
        @Override
        protected void setupMatrixRow(double[][] A, int idx, int i, int j, int N, double h) {
            double coeff = 1.0 / (h * h);
            double derivCoeff = 1.0 / (2 * h);
            
            // Point central
            A[idx][idx] = -4 * coeff;
            
            // Voisins pour le Laplacien
            if (i > 0) A[idx][idx - (N-1)] = coeff;
            if (i < N-2) A[idx][idx + (N-1)] = coeff;
            if (j > 0) A[idx][idx - 1] = coeff;
            if (j < N-2) A[idx][idx + 1] = coeff;
            
            // Terme dérivé première (u')
            if (i < N-2) A[idx][idx + (N-1)] += derivCoeff;
            if (i > 0) A[idx][idx - (N-1)] -= derivCoeff;
        }
        
        @Override
        public Result solve(int N) {
            double[][] A = buildMatrix(N);
            double[] b = buildRHS(N);
            double[] solution = solveLinearSystem(A, b);
            double error = calculateError(solution, N);
            
            return new Result(solution, error, N);
        }
    }
    
    // Solveur pour -u'' = f
    static class EDPSolver3 extends EDPSolver {
        @Override
        protected double exactSolution(double x, double y) {
            return exactSolution3(x, y);
        }
        
        @Override
        protected void setupMatrixRow(double[][] A, int idx, int i, int j, int N, double h) {
            double coeff = 1.0 / (h * h);
            
            // Point central
            A[idx][idx] = -4 * coeff;
            
            // Voisins
            if (i > 0) A[idx][idx - (N-1)] = coeff;
            if (i < N-2) A[idx][idx + (N-1)] = coeff;
            if (j > 0) A[idx][idx - 1] = coeff;
            if (j < N-2) A[idx][idx + 1] = coeff;
        }
        
        @Override
        public Result solve(int N) {
            double[][] A = buildMatrix(N);
            double[] b = buildRHS(N);
            double[] solution = solveLinearSystem(A, b);
            double error = calculateError(solution, N);
            
            return new Result(solution, error, N);
        }
    }
    
    // Classe pour stocker les résultats
    static class Result {
        double[] solution;
        double error;
        int N;
        
        Result(double[] solution, double error, int N) {
            this.solution = solution;
            this.error = error;
            this.N = N;
        }
    }
    
    // Panel pour afficher la solution
    class SolutionPanel extends JPanel implements MouseMotionListener {
        private Result result;
        private int N;
        private double minVal, maxVal;
        private String mouseInfo = "";
        
        public SolutionPanel(Result result, int N) {
            this.result = result;
            this.N = N;
            this.addMouseMotionListener(this);
            
            // Calculer min/max pour la colormap
            minVal = Double.MAX_VALUE;
            maxVal = Double.MIN_VALUE;
            for (double val : result.solution) {
                minVal = Math.min(minVal, val);
                maxVal = Math.max(maxVal, val);
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int startX = 50;
            int startY = 50;
            
            // Dessiner la grille colorée
            for (int i = 0; i < N-1; i++) {
                for (int j = 0; j < N-1; j++) {
                    double val = result.solution[i * (N-1) + j];
                    Color color = getColor(val);
                    
                    int x = startX + (i * width) / (N-1);
                    int y = startY + (j * height) / (N-1);
                    int w = width / (N-1) + 1;
                    int h = height / (N-1) + 1;
                    
                    g2d.setColor(color);
                    g2d.fillRect(x, y, w, h);
                }
            }
            
            // Dessiner les axes et labels
            g2d.setColor(Color.BLACK);
            g2d.drawRect(startX, startY, width, height);
            
            // Labels
            g2d.setFont(new Font("Arial", Font.BOLD, 12));
            g2d.drawString("x", startX + width + 5, startY + height/2);
            g2d.drawString("y", startX + width/2, startY - 5);
            
            // Colorbar
            drawColorbar(g2d, startX + width + 20, startY, 20, height);
            
            // Informations
            g2d.drawString(String.format("Erreur: %.6e", result.error), 10, 20);
            g2d.drawString(mouseInfo, 10, getHeight() - 20);
        }
        
        private Color getColor(double val) {
            double normalized = (val - minVal) / (maxVal - minVal);
            normalized = Math.max(0, Math.min(1, normalized));
            
            // Colormap bleu -> rouge
            float r = (float) normalized;
            float g = 0.5f;
            float b = (float) (1 - normalized);
            
            return new Color(r, g, b);
        }
        
        private void drawColorbar(Graphics2D g2d, int x, int y, int w, int h) {
            for (int i = 0; i < h; i++) {
                double val = maxVal - (maxVal - minVal) * i / h;
                Color color = getColor(val);
                g2d.setColor(color);
                g2d.fillRect(x, y + i, w, 1);
            }
            
            g2d.setColor(Color.BLACK);
            g2d.drawRect(x, y, w, h);
            
            // Labels de la colorbar
            g2d.setFont(new Font("Arial", Font.PLAIN, 10));
            g2d.drawString(String.format("%.3f", maxVal), x + w + 5, y + 10);
            g2d.drawString(String.format("%.3f", (maxVal + minVal)/2), x + w + 5, y + h/2);
            g2d.drawString(String.format("%.3f", minVal), x + w + 5, y + h - 5);
        }
        
        @Override
        public void mouseMoved(MouseEvent e) {
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int startX = 50;
            int startY = 50;
            
            if (e.getX() >= startX && e.getX() <= startX + width && 
                e.getY() >= startY && e.getY() <= startY + height) {
                
                int i = ((e.getX() - startX) * (N-1)) / width;
                int j = ((e.getY() - startY) * (N-1)) / height;
                
                if (i >= 0 && i < N-1 && j >= 0 && j < N-1) {
                    double val = result.solution[i * (N-1) + j];
                    double x = (double)(i + 1) / N;
                    double y = (double)(j + 1) / N;
                    mouseInfo = String.format("Position: (%.3f, %.3f), Valeur: %.6f", x, y, val);
                    repaint();
                }
            }
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            mouseMoved(e);
        }
    }
    
    // Panel pour les courbes de convergence
    class ConvergencePlotPanel extends JPanel {
        private EDPSolver solver1, solver2, solver3;
        private int[] meshSizes;
        private List<Double> errors1, errors2, errors3;
        
        public ConvergencePlotPanel(EDPSolver solver1, EDPSolver solver2, EDPSolver solver3, int[] meshSizes) {
            this.solver1 = solver1;
            this.solver2 = solver2;
            this.solver3 = solver3;
            this.meshSizes = meshSizes;
            
            // Calculer les erreurs
            errors1 = new ArrayList<>();
            errors2 = new ArrayList<>();
            errors3 = new ArrayList<>();
            
            for (int N : meshSizes) {
                errors1.add(solver1.solve(N).error);
                errors2.add(solver2.solve(N).error);
                errors3.add(solver3.solve(N).error);
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int startX = 50;
            int startY = 50;
            
            // Axes
            g2d.setColor(Color.BLACK);
            g2d.drawLine(startX, startY + height, startX + width, startY + height);
            g2d.drawLine(startX, startY, startX, startY + height);
            
            // Échelles logarithmiques
            double minH = 1.0 / meshSizes[meshSizes.length - 1];
            double maxH = 1.0 / meshSizes[0];
            double minError = Math.min(Math.min(getMinError(errors1), getMinError(errors2)), getMinError(errors3));
            double maxError = Math.max(Math.max(getMaxError(errors1), getMaxError(errors2)), getMaxError(errors3));
            
            // Dessiner les courbes
            drawCurve(g2d, errors1, Color.RED, "Eq1: -u'' + u = f", startX, startY, width, height, minH, maxH, minError, maxError);
            drawCurve(g2d, errors2, Color.BLUE, "Eq2: -u'' + u' = f", startX, startY, width, height, minH, maxH, minError, maxError);
            drawCurve(g2d, errors3, Color.GREEN, "Eq3: -u'' = f", startX, startY, width, height, minH, maxH, minError, maxError);
            
            // Légende
            g2d.setColor(Color.RED);
            g2d.drawString("Eq1: -u'' + u = f", startX + 10, startY + 20);
            g2d.setColor(Color.BLUE);
            g2d.drawString("Eq2: -u'' + u' = f", startX + 10, startY + 40);
            g2d.setColor(Color.GREEN);
            g2d.drawString("Eq3: -u'' = f", startX + 10, startY + 60);
            
            // Labels des axes
            g2d.setColor(Color.BLACK);
            g2d.drawString("h (pas de maillage)", startX + width/2, startY + height + 30);
            
            // Rotation pour l'axe Y
            g2d.rotate(-Math.PI/2);
            g2d.drawString("Erreur L2", -(startY + height/2), 20);
            g2d.rotate(Math.PI/2);
        }
        
        private void drawCurve(Graphics2D g2d, List<Double> errors, Color color, String label, 
                              int startX, int startY, int width, int height, 
                              double minH, double maxH, double minError, double maxError) {
            g2d.setColor(color);
            
            for (int i = 0; i < errors.size() - 1; i++) {
                double h1 = 1.0 / meshSizes[i];
                double h2 = 1.0 / meshSizes[i + 1];
                double error1 = errors.get(i);
                double error2 = errors.get(i + 1);
                
                int x1 = startX + (int) ((Math.log(h1) - Math.log(minH)) * width / (Math.log(maxH) - Math.log(minH)));
                int y1 = startY + height - (int) ((Math.log(error1) - Math.log(minError)) * height / (Math.log(maxError) - Math.log(minError)));
                int x2 = startX + (int) ((Math.log(h2) - Math.log(minH)) * width / (Math.log(maxH) - Math.log(minH)));
                int y2 = startY + height - (int) ((Math.log(error2) - Math.log(minError)) * height / (Math.log(maxError) - Math.log(minError)));
                
                g2d.drawLine(x1, y1, x2, y2);
                g2d.fillOval(x1 - 2, y1 - 2, 4, 4);
            }
        }
        
        private double getMinError(List<Double> errors) {
            return errors.stream().mapToDouble(Double::doubleValue).min().orElse(1e-10);
        }
        
        private double getMaxError(List<Double> errors) {
            return errors.stream().mapToDouble(Double::doubleValue).max().orElse(1.0);
        }
    }
}