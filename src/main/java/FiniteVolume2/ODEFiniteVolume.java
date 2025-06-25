package FiniteVolume2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

public class ODEFiniteVolume extends JFrame {
    
    // Paramètres du problème
    private static final double PI = Math.PI;
    private static final double DOMAIN_SIZE = 1.0; // Domaine [0,1] x [0,1]
    
    // Données pour les graphiques
    private List<Integer> meshSizes = new ArrayList<>();
    private List<Double> errors = new ArrayList<>();
    private List<Double> convergenceOrders = new ArrayList<>();
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            ODEFiniteVolume solver = new ODEFiniteVolume();
            solver.solve();
        });
    }
    
    public void solve() {
        System.out.println("=== Résolution de -u'' = f par Volumes Finis 2D ===");
        System.out.println("f(x,y) = sin(π*x) * sin(π*y)");
        System.out.println("Solution exacte: u(x,y) = (1/(2*π²)) * sin(π*x) * sin(π*y)");
        System.out.println();
        
        int[] nValues = {10, 20, 40};
        
        for (int n : nValues) {
            System.out.println("--- Résolution avec " + n + "x" + n + " mailles ---");
            
            // Résolution numérique
            double[][] numericalSolution = solveFiniteVolume(n);
            
            // Calcul de l'erreur
            double error = calculateError(numericalSolution, n);
            
            // Stockage des résultats
            meshSizes.add(n);
            errors.add(error);
            
            System.out.printf("Erreur L2: %.6e%n", error);
            
            // Calcul de l'ordre de convergence
            if (meshSizes.size() > 1) {
                double order = calculateConvergenceOrder();
                convergenceOrders.add(order);
                System.out.printf("Ordre de convergence: %.3f%n", order);
            }
            
            System.out.println();
        }
        
        // Affichage des graphiques
        displayResults(nValues);
    }
    
    // Fonction source f(x,y)
    private double sourceFunction(double x, double y) {
        return Math.sin(PI * x) * Math.sin(PI * y);
    }
    
    // Solution exacte u(x,y)
    private double exactSolution(double x, double y) {
        return (1.0 / (2.0 * PI * PI)) * Math.sin(PI * x) * Math.sin(PI * y);
    }
    
    // Résolution par méthode des volumes finis
    private double[][] solveFiniteVolume(int n) {
        double h = DOMAIN_SIZE / n;
        int totalNodes = (n + 1) * (n + 1);
        
        // Matrice du système et vecteur second membre
        double[][] A = new double[totalNodes][totalNodes];
        double[] b = new double[totalNodes];
        
        // Construction du système linéaire
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                int index = i * (n + 1) + j;
                
                // Conditions aux limites (Dirichlet homogène)
                if (i == 0 || i == n || j == 0 || j == n) {
                    A[index][index] = 1.0;
                    b[index] = 0.0;
                } else {
                    // Points intérieurs - schéma volumes finis
                    double x = i * h;
                    double y = j * h;
                    
                    // Coefficient central
                    A[index][index] = 4.0 / (h * h);
                    
                    // Coefficients voisins
                    A[index][index - 1] = -1.0 / (h * h); // gauche
                    A[index][index + 1] = -1.0 / (h * h); // droite
                    A[index][index - (n + 1)] = -1.0 / (h * h); // bas
                    A[index][index + (n + 1)] = -1.0 / (h * h); // haut
                    
                    // Terme source
                    b[index] = sourceFunction(x, y);
                }
            }
        }
        
        // Résolution du système linéaire par méthode directe (Gauss)
        double[] solution = solveLinearSystem(A, b);
        
        // Conversion en matrice 2D
        double[][] u = new double[n + 1][n + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                u[i][j] = solution[i * (n + 1) + j];
            }
        }
        
        return u;
    }
    
    // Résolution système linéaire par élimination de Gauss
    private double[] solveLinearSystem(double[][] A, double[] b) {
        int n = A.length;
        double[][] augmented = new double[n][n + 1];
        
        // Matrice augmentée
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, augmented[i], 0, n);
            augmented[i][n] = b[i];
        }
        
        // Élimination de Gauss avec pivotage partiel
        for (int k = 0; k < n - 1; k++) {
            // Recherche du pivot
            int maxRow = k;
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(augmented[i][k]) > Math.abs(augmented[maxRow][k])) {
                    maxRow = i;
                }
            }
            
            // Échange de lignes
            if (maxRow != k) {
                double[] temp = augmented[k];
                augmented[k] = augmented[maxRow];
                augmented[maxRow] = temp;
            }
            
            // Élimination
            for (int i = k + 1; i < n; i++) {
                if (Math.abs(augmented[k][k]) > 1e-14) {
                    double factor = augmented[i][k] / augmented[k][k];
                    for (int j = k; j <= n; j++) {
                        augmented[i][j] -= factor * augmented[k][j];
                    }
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
            if (Math.abs(augmented[i][i]) > 1e-14) {
                x[i] /= augmented[i][i];
            }
        }
        
        return x;
    }
    
    // Calcul de l'erreur L2
    private double calculateError(double[][] numerical, int n) {
        double h = DOMAIN_SIZE / n;
        double error = 0.0;
        
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= n; j++) {
                double x = i * h;
                double y = j * h;
                double exact = exactSolution(x, y);
                double diff = numerical[i][j] - exact;
                error += diff * diff;
            }
        }
        
        return Math.sqrt(error);
    }
    
    // Calcul de l'ordre de convergence
    private double calculateConvergenceOrder() {
        int size = errors.size();
        if (size < 2) return 0.0;
        
        double e1 = errors.get(size - 2);
        double e2 = errors.get(size - 1);
        double h1 = DOMAIN_SIZE / meshSizes.get(size - 2);
        double h2 = DOMAIN_SIZE / meshSizes.get(size - 1);
        
        return Math.sqrt(Math.log(e2) / (2 * Math.log(h2)) * Math.log(e1) / (2 * Math.log(h1)))/* / Math.log(h2 / h1)*/;
    }
    
    // Affichage des résultats graphiques
    private void displayResults(int[] nValues) {
        setTitle("Résolution EDP par Volumes Finis 2D");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setExtendedState(JFrame.MAXIMIZED_BOTH);
        
        JTabbedPane tabbedPane = new JTabbedPane();
        
        // Onglets pour les solutions
        for (int n : nValues) {
            double[][] solution = solveFiniteVolume(n);
            SolutionPanel panel = new SolutionPanel(solution, n, "Solution numérique " + n + "x" + n);
            tabbedPane.addTab(n + "x" + n, panel);
        }
        
        // Onglet pour la solution exacte
        double[][] exactSol = new double[81][81];
        for (int i = 0; i < 81; i++) {
            for (int j = 0; j < 81; j++) {
                double x = i / 80.0;
                double y = j / 80.0;
                exactSol[i][j] = exactSolution(x, y);
            }
        }
        SolutionPanel exactPanel = new SolutionPanel(exactSol, 80, "Solution exacte");
        tabbedPane.addTab("Exacte", exactPanel);
        
        // Onglet pour la courbe d'erreur
        ErrorPanel errorPanel = new ErrorPanel();
        tabbedPane.addTab("Convergence", errorPanel);
        
        add(tabbedPane);
        setVisible(true);
    }
    
    // Panel pour afficher une solution avec carte de couleurs
    class SolutionPanel extends JPanel implements MouseMotionListener {
        private double[][] data;
        private int n;
        private String title;
        private double minVal, maxVal;
        private Point mousePos = new Point();
        
        public SolutionPanel(double[][] data, int n, String title) {
            this.data = data;
            this.n = n;
            this.title = title;
            addMouseMotionListener(this);
            
            // Calcul min/max
            minVal = Double.MAX_VALUE;
            maxVal = Double.MIN_VALUE;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    minVal = Math.min(minVal, data[i][j]);
                    maxVal = Math.max(maxVal, data[i][j]);
                }
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g.create();
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int startX = 50;
            int startY = 50;
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 16));
            g2d.drawString(title, startX, 30);
            
            // Carte de couleurs
            double cellWidth = (double) width / data.length;
            double cellHeight = (double) height / data[0].length;
            
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    double normalized = (data[i][j] - minVal) / (maxVal - minVal);
                    Color color = getColorForValue(normalized);
                    g2d.setColor(color);
                    
                    int x = startX + (int) (i * cellWidth);
                    int y = startY + (int) ((data[i].length - 1 - j) * cellHeight);
                    int w = (int) Math.ceil(cellWidth);
                    int h = (int) Math.ceil(cellHeight);
                    
                    g2d.fillRect(x, y, w, h);
                }
            }
            
            // Barre de couleurs
            drawColorBar(g2d, startX + width + 20, startY, 30, height);
            
            // Valeur au curseur
            if (mousePos.x > startX && mousePos.x < startX + width &&
                mousePos.y > startY && mousePos.y < startY + height) {
                
                int i = (int) ((mousePos.x - startX) / cellWidth);
                int j = data[0].length - 1 - (int) ((mousePos.y - startY) / cellHeight);
                
                if (i >= 0 && i < data.length && j >= 0 && j < data[0].length) {
                    double x = (double) i / (data.length - 1);
                    double y = (double) j / (data[0].length - 1);
                    
                    g2d.setColor(Color.BLACK);
                    g2d.setFont(new Font("Arial", Font.PLAIN, 12));
                    String info = String.format("x=%.3f, y=%.3f, u=%.6f", x, y, data[i][j]);
                    g2d.drawString(info, mousePos.x + 10, mousePos.y - 10);
                }
            }
            
            g2d.dispose();
        }
        
        private Color getColorForValue(double normalized) {
            // Palette bleu -> cyan -> vert -> jaune -> rouge
            if (normalized < 0.25) {
                float t = (float) (normalized * 4);
                return new Color(0, (int) (255 * t), 255);
            } else if (normalized < 0.5) {
                float t = (float) ((normalized - 0.25) * 4);
                return new Color(0, 255, (int) (255 * (1 - t)));
            } else if (normalized < 0.75) {
                float t = (float) ((normalized - 0.5) * 4);
                return new Color((int) (255 * t), 255, 0);
            } else {
                float t = (float) ((normalized - 0.75) * 4);
                return new Color(255, (int) (255 * (1 - t)), 0);
            }
        }
        
        private void drawColorBar(Graphics2D g2d, int x, int y, int width, int height) {
            int steps = 100;
            double stepHeight = (double) height / steps;
            
            for (int i = 0; i < steps; i++) {
                double normalized = (double) i / (steps - 1);
                Color color = getColorForValue(normalized);
                g2d.setColor(color);
                g2d.fillRect(x, y + (int) ((steps - 1 - i) * stepHeight), width, (int) Math.ceil(stepHeight));
            }
            
            // Bordure
            g2d.setColor(Color.BLACK);
            g2d.drawRect(x, y, width, height);
            
            // Étiquettes
            g2d.setFont(new Font("Arial", Font.PLAIN, 10));
            g2d.drawString(String.format("%.3e", maxVal), x + width + 5, y + 5);
            g2d.drawString(String.format("%.3e", minVal), x + width + 5, y + height);
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            mousePos = e.getPoint();
            repaint();
        }
        
        @Override
        public void mouseMoved(MouseEvent e) {
            mousePos = e.getPoint();
            repaint();
        }
    }
    
    // Panel pour la courbe d'erreur et convergence
    class ErrorPanel extends JPanel {
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g.create();
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int startX = 50;
            int startY = 50;
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 16));
            g2d.drawString("Évolution de l'erreur L2 et ordre de convergence", startX, 30);
            
            if (errors.size() < 2) return;
            
            // Graphique log-log de l'erreur
            drawErrorGraph(g2d, startX, startY, width / 2 - 20, height);
            
            // Graphique de l'ordre de convergence
            drawConvergenceGraph(g2d, startX + width / 2 + 20, startY, width / 2 - 20, height);
            
            g2d.dispose();
        }
        
        private void drawErrorGraph(Graphics2D g2d, int x, int y, int width, int height) {
            // Axes
            g2d.setColor(Color.BLACK);
            g2d.drawLine(x, y + height, x + width, y + height); // axe x
            g2d.drawLine(x, y, x, y + height); // axe y
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 12));
            g2d.drawString("Erreur L2 vs Taille de maille", x, y - 10);
            
            // Calcul des échelles logarithmiques
            double minH = DOMAIN_SIZE / meshSizes.get(meshSizes.size() - 1);
            double maxH = DOMAIN_SIZE / meshSizes.get(0);
            double minE = errors.stream().min(Double::compare).orElse(1e-10);
            double maxE = errors.stream().max(Double::compare).orElse(1e-5);
            
            // Points de données
            g2d.setColor(Color.RED);
            for (int i = 0; i < meshSizes.size(); i++) {
                double h = DOMAIN_SIZE / meshSizes.get(i);
                double error = errors.get(i);
                
                int px = x + (int) (width * (Math.log(h) - Math.log(minH)) / (Math.log(maxH) - Math.log(minH)));
                int py = y + height - (int) (height * (Math.log(error) - Math.log(minE)) / (Math.log(maxE) - Math.log(minE)));
                
                g2d.fillOval(px - 3, py - 3, 6, 6);
                
                if (i > 0) {
                    double h_prev = DOMAIN_SIZE / meshSizes.get(i - 1);
                    double error_prev = errors.get(i - 1);
                    int px_prev = x + (int) (width * (Math.log(h_prev) - Math.log(minH)) / (Math.log(maxH) - Math.log(minH)));
                    int py_prev = y + height - (int) (height * (Math.log(error_prev) - Math.log(minE)) / (Math.log(maxE) - Math.log(minE)));
                    g2d.drawLine(px_prev, py_prev, px, py);
                }
            }
            
            // Ligne théorique de pente 1
            g2d.setColor(Color.BLUE);
            g2d.setStroke(new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{5}, 0));
            double refError = errors.get(0);
            double refH = DOMAIN_SIZE / meshSizes.get(0);
            
            for (int i = 0; i < meshSizes.size(); i++) {
                double h = DOMAIN_SIZE / meshSizes.get(i);
                double theoreticalError = refError * (h / refH); // Pente 1
                
                int px = x + (int) (width * (Math.log(h) - Math.log(minH)) / (Math.log(maxH) - Math.log(minH)));
                int py = y + height - (int) (height * (Math.log(theoreticalError) - Math.log(minE)) / (Math.log(maxE) - Math.log(minE)));
                
                if (i == 0) {
                    g2d.fillOval(px - 2, py - 2, 4, 4);
                } else {
                    double h_prev = DOMAIN_SIZE / meshSizes.get(i - 1);
                    double theoreticalError_prev = refError * (h_prev / refH);
                    int px_prev = x + (int) (width * (Math.log(h_prev) - Math.log(minH)) / (Math.log(maxH) - Math.log(minH)));
                    int py_prev = y + height - (int) (height * (Math.log(theoreticalError_prev) - Math.log(minE)) / (Math.log(maxE) - Math.log(minE)));
                    g2d.drawLine(px_prev, py_prev, px, py);
                }
            }
            
            // Légende
            g2d.setStroke(new BasicStroke(1));
            g2d.setFont(new Font("Arial", Font.PLAIN, 10));
            g2d.setColor(Color.RED);
            g2d.drawString("Erreur calculée", x + 10, y + 20);
            g2d.setColor(Color.BLUE);
            g2d.drawString("Pente théorique = 1", x + 10, y + 35);
        }
        
        private void drawConvergenceGraph(Graphics2D g2d, int x, int y, int width, int height) {
            if (convergenceOrders.isEmpty()) return;
            
            // Axes
            g2d.setColor(Color.BLACK);
            g2d.drawLine(x, y + height, x + width, y + height); // axe x
            g2d.drawLine(x, y, x, y + height); // axe y
            
            // Titre
            g2d.setFont(new Font("Arial", Font.BOLD, 12));
            g2d.drawString("Ordre de convergence", x, y - 10);
            
            // Ligne de référence à y = 1
            g2d.setColor(Color.GRAY);
            g2d.setStroke(new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0, new float[]{3}, 0));
            int refY = y + height - (int) (height * 0.5); // Ordre 1 au milieu
            g2d.drawLine(x, refY, x + width, refY);
            g2d.drawString("Ordre = 1", x + 5, refY - 5);
            
            // Points de convergence
            g2d.setColor(Color.GREEN);
            g2d.setStroke(new BasicStroke(2));
            
            double maxOrder = Math.max(2.0, convergenceOrders.stream().max(Double::compare).orElse(1.5));
            double minOrder = Math.min(0.0, convergenceOrders.stream().min(Double::compare).orElse(0.5));
            
            for (int i = 0; i < convergenceOrders.size(); i++) {
                double order = convergenceOrders.get(i);
                int px = x + (i + 1) * width / (convergenceOrders.size() + 1);
                int py = y + height - (int) (height * (order - minOrder) / (maxOrder - minOrder));
                
                g2d.fillOval(px - 4, py - 4, 8, 8);
                
                // Étiquette
                g2d.setFont(new Font("Arial", Font.PLAIN, 10));
                g2d.drawString(String.format("%.2f", order), px - 10, py - 10);
                
                if (i > 0) {
                    double order_prev = convergenceOrders.get(i - 1);
                    int px_prev = x + i * width / (convergenceOrders.size() + 1);
                    int py_prev = y + height - (int) (height * (order_prev - minOrder) / (maxOrder - minOrder));
                    g2d.drawLine(px_prev, py_prev, px, py);
                }
            }
        }
    }
}