/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;
import java.util.function.DoubleUnaryOperator;

/**
 * FiniteVolumeDeduction calcule et trace la solution exacte et les solutions FV (décentrée et centrée)
 * pour un problème déduit à partir d'une solution exacte donnée.
 * Il intègre également une analyse de convergence (erreur et ordre de convergence)
 * pour différents nombres de volumes internes.
 * 
 * Exemple 1: u(x)=sin(pi*x)  => f(x)=pi^2*sin(pi*x), u(0)=u(1)=0
 * Exemple 2: u(x)=x^3        => f(x)=-6x, u(0)=0, u(1)=1
 * 
 * @author Roddier
 */
public class FiniteVolumeDeduction extends JPanel {

    // Nombre de volumes internes (la grille complète comporte N+2 points, incluant les bornes)
    private int N;
    private double h;  // pas de la grille: h = 1/(N+1)

    // Grille complète [0,1]
    private double[] x; // positions des points : x[0]=0, ..., x[N+1]=1

    // Solutions complètes (de taille N+2) :
    // uExact : solution exacte
    // uDecentered : solution numérique obtenue avec la méthode décentrée
    // uCentered : solution numérique obtenue avec la méthode centrée
    private double[] uExact, uDecentered, uCentered;

    // Choix de l'exemple : 1 pour u(x)=sin(pi*x), 2 pour u(x)=x^3.
    private int example;

    // Constructeur qui prend l'exemple et le nombre de volumes internes.
    public FiniteVolumeDeduction(int example, int N) {
        this.example = example;
        this.N = N;
        this.h = 1.0 / (N + 1);
        initGrid();
        computeSolutions();
    }

    // Constructeur par défaut pour affichage graphique (N fixé par défaut, ici 50)
    public FiniteVolumeDeduction(int example) {
        this(example, 50);
    }

    // Construction de la grille sur [0,1] avec N+2 points.
    private void initGrid() {
        x = new double[N + 2];
        for (int i = 0; i < N + 2; i++) {
            x[i] = i * h;
        }
    }

    // Définir u exacte et en déduire f = -u'' (à utiliser pour construire le second membre)
    // Puis résoudre numériquement par la méthode des volumes finis avec deux variantes.
    private void computeSolutions() {
        // Définition de uExact et f en fonction de l'exemple.
        DoubleUnaryOperator uExactOp;
        DoubleUnaryOperator fOp;
        if (example == 1) {
            // u(x) = sin(pi*x) avec u(0)=0, u(1)=0
            uExactOp = (x) -> Math.sin(Math.PI * x);
            // u''(x) = -pi^2*sin(pi*x) donc f(x) = -u''(x) = pi^2*sin(pi*x)
            fOp = (x) -> Math.PI * Math.PI * Math.sin(Math.PI * x);
        } else {
            // u(x) = x^3 avec u(0)=0, u(1)=1
            uExactOp = (x) -> Math.pow(x, 3);
            // u''(x)= 6x donc f(x) = -6x
            fOp = (x) -> -6 * x;
        }

        // Construction de la solution exacte sur la grille complète
        uExact = new double[N + 2];
        for (int i = 0; i < N + 2; i++) {
            uExact[i] = uExactOp.applyAsDouble(x[i]);
        }

        // Les inconnues numériques seront calculées pour les points internes
        double[] xInternal = new double[N];
        for (int i = 0; i < N; i++) {
            // On prend ici xInternal[i] = x[i+1] (centre de la cellule)
            xInternal[i] = x[i + 1];
        }

        // Construction des seconds membres pour chaque méthode.
        // Pour une discrétisation FV de -u'' = f, le schéma usuel sur un maillage uniforme donne :
        //   (2/h^2) u_i - (1/h^2)(u_{i-1}+u_{i+1}) = f_i   (avec f_i défini par la méthode choisie)
        // Pour tenir compte des conditions aux limites, on ajoute la contribution de u(0) et u(1)
        double[] bDecentered = new double[N];
        double[] bCentered = new double[N];

        // Méthode décentrée : on évalue f au centre de la cellule
        for (int i = 0; i < N; i++) {
            bDecentered[i] = fOp.applyAsDouble(xInternal[i]);
        }

        // Méthode centrée : on évalue f sur les faces et on en fait la moyenne
        for (int i = 0; i < N; i++) {
            double xLeft = xInternal[i] - h / 2;
            double xRight = xInternal[i] + h / 2;
            bCentered[i] = 0.5 * (fOp.applyAsDouble(xLeft) + fOp.applyAsDouble(xRight));
        }

        // Ajout des contributions des conditions aux limites dans le second membre
        // Pour le premier point interne (i=0) : contribution de u[0]
        // Pour le dernier (i=N-1) : contribution de u[N+1]
        bDecentered[0] += (1.0 / (h * h)) * uExact[0];
        bCentered[0] += (1.0 / (h * h)) * uExact[0];
        bDecentered[N - 1] += (1.0 / (h * h)) * uExact[N + 1];
        bCentered[N - 1] += (1.0 / (h * h)) * uExact[N + 1];

        // Construction de la matrice tridiagonale (taille N) associée à l'opérateur FV.
        double diag = 2.0 / (h * h);
        double offDiag = -1.0 / (h * h);

        // Résolution des systèmes linéaires pour chacune des méthodes
        double[] uInternalDecentered = solveTridiagonal(N, offDiag, diag, offDiag, bDecentered);
        double[] uInternalCentered = solveTridiagonal(N, offDiag, diag, offDiag, bCentered);

        // Constitution des vecteurs solution complets (avec conditions aux limites)
        uDecentered = new double[N + 2];
        uCentered = new double[N + 2];
        // Les conditions aux limites proviennent de uExact (déjà imposées)
        uDecentered[0] = uExact[0];
        uCentered[0] = uExact[0];
        uDecentered[N + 1] = uExact[N + 1];
        uCentered[N + 1] = uExact[N + 1];
        for (int i = 0; i < N; i++) {
            uDecentered[i + 1] = uInternalDecentered[i];
            uCentered[i + 1] = uInternalCentered[i];
        }
    }

    // Résolution d'un système tridiagonal à coefficients constants.
    // On considère le système pour i=0,...,n-1 :
    //    a*u[i-1] + b*u[i] + c*u[i+1] = d[i]
    // (avec u[-1] et u[n] connus via les conditions aux limites, et intégrés dans d)
    private double[] solveTridiagonal(int n, double a, double b, double c, double[] d) {
        double[] cp = new double[n]; // coefficients modifiés pour la sur-diagonale
        double[] dp = new double[n]; // second membre modifié
        double[] sol = new double[n];

        // Élimination de Gauss
        cp[0] = c / b;
        dp[0] = d[0] / b;
        for (int i = 1; i < n; i++) {
            double m = b - a * cp[i - 1];
            cp[i] = c / m;
            dp[i] = (d[i] - a * dp[i - 1]) / m;
        }
        // Retour arrière
        sol[n - 1] = dp[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            sol[i] = dp[i] - cp[i] * sol[i + 1];
        }
        return sol;
    }

    // Méthode utilitaire pour tracer une courbe à partir d'un tableau de points.
    private void drawCurve(Graphics2D g2, DoubleUnaryOperator toScreenX,
                           DoubleUnaryOperator toScreenY, double[] xVals, double[] yVals) {
        Path2D path = new Path2D.Double();
        if (xVals.length == 0) return;
        path.moveTo(toScreenX.applyAsDouble(xVals[0]), toScreenY.applyAsDouble(yVals[0]));
        for (int i = 1; i < xVals.length; i++) {
            path.lineTo(toScreenX.applyAsDouble(xVals[i]), toScreenY.applyAsDouble(yVals[i]));
        }
        g2.draw(path);
    }

    // Affichage graphique des trois solutions (exacte, méthode décentrée et centrée)
    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        int width = getWidth();
        int height = getHeight();
        int margin = 40;

        // Dessin des axes
        g2.drawLine(margin, height - margin, width - margin, height - margin); // axe x
        g2.drawLine(margin, margin, margin, height - margin); // axe y

        // Détermination des bornes en y pour l'affichage
        double yMin = Double.POSITIVE_INFINITY, yMax = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < x.length; i++) {
            yMin = Math.min(yMin, Math.min(Math.min(uDecentered[i], uCentered[i]), uExact[i]));
            yMax = Math.max(yMax, Math.max(Math.max(uDecentered[i], uCentered[i]), uExact[i]));
        }
        double yPadding = 0.1 * (yMax - yMin);
        yMin -= yPadding;
        yMax += yPadding;
        final double yMi = yMin, yMa = yMax;

        // Fonctions de conversion des coordonnées physiques en coordonnées écran
        DoubleUnaryOperator toScreenX = (xval) -> margin + (xval - 0) / (1 - 0) * (width - 2 * margin);
        DoubleUnaryOperator toScreenY = (yval) -> height - margin - (yval - yMi) / (yMa - yMi) * (height - 2 * margin);

        // Tracé des courbes :
        // - Exacte : rouge
        // - Méthode décentrée : bleu
        // - Méthode centrée : vert
        Stroke defaultStroke = g2.getStroke();

        g2.setStroke(new BasicStroke(2));
        g2.setColor(Color.RED);
        drawCurve(g2, toScreenX, toScreenY, x, uExact);

        g2.setStroke(new BasicStroke(2));
        g2.setColor(Color.BLUE);
        drawCurve(g2, toScreenX, toScreenY, x, uDecentered);

        g2.setStroke(new BasicStroke(2));
        g2.setColor(Color.GREEN.darker());
        drawCurve(g2, toScreenX, toScreenY, x, uCentered);

        // Légende
        g2.setStroke(defaultStroke);
        g2.setColor(Color.BLACK);
        g2.drawString("Exacte (rouge)", margin + 10, margin + 15);
        g2.drawString("Decentree (bleu)", margin + 10, margin + 30);
        g2.drawString("Centrée (vert)", margin + 10, margin + 45);
    }

    // Création de la fenêtre d'affichage pour une solution donnée.
    private static void createAndShowGUI(int example) {
        String title;
        if (example == 1) {
            title = "FV - u(x)=sin(pi*x) et f(x)=pi^2*sin(pi*x)";
        } else {
            title = "FV - u(x)=x^3 et f(x)=-6x";
        }
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        FiniteVolumeDeduction panel = new FiniteVolumeDeduction(example);
        panel.setPreferredSize(new Dimension(600, 400));
        frame.getContentPane().add(panel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    // Méthode statique pour calculer l'erreur en norme max entre deux solutions (tableaux de même taille)
    public static double computeError(double[] approx, double[] exact) {
        double error = 0.0;
        for (int i = 0; i < approx.length; i++) {
            error = Math.max(error, Math.abs(approx[i] - exact[i]));
        }
        return error;
    }

    /**
     * Analyse de convergence pour un exemple donné.
     * Pour des maillages de 20, 40, 80, 160, 320 volumes internes,
     * la méthode calcule et affiche l'erreur (norme max) et l'ordre de convergence
     * pour la méthode décentrée et centrée.
     */
    public static void runConvergenceAnalysis(int example) {
        int[] meshSizes = {20, 40, 80, 160, 320};
        System.out.println("Analyse de convergence pour l'exemple " + example);
        if (example == 1) {
            System.out.println("Exemple: u(x)=sin(pi*x), f(x)=pi^2*sin(pi*x), u(0)=u(1)=0");
        } else {
            System.out.println("Exemple: u(x)=x^3, f(x)=-6x, u(0)=0, u(1)=1");
        }
        System.out.println("--------------------------------------------------------------------------");
        System.out.println("Mesh\t h\t\t Error (Déc.)\t Order\t Error (Centrée)\t Order");

        double prevErrorDec = 0.0, prevErrorCent = 0.0, prevH = 0.0;
        for (int i = 0; i < meshSizes.length; i++) {
            int N = meshSizes[i];
            double h = 1.0 / (N + 1);

            // Instanciation de la solution FV pour le maillage courant
            FiniteVolumeDeduction fv = new FiniteVolumeDeduction(example, N);
            // On récupère la solution exacte calculée dans l'instance
            double[] uEx = fv.uExact;
            double errDec = computeError(fv.uDecentered, uEx);
            double errCent = computeError(fv.uCentered, uEx);

            String orderDec = (i == 0) ? "-" :
                    String.format("%.2f", Math.log(prevErrorDec / errDec) / Math.log(prevH / h));
            String orderCent = (i == 0) ? "-" :
                    String.format("%.2f", Math.log(prevErrorCent / errCent) / Math.log(prevH / h));

            System.out.printf("%3d\t%.5f\t%.5e\t%s\t%.5e\t%s\n", N, h, errDec, orderDec, errCent, orderCent);

            prevErrorDec = errDec;
            prevErrorCent = errCent;
            prevH = h;
        }
        System.out.println();
    }

    // Méthode main : lance l'analyse de convergence et ouvre les fenêtres graphiques pour les deux exemples.
    public static void main(String[] args) {
        // Lancer l'analyse de convergence pour les deux exemples dans la console.
        runConvergenceAnalysis(1); // u(x)=sin(pi*x)
        runConvergenceAnalysis(2); // u(x)=x^3

        // Lancer les fenêtres graphiques pour visualiser les solutions.
        SwingUtilities.invokeLater(() -> createAndShowGUI(1)); // Exemple 1
        SwingUtilities.invokeLater(() -> createAndShowGUI(2)); // Exemple 2
    }
}
